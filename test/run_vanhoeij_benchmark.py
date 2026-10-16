"""Build/run fixed pre-safety, safety-net and repaired CLPoly/FLINT comparisons.

First run prepare_vanhoeij_benchmark.py. Supply CXX, CPPFLAGS (dependency
include paths), and LDFLAGS (FLINT/GMP library paths and link flags).
Use --build-only while other tests run, then --run-only on an idle machine.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "_build/vanhoeij-benchmark"
VERSIONS = ["presafety", "safetynet", "repaired"]


def build():
    commands = {}
    for version in VERSIONS:
        snapshot = OUT / version
        command = [os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
                   "-DCLPOLY_PROFILE", "-DCLPOLY_PHASE_BENCH", "-DCLPOLY_FLINT_REFERENCE"]
        if version == "repaired":
            command += ["-DCLPOLY_CERTIFICATE_BENCH"]
        command += ["-I" + str(snapshot)] + shlex.split(os.environ.get("CPPFLAGS", ""))
        command += [str(ROOT / "test/bench_vanhoeij_certificate.cc"),
                    str(snapshot / "clpoly/upolynomial.cc"), str(snapshot / "clpoly/variable.cc")]
        command += shlex.split(os.environ.get("LDFLAGS", "-lflint -lgmpxx -lgmp"))
        command += ["-o", str(snapshot / "benchmark")]
        print(f"Building {version}", flush=True)
        with (snapshot / "build.log").open("w") as output:
            subprocess.run(command, cwd=ROOT, stdout=output, stderr=subprocess.STDOUT, check=True)
        commands[version] = {"command": command,
            "binary_sha256": hashlib.sha256((snapshot / "benchmark").read_bytes()).hexdigest()}
    (OUT / "build-manifest.json").write_text(json.dumps(commands, indent=2) + "\n")


def run(timeout):
    records = []
    for mode in ["quadratic", "nonmonic", "hard"]:
        order = VERSIONS if mode != "nonmonic" else list(reversed(VERSIONS))
        for version in order:
            for selected_r in [24, 32]:
                logfile = OUT / f"{mode}-{version}-r{selected_r}.log"
                print(f"Running {mode} {version} r={selected_r}", flush=True)
                started = time.monotonic()
                with logfile.open("w") as output:
                    try:
                        result = subprocess.run([str(OUT / version / "benchmark"), mode, str(selected_r)], cwd=ROOT,
                            stdout=output, stderr=subprocess.STDOUT, timeout=timeout)
                        status = "passed" if result.returncode == 0 else f"exit_{result.returncode}"
                    except subprocess.TimeoutExpired:
                        status = "timeout"
                elapsed = time.monotonic() - started
                text = logfile.read_text()
                samples = re.findall(r"^r=(24|32) repeat=(\d+) ms=([\d.]+)", text, re.M)
                reference = re.findall(r"^flint=3.1.2 r=(24|32) repeat=(\d+) ms=([\d.]+)", text, re.M)
                phases = {}
                for count, repeat, fields in re.findall(r"^phases r=(24|32) repeat=(\d+) (.+)$", text, re.M):
                    phases[(count, repeat)] = {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", fields)}
                if status == "passed":
                    assert len(samples) == len(reference) == len(phases) == 3, text
                    expected = {(str(selected_r), str(i)) for i in range(3)}
                    assert {(r, i) for r, i, _ in samples} == expected, text
                    assert {(r, i) for r, i, _ in reference} == expected, text
                    assert set(phases) == expected, text
                    if version == "repaired":
                        assert text.count("certified=1 fallback=0") == 3, text
                        assert all(p["lll_calls"] > 0 and p["inner_zass_calls"] == 0 and
                                   p["outer_zass_calls"] == 0 for p in phases.values()), text
                for count in [selected_r]:
                    values = [float(ms) for r, _, ms in samples if int(r) == count]
                    refs = [float(ms) for r, _, ms in reference if int(r) == count]
                    records.append({"mode": mode, "version": version, "r": count, "status": status,
                        "elapsed_seconds": elapsed, "timeout_seconds": timeout,
                        "samples_ms": values, "median_ms": statistics.median(values) if len(values) == 3 else None,
                        "flint_samples_ms": refs,
                        "phases": [fields for (r, _), fields in phases.items() if int(r) == count]})
                (OUT / "results.json").write_text(json.dumps(records, indent=2) + "\n")
                print(f"Finished {mode} {version} r={selected_r}: {status}, {len(samples)} validated samples", flush=True)
                if version == "repaired" and status != "passed":
                    raise SystemExit(f"repaired benchmark failed: {mode} {status}; inspect {logfile}")
    print(json.dumps(records, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--build-only", action="store_true")
    modes.add_argument("--run-only", action="store_true")
    parser.add_argument("--timeout", type=float, default=180,
                        help="seconds per version/input family/r; timed-out runs are explicitly censored")
    args = parser.parse_args()
    if not args.run_only:
        build()
    if not args.build_only:
        run(args.timeout)
