"""Create isolated, identically instrumented headers for the three fixed baselines.

Only benchmark copies get timing scopes. Production sources remain untouched.
Build bench_vanhoeij_certificate.cc with CLPOLY_PROFILE, CLPOLY_PHASE_BENCH,
and CLPOLY_FLINT_REFERENCE, using the respective snapshot include directory.
The repaired version also supports CLPOLY_CERTIFICATE_BENCH.
"""
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "_build/vanhoeij-benchmark"
GIT = os.environ.get("GIT", "git")
VERSIONS = {
    "presafety": "eb45ca4845640bd38d0a4d507cda86c3988f1a83",
    "safetynet": "428a27fe2c22dcbeee100d1100e310d8a08f1e77",
    "repaired": None,
}
TIMERS = r'''
    // Benchmark-only scopes; this block never enters production sources.
    struct __bench_phase_stats {
        long long hensel_ns = 0, lll_ns = 0, inner_zass_ns = 0, outer_zass_ns = 0;
        long long hensel_calls = 0, lll_calls = 0, inner_zass_calls = 0, outer_zass_calls = 0;
    };
    inline __bench_phase_stats __bench_phases;
    inline unsigned __bench_vh_depth = 0;
    struct __bench_phase_timer {
        long long& elapsed;
        std::chrono::steady_clock::time_point start;
        __bench_phase_timer(long long& value, long long& calls)
            : elapsed(value), start(std::chrono::steady_clock::now()) { ++calls; }
        ~__bench_phase_timer() {
            elapsed += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - start).count();
        }
    };
    struct __bench_vh_context {
        __bench_vh_context() { ++__bench_vh_depth; }
        ~__bench_vh_context() { --__bench_vh_depth; }
    };
'''


def instrument_function(source, name, statement):
    definitions = []
    for match in re.finditer(r"\n\s*" + re.escape(name) + r"\s*\(", source):
        depth, pos = 1, match.end()
        while depth:
            if source[pos] == "(":
                depth += 1
            elif source[pos] == ")":
                depth -= 1
            pos += 1
        while source[pos].isspace():
            pos += 1
        if source[pos] == "{":
            definitions.append(pos + 1)
    assert len(definitions) == 1, (name, len(definitions))
    pos = definitions[0]
    return source[:pos] + "\n        " + statement + source[pos:]


def main():
    manifest = {}
    original_realroot = None
    for name, revision in VERSIONS.items():
        destination = OUT / name
        destination.mkdir(parents=True, exist_ok=True)
        if revision:
            archive = subprocess.check_output([GIT, "archive", revision, "clpoly"], cwd=ROOT)
            with tarfile.open(fileobj=io.BytesIO(archive)) as tf:
                tf.extractall(destination, filter="data")
            realroot = (destination / "clpoly/realroot.hh").read_bytes()
            if original_realroot is None:
                original_realroot = realroot
            else:
                assert realroot == original_realroot, "baseline realroot sources differ"
            # Same template-definition ordering fix in all builds, needed only
            # for Clang to parse the common FLINT conversion test header.
            shutil.copy2(ROOT / "clpoly/realroot.hh", destination / "clpoly/realroot.hh")
        else:
            shutil.copytree(ROOT / "clpoly", destination / "clpoly", dirs_exist_ok=True)
        header = destination / "clpoly/polynomial_factorize_univar.hh"
        source = header.read_text()
        original_hash = hashlib.sha256(source.encode()).hexdigest()
        assert source.count("namespace clpoly{") == 1
        source = source.replace("namespace clpoly{", "namespace clpoly{" + TIMERS, 1)
        for function, phase in [("__hensel_lift", "hensel"), ("__lll_reduce", "lll")]:
            source = instrument_function(source, function,
                f"__bench_phase_timer phase_timer(__bench_phases.{phase}_ns, __bench_phases.{phase}_calls);")
        source = instrument_function(source, "__vanhoeij_recombine", "__bench_vh_context phase_context;")
        source = instrument_function(source, "__zassenhaus_recombine", """__bench_phase_timer phase_timer(
            __bench_vh_depth ? __bench_phases.inner_zass_ns : __bench_phases.outer_zass_ns,
            __bench_vh_depth ? __bench_phases.inner_zass_calls : __bench_phases.outer_zass_calls);""")
        header.write_text(source)
        manifest[name] = {"revision": revision or "worktree",
                          "source_header_sha256": original_hash,
                          "instrumented_header_sha256": hashlib.sha256(source.encode()).hexdigest()}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
