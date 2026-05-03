#!/usr/bin/env python3
"""B2B 测试编排器

Spawn C++ driver + Lean driver 作为长驻 subprocess，逐条发 vector，对比响应。

用法：
    python3 runner.py <vector_file_or_dir> [...]

退出码：0=全 PASS，1=有 FAIL，2=runner 自身崩。
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import diff as diff_mod

B2B_ROOT = Path(__file__).resolve().parent.parent
LEAN_PROJECT = B2B_ROOT.parent / "lean"
CPP_DRIVER = B2B_ROOT / "cpp_driver"


class Driver:
    """长驻 subprocess 封装"""

    def __init__(self, label: str, cmd: list[str], cwd: Path):
        self.label = label
        self.proc = subprocess.Popen(
            cmd, cwd=cwd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,  # line-buffered
        )

    def call(self, req: dict) -> dict:
        """发一条请求，等响应。"""
        line = json.dumps(req)
        self.proc.stdin.write(line + "\n")
        self.proc.stdin.flush()
        resp_line = self.proc.stdout.readline()
        if not resp_line:
            err = self.proc.stderr.read()
            raise RuntimeError(
                f"{self.label} driver closed unexpectedly. stderr:\n{err}"
            )
        try:
            return json.loads(resp_line)
        except json.JSONDecodeError as e:
            raise RuntimeError(
                f"{self.label} bad response: {resp_line!r} ({e})"
            )

    def close(self):
        try:
            self.proc.stdin.close()
            self.proc.wait(timeout=5)
        except Exception:
            self.proc.kill()


def collect_vectors(paths: list[Path]) -> list[tuple[Path, dict]]:
    """递归展开 paths 为 vector 文件列表。"""
    out: list[tuple[Path, dict]] = []
    for p in paths:
        if p.is_dir():
            for f in sorted(p.rglob("*.json")):
                with open(f) as fh:
                    out.append((f, json.load(fh)))
        else:
            with open(p) as fh:
                out.append((p, json.load(fh)))
    return out


def run(args):
    if not CPP_DRIVER.exists():
        print(f"FATAL: {CPP_DRIVER} not found. Run `make -C proof/b2b cpp_driver` first.",
              file=sys.stderr)
        return 2

    print(f"[runner] starting C++ driver: {CPP_DRIVER}")
    cpp = Driver("cpp", [str(CPP_DRIVER)], cwd=B2B_ROOT)

    elan_path = os.path.expanduser("~/.elan/bin")
    env = os.environ.copy()
    env["PATH"] = f"{elan_path}:{env.get('PATH','')}"

    print(f"[runner] lake build B2B.Registry B2B.Driver (ensure freshness)")
    bres = subprocess.run(
        ["lake", "build", "B2B.Registry", "B2B.Driver"],
        cwd=LEAN_PROJECT, env=env, capture_output=True, text=True,
    )
    if bres.returncode != 0:
        print(f"FATAL: lake build failed:\n{bres.stderr}", file=sys.stderr)
        cpp.close()
        return 2

    print(f"[runner] starting Lean driver via interpreter (~2s startup)")
    lean_proc = subprocess.Popen(
        ["lake", "env", "lean", "--run", "B2B/Driver.lean"],
        cwd=LEAN_PROJECT,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
        env=env,
    )
    lean = Driver.__new__(Driver)
    lean.label = "lean"
    lean.proc = lean_proc

    vectors = collect_vectors([Path(p) for p in args.vectors])
    if not vectors:
        print("[runner] no vector files found", file=sys.stderr)
        cpp.close(); lean.close()
        return 2

    total = pass_n = fail_n = 0
    failures: list[str] = []

    for vec_path, vec in vectors:
        fn = vec.get("fn")
        cases = vec.get("cases", [])
        print(f"\n=== {vec_path.name} ({fn}) — {len(cases)} cases ===")
        for i, case in enumerate(cases):
            name = case.get("name", f"case_{i}")
            req_id = f"{fn}/{name}"
            req = {"id": req_id, "fn": fn, "args": case["args"]}
            try:
                cpp_resp = cpp.call(req)
                lean_resp = lean.call(req)
            except Exception as e:
                print(f"  ERROR  {name}: {e}")
                failures.append(f"{req_id}: driver error {e}")
                fail_n += 1
                total += 1
                continue
            status, detail = diff_mod.diff_one(
                cpp_resp, lean_resp, case.get("expected")
            )
            total += 1
            if status == "PASS":
                pass_n += 1
                print(f"  PASS   {name}")
            else:
                fail_n += 1
                print(f"  {status:<13} {name} — {detail}")
                failures.append(f"{req_id} [{status}]: {detail}")

    cpp.close()
    lean.close()

    print(f"\n=== Summary ===")
    print(f"  Total:  {total}")
    print(f"  PASS:   {pass_n}")
    print(f"  FAIL:   {fail_n}")
    if failures and args.verbose:
        print("\n--- Failures ---")
        for f in failures:
            print(f"  {f}")

    return 0 if fail_n == 0 else 1


def main():
    ap = argparse.ArgumentParser(description="B2B test runner")
    ap.add_argument("vectors", nargs="+", help="vector file(s) or dir(s)")
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="print failure details at end")
    args = ap.parse_args()
    sys.exit(run(args))


if __name__ == "__main__":
    main()
