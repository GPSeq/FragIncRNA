#!/usr/bin/env python3
"""Run one command and record wall time and peak child RSS as JSON.

Example: python3 scripts/measure_command.py --output /tmp/build-metrics.json -- cmake --build build -j 4
"""
import argparse
import json
import pathlib
import resource
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("A command is required after --")

    started = time.perf_counter()
    result = subprocess.run(command, check=False)
    seconds = time.perf_counter() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    # macOS reports bytes; Linux reports KiB.
    rss_bytes = int(usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024))
    report = {"command": command, "exit_code": result.returncode,
              "wall_seconds": seconds, "peak_child_resident_bytes": rss_bytes}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wall={seconds:.3f}s peak_child_rss={rss_bytes}B exit={result.returncode}",
          file=sys.stderr)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
