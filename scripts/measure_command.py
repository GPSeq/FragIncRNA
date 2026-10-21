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


def linux_descendant_peak(pid):
    """Largest sampled per-process HWM, not a sum of process-tree RSS.

    A killed timeout process may not reap its child, so RUSAGE_CHILDREN alone
    can omit the actual benchmark. Poll descendants while they still exist.
    """
    pending, seen, peak = [pid], set(), 0
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        root = pathlib.Path(f"/proc/{current}")
        try:
            for line in (root / "status").read_text().splitlines():
                if line.startswith("VmHWM:"):
                    peak = max(peak, int(line.split()[1]) * 1024)
            for task in (root / "task").iterdir():
                try:
                    pending.extend(map(int, (task / "children").read_text().split()))
                except (OSError, ValueError):
                    pass  # exiting thread/process
        except (OSError, ValueError):
            pass
    return peak


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("A command is required after --")

    started = time.perf_counter()
    process = subprocess.Popen(command)
    sampled_peak = 0
    while True:
        if sys.platform.startswith("linux"):
            sampled_peak = max(sampled_peak, linux_descendant_peak(process.pid))
        try:
            returncode = process.wait(timeout=0.05)
            break
        except subprocess.TimeoutExpired:
            pass
    seconds = time.perf_counter() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    # macOS reports bytes; Linux reports KiB.
    reaped_peak = int(usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024))
    rss_bytes = max(reaped_peak, sampled_peak)
    report = {"command": command, "exit_code": returncode,
              "wall_seconds": seconds, "peak_child_resident_bytes": rss_bytes,
              "reaped_child_peak_resident_bytes": reaped_peak,
              "linux_sampled_descendant_peak_bytes": sampled_peak if sys.platform.startswith("linux") else None,
              "rss_sampling_interval_seconds": 0.05 if sys.platform.startswith("linux") else None,
              "rss_measurement_scope": "maximum single-process HWM, not summed tree RSS; sampling can miss short-lived or orphaned children"}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wall={seconds:.3f}s peak_child_rss={rss_bytes}B exit={returncode}",
          file=sys.stderr)
    return returncode


if __name__ == "__main__":
    raise SystemExit(main())
