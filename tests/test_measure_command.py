#!/usr/bin/env python3
"""The resource wrapper must include a force-killed timeout descendant."""
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile

measure = pathlib.Path(__file__).resolve().parents[1] / "scripts/measure_command.py"
with tempfile.TemporaryDirectory(prefix="frag-measure-") as directory:
    root = pathlib.Path(directory)
    normal = root / "normal.json"
    subprocess.run([sys.executable,str(measure),"--output",str(normal),"--",
                    sys.executable,"-c","x=bytearray(64*1024*1024)"],check=True)
    result = json.loads(normal.read_text())
    assert result["exit_code"] == 0 and result["peak_child_resident_bytes"] >= 64*1024*1024
    if sys.platform.startswith("linux") and shutil.which("timeout"):
        killed = root / "killed.json"
        child = "import signal,time; signal.signal(signal.SIGTERM,lambda *a:None); x=bytearray(64*1024*1024); time.sleep(10)"
        process = subprocess.run([sys.executable,str(measure),"--output",str(killed),"--",
            "timeout","--kill-after=0.3","0.4",sys.executable,"-c",child],check=False)
        result = json.loads(killed.read_text())
        assert process.returncode != 0 and result["exit_code"] != 0
        assert result["linux_sampled_descendant_peak_bytes"] >= 64*1024*1024, result
        assert result["peak_child_resident_bytes"] >= 64*1024*1024
        assert result["wall_seconds"] < 5
print("Resource measurement normal/forced-kill checks pass")
