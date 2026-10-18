"""The benchmark runner must reject a mislabeled full-sized query file."""
import pathlib
import subprocess
import sys
import tempfile

script = pathlib.Path(__file__).resolve().parents[1]/"scripts/run_per_bin_production_benchmark.py"
with tempfile.TemporaryDirectory(prefix="fraglncrna-runner-cap-") as directory:
    root = pathlib.Path(directory)
    (root/"reference_first_contig.fa").write_text(">r\nACGT\n")
    (root/"queries_100.fa").write_text("".join(f">q{i}\nACGT\n" for i in range(101)))
    result = subprocess.run([sys.executable,str(script),"--program",sys.executable,"--directory",str(root),
                             "--stage","smoke"],capture_output=True,text=True)
    assert result.returncode != 0
    assert "query count cap of 100" in result.stderr
    assert not (root/"cold_validate100.command.json").exists()
print("benchmark query-count safety check passed")
