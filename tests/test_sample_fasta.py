"""Check deterministic, order-preserving FASTA reservoir sampling."""
import gzip
import json
import pathlib
import subprocess
import sys
import tempfile


sampler = pathlib.Path(sys.argv[1])
with tempfile.TemporaryDirectory(prefix="fraglncrna-sampler-test-") as directory:
    root = pathlib.Path(directory)
    source = root / "input.fa.gz"
    with gzip.open(source, "wt") as out:
        for i in range(20):
            out.write(f">record_{i} description\nAcGTNN{i % 10}\n")
    samples = []
    for name in ("a", "b"):
        output = root / f"{name}.fa.gz"
        subprocess.run([sys.executable, str(sampler), "--input", str(source),
                        "--output", str(output), "--count", "5", "--seed", "42"],
                       check=True, capture_output=True, text=True)
        with gzip.open(output, "rt") as handle:
            samples.append(handle.read())
        metrics = json.loads(pathlib.Path(str(output) + ".metrics.json").read_text())
        assert metrics["sampled_records"] == 5
        assert metrics["total_records"] == 20
        assert metrics["peak_resident_bytes"] > 0
        assert metrics["wall_seconds"] >= 0
    assert samples[0] == samples[1]
    headers = [line for line in samples[0].splitlines() if line.startswith(">")]
    indices = [int(header.split()[0].split("_")[1]) for header in headers]
    assert len(indices) == 5 and indices == sorted(indices)
    assert "AcGTNN" in samples[0]

print("FASTA sampling checks passed")
