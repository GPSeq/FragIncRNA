"""Schema-3 merge compatibility: reject mixed scientific scoring policies."""
import csv
import json
import pathlib
import subprocess
import sys
import tempfile

program = pathlib.Path(sys.argv[1]).resolve()
script = pathlib.Path(__file__).resolve().parents[1]/"scripts"/"merge_kmer_v2.py"
with tempfile.TemporaryDirectory(prefix="fraglncrna-merge-production-") as directory:
    root = pathlib.Path(directory)
    queries = root/"q.fa"; queries.write_text(">q\nACGTACGTACGT\n>N\nNNNN\n")
    inputs = []
    for label in ("a","b"):
        reference = root/(label+".fa"); reference.write_text(">chr\nACGTACGTACGT\n")
        output = root/(label+".tsv")
        subprocess.run([str(program),"--reference",str(reference),"--queries",str(queries),"--output",str(output),
            "--mode","ibf","--scoring-method","per_bin_verified","--k","4","--fragment-size","8"],
            capture_output=True,text=True,check=True)
        inputs.append(output)
    command = [sys.executable,str(script),"--input",str(inputs[0]),"--input",str(inputs[1]),
               "--output",str(root/"merged.tsv"),"--summary",str(root/"summary.json")]
    completed = subprocess.run(command,capture_output=True,text=True)
    assert completed.returncode == 0, completed.stderr
    with (root/"merged.tsv").open() as handle:
        rows = list(csv.DictReader(handle,delimiter="\t"))
    assert len(rows) == 4 and all(row["schema_version"] == "3" for row in rows)
    assert all(row["primary_match_source"] == "candidate_exact" for row in rows)
    manifest_path = pathlib.Path(str(inputs[1])+".metrics.json")
    original = manifest_path.read_text()
    for key,value in (("per_bin_top_k",32),("scoring_method","per_bin_ibf"),("schema_version",2)):
        changed = json.loads(original); changed[key] = value
        manifest_path.write_text(json.dumps(changed))
        failed = subprocess.run(command,capture_output=True,text=True)
        assert failed.returncode != 0, (key,failed.stdout,failed.stderr)
        assert "differ" in failed.stderr
    manifest_path.write_text(original)
    text = inputs[1].read_text()
    inputs[1].write_text(text.replace("candidate_exact","ibf_candidates",1))
    failed = subprocess.run(command,capture_output=True,text=True)
    assert failed.returncode != 0 and "provenance differ" in failed.stderr
    inputs[1].write_text(text)
print("schema-3 merge compatibility checks passed")
