"""Black-box TOML coverage, CLI equivalence, and strict error checks for v2."""

import csv
import hashlib
import json
import pathlib
import re
import subprocess
import sys
import tempfile


PROGRAM = pathlib.Path(sys.argv[1])
ROOT = pathlib.Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "kmer_v2"
EXAMPLE = ROOT / "configs" / "example.toml"


example_text = EXAMPLE.read_text()
assert '[v2]' in example_text and 'mode = "exact"' in example_text
source = (ROOT / "src" / "kmer_v2.cpp").read_text()
cli_keys = {key for key in re.findall(r'(?:if|else if) \(key == "--([^"]+)"\)', source)
            if key not in {"help", "config"} and not key.startswith("no-")}
documented_keys = set(re.findall(r"^\s*#?\s*([a-z][a-z0-9_]*)\s*=", example_text, re.M))
assert {key.replace("-", "_") for key in cli_keys} <= documented_keys


def run(*args, good=True):
    result = subprocess.run([str(PROGRAM), *map(str, args)], capture_output=True, text=True)
    if good:
        assert result.returncode == 0, result.stderr
    else:
        assert result.returncode != 0
    return result


with tempfile.TemporaryDirectory(prefix="fraglncrna-v2-toml-") as directory:
    tmp = pathlib.Path(directory)
    config = tmp / "complete.toml"
    out_config = tmp / "from_config.tsv"
    out_cli = tmp / "from_cli.tsv"
    paths = {
        "positions": (tmp / "config_positions.tsv", tmp / "cli_positions.tsv"),
        "fragments": (tmp / "config_fragments.tsv", tmp / "cli_fragments.tsv"),
        "repeat_classes": (tmp / "config_query_classes.tsv", tmp / "cli_query_classes.tsv"),
        "reference_hit_classes": (tmp / "config_reference_classes.tsv", tmp / "cli_reference_classes.tsv"),
        "discover_loci": (tmp / "config_loci.tsv", tmp / "cli_loci.tsv"),
    }
    lines = ["[v2]", f'reference = "{FIXTURE / "reference.fa"}"',
             f'queries = "{FIXTURE / "queries.fa"}"', f'output = "{out_config}"',
             'k = 4', 'fragment_size = 4', 'mode = "exact"', 'repeat = "annotate"',
             'local_window = 4', 'batch_windows = 1_000']
    for name, (config_path, _) in paths.items():
        lines.append(f'{name} = "{config_path}"')
    config.write_text("\n".join(lines) + "\n")
    run("--config", config)
    cli = ["--reference", FIXTURE / "reference.fa", "--queries", FIXTURE / "queries.fa",
           "--output", out_cli, "--k", "4", "--fragment-size", "4", "--mode", "exact",
           "--repeat", "annotate", "--local-window", "4", "--batch-windows", "1000"]
    for name, (_, cli_path) in paths.items():
        cli.extend(("--" + name.replace("_", "-"), cli_path))
    run(*cli)
    assert out_config.read_bytes() == out_cli.read_bytes()
    for config_path, cli_path in paths.values():
        assert config_path.read_bytes() == cli_path.read_bytes()
    metrics = json.loads(pathlib.Path(str(out_config) + ".metrics.json").read_text())
    assert metrics["config_path"] == str(config)
    assert metrics["config_sha256"] == hashlib.sha256(config.read_bytes()).hexdigest()
    assert metrics["total_seconds"] >= 0 and metrics["peak_resident_bytes"] > 0

    override = tmp / "overridden.tsv"
    run("--config", config, "--output", override, "--strand", "forward")
    overridden = json.loads(pathlib.Path(str(override) + ".metrics.json").read_text())
    assert overridden["strand_mode"] == "forward"
    assert overridden["config_path"] == str(config)
    with override.open(newline="") as handle:
        assert all(row["strand_mode"] == "forward" for row in csv.DictReader(handle, delimiter="\t"))

    calibration = tmp / "calibration.toml"
    bins = tmp / "calibration_bins.tsv"
    calibration.write_text(
        f'[v2]\nreference = "{FIXTURE / "reference.fa"}"\n'
        f'queries = "{FIXTURE / "queries.fa"}"\noutput = "{tmp / "calibration.tsv"}"\n'
        f'calibration_bins = "{bins}"\nk = 4\nfragment_size = 4\n'
        'mode = "calibrate"\ngenome_fpr = 0.01\n')
    run(calibration)  # A lone TOML path is accepted, like the legacy executable.
    calibrated = json.loads((tmp / "calibration.tsv.metrics.json").read_text())
    assert calibrated["target_genome_wide_fpr"] == 0.01
    assert calibrated["configured_per_bin_fpr"] < 0.01
    assert bins.exists()

    per_bin_config = tmp / "per_bin.toml"
    per_bin_config.write_text(
        f'[v2]\nreference = "{FIXTURE / "reference.fa"}"\n'
        f'queries = "{FIXTURE / "queries.fa"}"\noutput = "{tmp / "per_bin_config.tsv"}"\n'
        'k = 4\nfragment_size = 4\nmode = "validate"\n'
        'fdr_control_per_bin = true\nper_bin_exact_verify = true\n'
        'per_bin_min_hits = 1\nper_bin_min_fraction = 0.25\n')
    run("--config", per_bin_config)
    per_bin_cli = tmp / "per_bin_cli.tsv"
    run("--reference", FIXTURE / "reference.fa", "--queries", FIXTURE / "queries.fa",
        "--output", per_bin_cli, "--k", "4", "--fragment-size", "4", "--mode", "validate",
        "--fdr-control-per-bin", "--per-bin-min-hits", "1", "--per-bin-min-fraction", "0.25")
    assert (tmp / "per_bin_config.tsv.per_bin.tsv").read_bytes() == pathlib.Path(
        str(per_bin_cli) + ".per_bin.tsv").read_bytes()
    assert (tmp / "per_bin_config.tsv.per_bin.tsv.candidates.tsv").read_bytes() == pathlib.Path(
        str(per_bin_cli) + ".per_bin.tsv.candidates.tsv").read_bytes()
    disabled = tmp / "per_bin_disabled.tsv"
    run("--config", per_bin_config, "--no-fdr-control-per-bin", "--output", disabled)
    assert not pathlib.Path(str(disabled) + ".per_bin.tsv").exists()
    false_config = tmp / "per_bin_false.toml"
    false_config.write_text(per_bin_config.read_text().replace("fdr_control_per_bin = true",
                                                          "fdr_control_per_bin = false"))
    enabled = tmp / "per_bin_enabled.tsv"
    run("--config", false_config, "--fdr-control-per-bin", "--output", enabled)
    assert pathlib.Path(str(enabled) + ".per_bin.tsv").exists()

    for bad, expected in (
        ('[v2]\nunknown = 1\n', "Unknown v2 TOML key"),
        ('[v2]\nk = 4\nk = 5\n', "Duplicate TOML key"),
        ('[v2]\nmode = 4\n', "requires a quoted string"),
        ('[v2]\nk = "4"\n', "requires a nonnegative number"),
        ('[v2]\nk = 1__000\n', "Invalid TOML number"),
        ('[v2]\nmode = "ex"act"\n', "Unescaped quote"),
        ('[v2]\nfdr_control_per_bin = "true"\n', "requires true or false"),
        ('[v2]\nfdr_control_per_bin = 1\n', "requires true or false"),
        ('[general]\nreference = "x"\n', "Expected one [v2]"),
    ):
        broken = tmp / "broken.toml"
        broken.write_text(bad)
        assert expected in run("--config", broken, good=False).stderr

print("kmer v2 TOML checks passed")
