# negative_set

Standalone executable for generating shuffled negative sequences from one genome FASTA and querying them with the FraglncRNA IBF implementation. The genome input can be plain FASTA or gzip-compressed FASTA, including `.fa`, `.fna`, `.fasta`, `.fa.gz`, `.fna.gz`, and `.fasta.gz`.

Build from the repository root, after building the FraglncRNA tool, the FraglncRNA supports automatic building of the negative set generation, which is then imported here:

```bash
cmake -S . -B build
cmake --build build --target negative_set
```

Run:

```bash
./build/negative_set/negative_set negative_set/config.toml
```

The tool reads the first `chromosome_count` records from `genome_file`, samples `sequence_count` windows with lengths between `min_length` and `max_length`, shuffles each sampled sequence, builds an IBF from the original sampled windows, searches the shuffled negatives against that IBF, and writes FASTA plus TSV outputs.
