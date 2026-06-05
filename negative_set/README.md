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

The tool reads the first `chromosome_count` records from `genome_file`, samples `sequence_count` windows with lengths between `min_length` and `max_length`, shuffles each sampled sequence, fragments the selected chromosome records using `fragment_size` and the same 3 bp overlap as the main FraglncRNA tool, builds an IBF from those chromosome fragments, searches the shuffled negatives against that IBF, and writes FASTA plus TSV outputs.

`sequence_count` controls only the number of generated negative query sequences; it does not control the number of IBF bins. The number of bins is determined by `chromosome_count`, chromosome lengths, and `fragment_size`.
