
python3 scripts/null_models/prepare_null_models.py \
  --input-dir . \
  --output-dir generated_3000 \
  --cache-dir /tmp/fraglncrna_null_models_cache \
  --n-transcripts 3000 \
  --shuffles-per-transcript 10 \
  --intergenic-candidates 4096 \
  --seed 42 \
  --threads 15 \
  --k-values 15,20,32 \
  --allow-partial


python3 run_collection_null/run_null_models_ibf.py \
  --genomes-dir ../../genomes \
  --datasets-dir generated_3000 \
  --output-dir null_ibf_results \
  --binary ../FragIncRNA/build/lncrna_mers_v2 \
  --threads 15 \
  --parallel-genomes 1 \
  --k-values 15,20,32 \
  --hash-functions 3

