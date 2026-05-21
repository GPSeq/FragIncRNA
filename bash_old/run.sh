mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j

./lncrna_mers \
  -r /mnt/e/PrimateGenomes/prim/references/ \
  -q /mnt/e/PrimateGenomes/prim/gencode.v37.lncRNA_transcripts.fa \
  -f 8000 \
  -k 15 \
  -H 3 \
  -t 13 \
  -O ./ibf_results \
  --output-file results.tsv \
  --log-file ibf_run.log \
  --store-ibf \
  --cleanup-ibf
  #--store-fragments


./lncrna_mers   -r /mnt/e/PrimateGenomes/prim/references/   -q /mnt/e/PrimateGenomes/prim/gencode.v37.lncRNA_transcripts.fa   -f 1000   -k 17   -H 2   -t 15   -O ./ibf_results   --output-file results.tsv   --log-file ibf_run.log   --store-ibf   --cleanup-ibf --single-results-writer false
./lncrna_mers   -r /mnt/e/PrimateGenomes/all_primates/   -q /mnt/e/PrimateGenomes/prim/gencode.v37.lncRNA_transcripts.fa   -f 10000   -k 20   -H 2   -t 15   -O ./ibf_results   --output-file results.tsv   --log-file ibf_run.log   --store-ibf   --cleanup-ibf --single-results-writer false


*_count = total number of IBF hits (summed over all bins) for that query against that reference.

*_pass = 1 if count ≥ --hit-threshold, otherwise 0.

*_pct = count divided by the number of k-mers in the query (hits per k-mer, can be > 1).