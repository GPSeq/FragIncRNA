// Compile against a frozen pre-optimization source snapshot. This driver uses
// that snapshot's encoder/archive reader/verifier without editing its sources.
// Scope: the documented k20/canonical/include/skip production benchmark only.
#ifndef FRAG_BASELINE_SOURCE
#error Define FRAG_BASELINE_SOURCE to the frozen src/kmer_v2.cpp
#endif
#define main frozen_application_main
#include FRAG_BASELINE_SOURCE
#undef main

int main(int argc,char ** argv) {
    if (argc != 6) { std::cerr << "manifest queries candidates raw_reference_cache output\n"; return 2; }
    try {
        Options o; o.k=20; o.fragment_size=80000; o.mode="ibf"; o.scoring_method="per_bin_verified";
        o.strand="canonical"; o.repeat="include"; o.fdr_control_per_bin=true; o.actual_threads=15;
        ReferenceInfo info; struct Archive { std::size_t begin; std::string file; }; std::vector<Archive> archives;
        std::ifstream manifest(argv[1]); std::string line;
        for (int i=0;i<4;++i) std::getline(manifest,line);
        while (std::getline(manifest,line) && !line.starts_with("shard_id\t")) {
            auto f=split_tabs(line); if (f.size()!=7) throw std::runtime_error("Invalid frozen manifest mapping");
            if (info.ids.empty() || info.ids.back()!=f[3]) { info.ids.push_back(f[3]); info.lengths.push_back(0); }
            info.lengths.back()=std::max(info.lengths.back(),static_cast<std::size_t>(std::stoull(f[5]))); ++info.bin_count;
        }
        while (std::getline(manifest,line)) { auto f=split_tabs(line); archives.push_back({std::stoull(f[1]),f[3]}); }
        std::vector<Query> batch; std::unordered_map<std::string,std::size_t> queries;
        std::string id,sequence; RepeatBed bed;
        auto flush=[&] { if (id.empty()) return; queries[id]=batch.size(); batch.push_back(make_query(id,batch.size(),sequence,o,nullptr,bed)); batch.back().candidate_support.resize(batch.back().windows.size()); sequence.clear(); };
        scan_fasta(argv[2],[&](auto const & next){flush();id=next;},[&](std::size_t,char c){sequence+=c;}); flush();
        std::ifstream candidates(argv[3]); std::getline(candidates,line); auto columns=split_tabs(line); std::map<std::string,std::size_t> ci;
        for (std::size_t i=0;i<columns.size();++i) ci[columns[i]]=i;
        while (std::getline(candidates,line)) {
            auto f=split_tabs(line); auto & q=batch.at(queries.at(f[ci.at("query_id")]));
            auto contig=std::find(info.ids.begin(),info.ids.end(),f[ci.at("contig")])-info.ids.begin();
            q.per_bin_regions.push_back({static_cast<std::size_t>(contig),std::stoull(f[ci.at("start")]),std::stoull(f[ci.at("end")]),0,0,0});
            auto bins=f[ci.at("source_bin_ids")]; std::replace(bins.begin(),bins.end(),',',' '); std::istringstream in(bins); std::size_t b;
            while (in>>b) q.per_bin_candidates.push_back({b,0,0,0,0});
        }
        for (auto const & shard:archives) {
            auto path=fs::path(argv[1]).parent_path()/shard.file; std::ifstream header(path,std::ios::binary);
            header.seekg(index_magic.size()); uint64_t size=0; header.read(reinterpret_cast<char *>(&size),sizeof(size)); header.seekg(64,std::ios::cur);
            if (size>64*1024*1024) throw std::runtime_error("Invalid metadata size");
            std::string bytes(size,'\0');header.read(bytes.data(),size);std::istringstream encoded(bytes,std::ios::binary);IndexMetadata stored;
            {cereal::BinaryInputArchive archive(encoded);archive(stored);}
            ReferenceInfo local;local.bin_count=stored.bin_count;
            auto ibf=load_index_file(path,stored,local,0);
            fraglncrna::parallel_records(batch.size(),15,[&](std::size_t qi){
                auto & q=batch[qi];auto agent=ibf->membership_agent();std::vector<std::size_t> selected;
                for(auto const & c:q.per_bin_candidates)if(c.bin>=shard.begin&&c.bin<shard.begin+local.bin_count)selected.push_back(c.bin-shard.begin);
                for(std::size_t wi=0;wi<q.windows.size();++wi)if(q.windows[wi].valid){
                    auto const & bits=agent.bulk_contains(index_key(q.windows[wi],o));q.bin_hits[wi]+=std::ranges::count(bits,true);
                    for(auto b:selected)if(bits[b]){q.candidate_support[wi]=1;break;}
                }
            });
        }
        IndexedReference reference{argv[4]};ExactResult result;
        verify_candidate_intervals(o,info,bed,batch,reference,result);
        std::ofstream out(argv[5]);out<<"query_index\tquery_id\tk_index\tvalid\tibf_any_bin_positive\tibf_candidate_positive\texact_verified_in_candidates\n";
        for(std::size_t qi=0;qi<batch.size();++qi){auto const & q=batch[qi];for(std::size_t wi=0;wi<q.windows.size();++wi)
            out<<q.index<<'\t'<<q.id<<'\t'<<wi<<'\t'<<q.windows[wi].valid<<'\t'<<bool(q.bin_hits[wi])<<'\t'<<unsigned(q.candidate_support[wi])<<'\t'<<bool(result.candidate_window_masks[qi][wi/64]&(uint64_t{1}<<(wi%64)))<<'\n';}
        out.close();if(!out)throw std::runtime_error("Cannot finalize oracle output");
    } catch(std::exception const & e){std::cerr<<e.what()<<'\n';return 1;}
}
