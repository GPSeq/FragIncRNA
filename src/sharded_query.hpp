#pragma once

void query_shards(ShardedIbf & shards, std::vector<Query> & batch, Options const & o, ReferenceInfo const & info) {
    // Block-major hybrid. Only sparse H_b and per-window union counts survive
    // each shard. Query order and window identities never depend on threads.
    for (std::size_t id = 0; id < shards.size(); ++id) {
        shards.visit(id,[&](Ibf const & ibf) {
            auto begin = shards.begin(id), bins = shards.bins(id);
            fraglncrna::parallel_records(batch.size(),static_cast<int>(o.actual_threads),[&](std::size_t qi) {
                auto & q = batch[qi]; auto agent = ibf.membership_agent();
                for (std::size_t wi = 0; wi < q.windows.size(); ++wi) {
                    auto const & w = q.windows[wi]; if (!w.valid) continue;
                    auto hits = positive_bins(agent,w,o);
                    for (auto bin : hits) {
                        if (bin >= bins) throw std::runtime_error("Nonzero padded shard bin");
                        ++q.bin_hits[wi];
                        if (o.mode == "calibrate" || o.fdr_control_per_bin) {
                            if (++q.positive_bin_pairs > o.max_ibf_bin_pairs)
                                throw std::runtime_error("IBF-positive bin pairs exceed --max-ibf-bin-pairs for one query");
                        }
                        if (o.fdr_control_per_bin) ++q.per_bin_counts[begin+bin];
                        if (o.mode == "calibrate") q.positive_bins[wi].push_back(begin+bin);
                    }
                    if ((wi+1)%65536 == 0) fraglncrna::check_stop();
                }
            });
        });
        fraglncrna::production_log.progress("shard-query","evidence_shards=" + std::to_string(id+1) +
            "/" + std::to_string(shards.size()) + " block_queries=" + std::to_string(batch.size()),true);
    }
    fraglncrna::parallel_records(batch.size(),static_cast<int>(o.actual_threads),[&](std::size_t qi) {
        prepare_per_bin_query(batch[qi],o,info,nullptr);
    });
    // One global selection precedes the second probe. Load each selected shard
    // once per block, never once per transcript. Store only a window bit/byte.
    if (!o.fdr_control_per_bin) return;
    for (std::size_t id = 0; id < shards.size(); ++id) {
        auto begin = shards.begin(id), end = begin+shards.bins(id);
        bool needed = false;
        for (auto const & q : batch) for (auto const & c : q.per_bin_candidates)
            needed |= c.bin >= begin && c.bin < end;
        if (!needed) continue;
        shards.visit(id,[&](Ibf const & ibf) {
            fraglncrna::parallel_records(batch.size(),static_cast<int>(o.actual_threads),[&](std::size_t qi) {
                auto & q = batch[qi]; std::vector<std::size_t> selected;
                for (auto const & c : q.per_bin_candidates) if (c.bin >= begin && c.bin < end) selected.push_back(c.bin-begin);
                if (selected.empty()) return;
                auto agent = ibf.membership_agent(); std::unordered_map<uint64_t,bool> memo;
                for (std::size_t wi = 0; wi < q.windows.size(); ++wi) {
                    auto const & w = q.windows[wi]; if (!w.valid || q.candidate_support[wi]) continue;
                    auto [it,inserted] = memo.try_emplace(index_key(w,o),false);
                    if (inserted) {
                        auto const & first = agent.bulk_contains(it->first);
                        for (auto bin : selected) if (first[bin]) { it->second = true; break; }
                        if (!it->second && o.strand == "both" && w.forward != w.reverse) {
                            auto const & second = agent.bulk_contains(w.reverse);
                            for (auto bin : selected) if (second[bin]) { it->second = true; break; }
                        }
                    }
                    q.candidate_support[wi] |= it->second;
                }
            });
        });
    }
}
