#pragma once

void query_shards(ShardedIbf & shards, std::vector<Query> & batch, Options const & o, ReferenceInfo const & info) {
    // Global H_b uses occurrence multiplicity; only membership work is deduplicated.
    auto phase = Clock::now();
    for (std::size_t id = 0; id < shards.size(); ++id) {
        shards.visit(id,[&](Ibf const & ibf) {
            auto begin = shards.begin(id), bins = shards.bins(id);
            fraglncrna::parallel_records(batch.size(),static_cast<int>(o.actual_threads),[&](std::size_t qi) {
                auto & q = batch[qi]; auto agent = ibf.membership_agent();
                std::vector<std::size_t> counts(bins,0);
                using Bits = std::remove_cvref_t<decltype(agent.bulk_contains(0))>;
                Bits union_bits;
                std::vector<std::size_t> hits; // Calibration only; capacity reused.
                for (std::size_t ui = 0; ui < q.unique_windows.size(); ++ui) {
                    auto wi = q.unique_windows[ui]; auto const & w = q.windows[wi];
                    auto const & first = agent.bulk_contains(index_key(w,o));
                    auto const * bits = &first;
                    if (o.strand == "both" && w.forward != w.reverse) {
                        union_bits = first;
                        auto const & second = agent.bulk_contains(w.reverse);
                        for (std::size_t b = 0; b < bins; ++b) union_bits[b] = union_bits[b] || second[b];
                        bits = &union_bits;
                    }
                    auto multiplicity = q.unique_multiplicity[ui]; std::size_t positive = 0;
                    if (o.mode == "calibrate") hits.clear();
                    for (std::size_t b = 0; b < bins; ++b) if ((*bits)[b]) {
                        ++positive;
                        if (o.fdr_control_per_bin) counts[b] += multiplicity;
                        if (o.mode == "calibrate") hits.push_back(begin+b);
                    }
                    q.bin_hits[wi] += positive;
                    if (positive && multiplicity > (std::numeric_limits<std::uint64_t>::max()-q.positive_bin_pairs)/positive)
                        throw std::runtime_error("Logical IBF pair overflow query=" + q.id + " shard=" + std::to_string(id));
                    q.positive_bin_pairs += std::uint64_t(positive)*multiplicity;
                    if (o.mode == "calibrate") {
                        if (q.positive_bin_pairs > o.max_ibf_bin_pairs)
                            throw std::runtime_error("IBF-positive bin pairs exceed --max-ibf-bin-pairs for calibration query=" + q.id +
                                " index=" + std::to_string(q.index) + " valid_windows=" + std::to_string(q.valid_windows) + " shard=" + std::to_string(id));
                        // Calibration intentionally retains each original occurrence.
                        for (std::size_t at = q.calibration_offsets[ui]; at < q.calibration_offsets[ui+1]; ++at) {
                            auto original = q.calibration_occurrences[at];
                            q.positive_bins[original].insert(q.positive_bins[original].end(),hits.begin(),hits.end());
                        }
                    }
                    if ((ui & 4095) == 0) fraglncrna::check_stop();
                }
                // Dense shard scratch removes a hash-table update per positive pair.
                if (o.fdr_control_per_bin) for (std::size_t b = 0; b < bins; ++b)
                    if (counts[b]) q.per_bin_counts[begin+b] += counts[b];
            });
        });
        fraglncrna::production_log.progress("shard-query","evidence_shards=" + std::to_string(id+1) +
            "/" + std::to_string(shards.size()) + " block_queries=" + std::to_string(batch.size()));
    }
    shards.first_pass_seconds += elapsed(phase); phase = Clock::now();
    fraglncrna::parallel_records(batch.size(),static_cast<int>(o.actual_threads),[&](std::size_t qi) {
        auto & q = batch[qi];
        for (std::size_t wi = 0; wi < q.windows.size(); ++wi) if (q.windows[wi].valid)
            q.bin_hits[wi] = q.bin_hits[q.unique_windows[q.window_to_unique[wi]]];
        if (q.positive_bin_pairs > o.max_ibf_bin_pairs && o.mode != "calibrate")
            fraglncrna::production_log.event("shard-query","workload","high_logical_pairs query=" + q.id +
                " index=" + std::to_string(q.index) + " valid=" + std::to_string(q.valid_windows) +
                " unique=" + std::to_string(q.unique_windows.size()) + " logical_pairs=" + std::to_string(q.positive_bin_pairs));
        prepare_per_bin_query(q,o,info,nullptr);
    });
    shards.selection_seconds += elapsed(phase); phase = Clock::now();
    if (!o.fdr_control_per_bin) return;
    // Reverse order reuses the most recently resident shards immediately.
    for (std::size_t reverse = shards.size(); reverse > 0; --reverse) {
        auto id = reverse-1, begin = shards.begin(id), end = begin+shards.bins(id);
        bool needed = false;
        for (auto const & q : batch) for (auto const & c : q.per_bin_candidates)
            needed |= c.bin >= begin && c.bin < end;
        if (!needed) continue;
        shards.visit(id,[&](Ibf const & ibf) {
            fraglncrna::parallel_records(batch.size(),static_cast<int>(o.actual_threads),[&](std::size_t qi) {
                auto & q = batch[qi]; std::vector<std::size_t> selected;
                for (auto const & c : q.per_bin_candidates) if (c.bin >= begin && c.bin < end) selected.push_back(c.bin-begin);
                if (selected.empty()) return;
                auto agent = ibf.membership_agent();
                for (auto wi : q.unique_windows) {
                    auto const & w = q.windows[wi]; if (q.candidate_support[wi]) continue;
                    bool positive = false;
                    auto const & first = agent.bulk_contains(index_key(w,o));
                    for (auto bin : selected) if (first[bin]) { positive = true; break; }
                    if (!positive && o.strand == "both" && w.forward != w.reverse) {
                        auto const & second = agent.bulk_contains(w.reverse);
                        for (auto bin : selected) if (second[bin]) { positive = true; break; }
                    }
                    q.candidate_support[wi] |= positive;
                }
            });
        });
    }
    for (auto & q : batch) for (std::size_t wi = 0; wi < q.windows.size(); ++wi) if (q.windows[wi].valid)
        q.candidate_support[wi] = q.candidate_support[q.unique_windows[q.window_to_unique[wi]]];
    shards.second_pass_seconds += elapsed(phase);
}
