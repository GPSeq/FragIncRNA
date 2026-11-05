#pragma once

void compare_sharded_validation_bins(Options const & o, std::vector<Query> const & batch,
                                      ExactResult const & exact, ShardedIbf & shards, Totals & totals) {
    for (auto const & q : batch) for (auto const & [bin,hits] : q.per_bin_counts) totals.bin_fp[bin] += hits;
    // Probe true exact pairs against their owning shard. No dense matrix and
    // no assumption that candidate intervals contain all exact truth.
    for (std::size_t id = 0; id < shards.size(); ++id) {
        auto begin = shards.begin(id), end = begin+shards.bins(id);
        shards.visit(id,[&](Ibf const & ibf) {
            auto agent = ibf.membership_agent();
            for (auto const & q : batch) for (auto const & w : q.windows) {
                if (!w.valid) continue;
                auto bins = exact_bins_for_window(w,o,exact);
                std::erase_if(bins,[&](auto bin) { return bin < begin || bin >= end; });
                if (bins.empty()) continue;
                std::vector<unsigned char> positive(bins.size());
                auto const & first = agent.bulk_contains(index_key(w,o));
                for (std::size_t i = 0; i < bins.size(); ++i) positive[i] = first[bins[i]-begin];
                if (o.strand == "both" && w.forward != w.reverse) {
                    auto const & second = agent.bulk_contains(w.reverse);
                    for (std::size_t i = 0; i < bins.size(); ++i) positive[i] |= second[bins[i]-begin];
                }
                for (std::size_t i = 0; i < bins.size(); ++i) {
                    auto bin = bins[i]; ++totals.bin_exact_present[bin];
                    if (positive[i]) {
                        if (!totals.bin_fp[bin]) throw std::runtime_error("Inconsistent sharded per-bin membership accounting");
                        --totals.bin_fp[bin];
                    } else ++totals.bin_fn[bin];
                }
            }
        });
    }
}
