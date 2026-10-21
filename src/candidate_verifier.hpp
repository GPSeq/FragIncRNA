#pragma once

void verify_candidate_intervals(Options const & o, ReferenceInfo const & info, RepeatBed const & reference_bed,
                                std::vector<Query> const & batch, IndexedReference const & reference,
                                ExactResult & result) {
    using fraglncrna::production_log;
    auto grouping_start = Clock::now();
    std::vector<std::vector<std::vector<uint64_t>>> masks(batch.size());
    result.per_bin_region_exact_windows.resize(batch.size());
    result.per_bin_candidate_union_windows.resize(batch.size());
    result.candidate_window_masks.resize(batch.size());
    std::size_t mask_words = 0;
    constexpr std::size_t max_words = 256ULL * 1024 * 1024 / sizeof(uint64_t);
    std::unordered_map<uint64_t, std::vector<std::pair<std::size_t,std::size_t>>> owners;
    struct Interval { std::size_t contig, start, end; };
    std::vector<Interval> intervals;
    std::vector<std::pair<std::size_t,std::size_t>> boundaries;
    std::vector<std::vector<std::size_t>> coordinate_regions(batch.size());
    std::vector<std::size_t> offsets(info.ids.size(), 0);
    for (std::size_t ci = 1; ci < offsets.size(); ++ci) offsets[ci] = offsets[ci-1] + info.lengths[ci-1];
    for (std::size_t qi = 0; qi < batch.size(); ++qi) {
        auto const & q = batch[qi];
        auto words = (q.windows.size()+63)/64;
        if (words && q.per_bin_regions.size()+1 > (max_words-mask_words)/words)
            throw std::runtime_error("Per-bin exact region masks exceed 256 MiB per batch; reduce --batch-windows or --per-bin-top-k");
        mask_words += words * (q.per_bin_regions.size()+1);
        result.candidate_window_masks[qi].resize(words);
        for (auto const & region : q.per_bin_regions) {
            masks[qi].emplace_back(words,0);
            intervals.push_back({region.contig,region.start,region.end});
            boundaries.emplace_back(region.contig,region.start);
            boundaries.emplace_back(region.contig,region.end);
        }
        auto & ordered = coordinate_regions[qi];
        ordered.resize(q.per_bin_regions.size()); std::iota(ordered.begin(),ordered.end(),0);
        std::sort(ordered.begin(),ordered.end(),[&](auto a, auto b) {
            auto const & x = q.per_bin_regions[a]; auto const & y = q.per_bin_regions[b];
            return std::tie(x.contig,x.start) < std::tie(y.contig,y.start);
        });
        for (std::size_t i = 1; i < ordered.size(); ++i) {
            auto const & previous = q.per_bin_regions[ordered[i-1]];
            auto const & current = q.per_bin_regions[ordered[i]];
            if (previous.contig == current.contig && previous.end > current.start)
                throw std::runtime_error("Overlapping candidate start regions violate normalized locus layout");
        }
        if (q.per_bin_regions.empty()) continue;
        for (std::size_t wi = 0; wi < q.windows.size(); ++wi) {
            auto const & w = q.windows[wi]; if (!w.valid) continue;
            owners[w.forward].emplace_back(qi,wi);
            if (o.strand != "forward" && w.reverse != w.forward) owners[w.reverse].emplace_back(qi,wi);
        }
    }
    // A cache-local rejection filter avoids two large hash-table probes for
    // almost every reference key. It can only admit extra exact lookups: every
    // owner key sets the same three bits used here, so it cannot reject a true
    // owner. It does not replace the exact map or change scientific FPR.
    auto filter_hash = [](uint64_t value) {
        value ^= value >> 30; value *= 0xbf58476d1ce4e5b9ULL;
        value ^= value >> 27; value *= 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    };
    auto filter_bits = [](uint64_t hash) {
        return (uint64_t{1} << (hash & 63)) | (uint64_t{1} << ((hash >> 6) & 63)) |
               (uint64_t{1} << ((hash >> 12) & 63));
    };
    std::size_t filter_words = 1;
    while (filter_words < owners.size()/4+1) filter_words *= 2;
    std::vector<uint64_t> rejection_filter(o.candidate_prefilter ? filter_words : 0,0);
    if (o.candidate_prefilter) for (auto const & [key,ignored] : owners) {
        auto hash = filter_hash(key); rejection_filter[(hash >> 18) & (filter_words-1)] |= filter_bits(hash);
    }
    // Fetch a reference start interval once per batch, even when many queries
    // or loci request it. Extension k-1 is bounded by the original contig.
    std::sort(intervals.begin(), intervals.end(), [](auto const & a, auto const & b) {
        return std::tie(a.contig,a.start,a.end) < std::tie(b.contig,b.start,b.end);
    });
    std::vector<Interval> merged;
    for (auto const & interval : intervals) {
        if (!merged.empty() && merged.back().contig == interval.contig && interval.start <= merged.back().end)
            merged.back().end = std::max(merged.back().end,interval.end);
        else merged.push_back(interval);
    }
    for (auto const & interval : merged) result.candidate_unique_start_bases += interval.end-interval.start;
    for (auto const & q : batch) for (auto const & region : q.per_bin_regions)
        result.candidate_logical_requested_bases += std::min(region.end+o.k-1,info.lengths[region.contig])-region.start;
    result.candidate_unique_intervals = merged.size();
    // Within each task, the set of query/locus owners is constant. Therefore
    // one valid occurrence of a key suffices for every owner in that task.
    // Splitting at every locus boundary makes bounded repeated-key caching
    // correct even for overlapping requests from different queries.
    std::sort(boundaries.begin(),boundaries.end());
    boundaries.erase(std::unique(boundaries.begin(),boundaries.end()),boundaries.end());
    std::vector<Interval> tasks;
    for (auto const & interval : merged)
        for (std::size_t start = interval.start; start < interval.end;) {
            auto end = start + std::min<std::size_t>(1048576,interval.end-start);
            auto next = std::upper_bound(boundaries.begin(),boundaries.end(),std::pair{interval.contig,start});
            if (next != boundaries.end() && next->first == interval.contig) end = std::min(end,next->second);
            tasks.push_back({interval.contig,start,end}); start = end;
        }
    for (auto const & task : tasks)
        result.candidate_reference_bases += std::min(task.end+o.k-1,info.lengths[task.contig])-task.start;
    result.candidate_interval_reads = tasks.size();
    result.candidate_grouping_seconds = std::chrono::duration<double>(Clock::now()-grouping_start).count();
    auto const & immutable_owners = std::as_const(owners);
    std::atomic<std::size_t> cache_hits{0}, cache_insertions{0}, prefilter_rejections{0}, full_lookups{0};
    // Private streams and rolling encoders. Only word-wise atomic OR touches
    // shared scientific state, then ordered aggregation after the barrier.
    // Wave size bounds progress latency to at most one 1-Mb task per worker.
    for (std::size_t base = 0; base < tasks.size(); base += o.actual_threads) {
        std::size_t count = std::min<std::size_t>(o.actual_threads,tasks.size()-base);
        fraglncrna::parallel_records(count,static_cast<int>(o.actual_threads),[&](std::size_t ti) {
            auto const & task = tasks[base+ti];
            std::ifstream input(reference.path(),std::ios::binary);
            if (!input) throw std::runtime_error("Cannot open indexed reference for candidate verification");
            input.seekg(static_cast<std::streamoff>(offsets[task.contig]+task.start));
            std::size_t end = std::min(task.end+o.k-1,info.lengths[task.contig]);
            WindowMaker maker{o,o.repeat == "exclude" && o.mask_target != "query"};
            BedCursor bed{reference_bed}; bed.reset(info.ids[task.contig]);
            std::array<char,65536> buffer;
            std::unordered_set<uint64_t> seen_keys;
            std::size_t local_cache_hits = 0, local_insertions = 0, local_rejections = 0, local_lookups = 0;
            constexpr std::size_t cache_limit = 8192;
            for (std::size_t pos = task.start; pos < end;) {
                fraglncrna::check_stop();
                auto n = std::min(buffer.size(),end-pos);
                reference.read(input,offsets[task.contig]+pos,buffer.data(),n);
                for (std::size_t j = 0; j < n; ++j,++pos) {
                    if ((pos & 4095) == 0) fraglncrna::check_stop();
                    auto w = maker.push(buffer[j],bed.at(pos));
                    if (!w || !w->valid) continue;
                    if (o.candidate_prefilter) {
                        auto hash = filter_hash(w->forward), bits = filter_bits(hash);
                        if ((rejection_filter[(hash >> 18) & (filter_words-1)] & bits) != bits) { ++local_rejections; continue; }
                    }
                    if (seen_keys.contains(w->forward)) { ++local_cache_hits; continue; }
                    ++local_lookups;
                    auto found = immutable_owners.find(w->forward);
                    if (found == immutable_owners.end()) continue;
                    std::size_t start = pos+1-o.k;
                    if (start >= task.end) continue;
                    auto const & key_owners = found->second;
                    for (std::size_t first = 0; first < key_owners.size();) {
                            auto qi = key_owners[first].first;
                            auto last = std::upper_bound(key_owners.begin()+first,key_owners.end(),qi,
                                [](auto query, auto const & owner) { return query < owner.first; });
                            std::size_t after = static_cast<std::size_t>(last-key_owners.begin());
                            std::size_t ri = std::numeric_limits<std::size_t>::max();
                            auto const & ordered = coordinate_regions[qi];
                            auto it = std::upper_bound(ordered.begin(),ordered.end(),std::pair{task.contig,start},
                                [&](auto const & location, auto index) {
                                    auto const & region = batch[qi].per_bin_regions[index];
                                    return location < std::pair{region.contig,region.start};
                                });
                            if (it != ordered.begin()) {
                                auto candidate = *--it;
                                auto const & region = batch[qi].per_bin_regions[candidate];
                                if (region.contig == task.contig && start >= region.start && start < region.end) ri = candidate;
                            }
                        if (ri != std::numeric_limits<std::size_t>::max()) {
                            auto wi = key_owners[first].second;
                            uint64_t marked;
#ifdef _OPENMP
#pragma omp atomic read
#endif
                            marked = masks[qi][ri][wi/64];
                            // Equal keys have equal window-owner sets, also
                            // for reverse complements. A thread marking the
                            // first window completes that whole group before
                            // the barrier; repeated hits need not revisit it.
                            if (!(marked & (uint64_t{1} << (wi%64))))
                            for (std::size_t owner = first; owner < after; ++owner) {
                            wi = key_owners[owner].second;
                            auto & word = masks[qi][ri][wi/64];
                            auto bit = uint64_t{1} << (wi%64);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                            word |= bit;
                            }
                        }
                        first = after;
                    }
                    if (seen_keys.size() < cache_limit) { seen_keys.insert(w->forward); ++local_insertions; }
                }
            }
            cache_hits.fetch_add(local_cache_hits,std::memory_order_relaxed);
            cache_insertions.fetch_add(local_insertions,std::memory_order_relaxed);
            prefilter_rejections.fetch_add(local_rejections,std::memory_order_relaxed);
            full_lookups.fetch_add(local_lookups,std::memory_order_relaxed);
        });
        production_log.progress("verify","completed_intervals=" + std::to_string(base+count) + "/" + std::to_string(tasks.size()));
    }
    result.verification_prefilter_rejections = prefilter_rejections.load();
    result.verifier_full_key_lookups = full_lookups.load();
    result.candidate_cache_hits = cache_hits.load();
    result.candidate_cache_insertions = cache_insertions.load();
    for (std::size_t qi = 0; qi < batch.size(); ++qi) {
        auto & united = result.candidate_window_masks[qi];
        if (!o.window_output.empty()) result.candidate_region_matches.emplace_back(batch[qi].windows.size(),0);
        for (auto const & mask : masks[qi]) {
            std::size_t count = 0;
            for (std::size_t word = 0; word < mask.size(); ++word) {
                count += std::popcount(mask[word]); united[word] |= mask[word];
                if (!o.window_output.empty()) {
                    auto bits = mask[word];
                    while (bits) { ++result.candidate_region_matches[qi][word*64+std::countr_zero(bits)]; bits &= bits-1; }
                }
            }
            result.per_bin_region_exact_windows[qi].push_back(count);
        }
        for (auto word : united) result.per_bin_candidate_union_windows[qi] += std::popcount(word);
    }
}
