#pragma once

int export_kmers(Options const & o) {
    std::signal(SIGTERM,fraglncrna::request_stop); std::signal(SIGINT,fraglncrna::request_stop);
    auto bed = load_repeat_bed(o.query_repeat_bed,o.max_repeat_intervals);
    WindowReport report{o}; std::unordered_set<std::string> seen;
    std::string id,sequence; std::size_t index = 0;
    auto flush = [&] {
        if (id.empty()) return;
        auto intervals = bed.intervals.find(id);
        if (intervals != bed.intervals.end()) for (auto const & interval : intervals->second)
            if (interval.end > sequence.size()) throw std::runtime_error("Repeat BED interval outside query: " + id);
        std::vector<Query> batch;
        batch.push_back(make_query(id,index++,sequence,o,nullptr,bed));
        report.write("NA",batch,ExactResult{}); sequence.clear();
    };
    scan_fasta(o.queries,[&](std::string const & next) {
        flush(); if (!seen.insert(next).second) throw std::runtime_error("Duplicate query FASTA identifier: " + next); id = next;
    },[&](std::size_t,char base) {
        if (sequence.size() >= o.max_query_bases) throw std::runtime_error("Query record exceeds --max-query-bases: " + id);
        sequence += base;
    });
    flush();
    for (auto const & [query,ignored] : bed.intervals)
        if (!seen.contains(query)) throw std::runtime_error("Repeat BED query absent from FASTA: " + query);
    report.finish();
    std::cerr << "Exported " << report.rows() << " original windows from " << index << " queries; genomic predicates unavailable (NA).\n";
    return 0;
}
