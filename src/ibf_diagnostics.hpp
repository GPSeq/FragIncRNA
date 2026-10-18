#pragma once


double diagnose_ibf(Options & o, ReferenceInfo const & info) {
    using fraglncrna::production_log;
    production_log.event("ibf-diagnostics", "stage_started", "IBF SETTINGS DIAGNOSTICS scoring=" + o.scoring_method);
    if (!info.bin_count) throw std::runtime_error("Reference has zero bins; provide a nonempty FASTA");
    double probes = o.strand == "both" ? 2.0 : 1.0;
    double proposed = std::min(o.requested_per_bin_fpr, o.genome_fpr / (double(info.bin_count)*probes));
    bool legacy = o.legacy_global_budget || (o.genome_fpr_given && !o.auto_tune_given);
    bool apply = o.ibf_auto_tune == "apply" && o.scoring_method == "global_ibf";
    if (legacy || apply) {
        if (!o.load_index.empty() && apply && !legacy && proposed != o.fpr)
            production_log.event("ibf-diagnostics", "warning", "Loaded archive will not be resized. Rebuild using --save-index NEW.ibf --ibf-auto-tune apply", 3);
        o.fpr = proposed;
    }
    std::size_t bits = planned_bin_bits(o.fragment_size, o.hash_functions, o.fpr);
    std::size_t bytes = planned_index_bytes(info.bin_count, bits);
    double union_estimate = fraglncrna::union_independent(info.bin_count,o.fpr);
    double bound = fraglncrna::union_bound(info.bin_count,o.fpr);
    bool global = o.scoring_method == "global_ibf";
    double window_bound = std::min(1.0,double(info.bin_count)*probes*o.fpr);
    bool unsafe = global && window_bound > o.genome_fpr * (1 + 1e-12);
    std::string recommendation = global ?
        "Use --ibf-auto-tune apply --genome-fpr " + std::to_string(o.genome_fpr) +
        "; compare --hash-functions 5 memory estimate before rebuilding with --save-index NEW.ibf" :
        "Keep per-bin FPR and fragment_size unchanged; validate candidate recall on disjoint controls and inspect false candidate burden before tuning thresholds";
    std::ostringstream message;
    message << std::setprecision(17) << "scoring=" << o.scoring_method << " B=" << info.bin_count
            << " requested_p_bin=" << o.requested_per_bin_fpr << " effective_p_bin=" << o.fpr
            << " p_union_independent=" << union_estimate << " union_bound=" << bound
            << " estimated_bytes=" << bytes << " hashes=" << o.hash_functions
            << " changed=" << (o.fpr != o.requested_per_bin_fpr) << " recommendation=" << recommendation;
    production_log.event("ibf-diagnostics", "stage_progress", message.str());
    if (union_estimate > 0.1)
        production_log.event("ibf-diagnostics", "warning",
            "nominal any-bin false-positive probability near " + std::to_string(union_estimate) +
            " (independent-bin approximation, not certainty or biological FDR). " + message.str(), 3);
    fs::create_directories(o.ibf_tuning_report.parent_path().empty() ? "." : o.ibf_tuning_report.parent_path());
    std::ofstream report(o.ibf_tuning_report);
    if (!report) throw std::runtime_error("Cannot write IBF tuning report");
    report << std::setprecision(17) << "{\n\"scoring_method\":" << json_string(o.scoring_method)
           << ",\n\"auto_tune\":" << json_string(o.ibf_auto_tune)
           << ",\n\"legacy_budget_rule\":" << (legacy ? "true" : "false")
           << ",\n\"changed\":" << (o.fpr != o.requested_per_bin_fpr ? "true" : "false")
           << ",\n\"bin_count\":" << info.bin_count << ",\n\"contig_count\":" << info.ids.size()
           << ",\n\"fragment_size\":" << o.fragment_size << ",\n\"requested_per_bin_fpr\":" << o.requested_per_bin_fpr
           << ",\n\"effective_per_bin_fpr\":" << o.fpr << ",\n\"genome_target\":" << o.genome_fpr
           << ",\n\"p_union_independent\":" << union_estimate << ",\n\"union_bound\":" << bound
           << ",\n\"maximum_lookup_orientations\":" << probes << ",\n\"query_window_union_bound\":" << window_bound
           << ",\n\"bin_bits\":" << bits << ",\n\"estimated_index_bytes\":" << bytes
           << ",\n\"max_index_memory_gb\":" << o.ibf_max_index_memory_gb
           << ",\n\"unsafe_global_settings\":" << (unsafe ? "true" : "false")
           << ",\n\"formula\":\"-expm1(B*log1p(-p_bin)); union bound min(1,B*p_bin)\""
           << ",\n\"assumptions\":\"B-bin formulas describe one indexed key; both-strand windows use up to two probes. Independence is approximate; union bounds require valid individual FPR bounds; not transcript error or biological FDR\""
           << ",\n\"load_estimator\":\"conservative fragment_size insertion upper bound; distinct counts measured during build\""
           << ",\n\"recommendation\":" << json_string(recommendation) << ",\n\"hash_sensitivity\":[";
    auto estimate = [](std::size_t bins, std::size_t fragment, std::size_t hashes, double fpr) {
        try { return std::pair{std::to_string(planned_index_bytes(bins,planned_bin_bits(fragment,hashes,fpr))),std::string{}}; }
        catch (std::runtime_error const & error) { return std::pair{std::string{"null"},std::string{error.what()}}; }
    };
    for (std::size_t h = 1; h <= 5; ++h) {
        if (h != 1) report << ',';
        auto current = estimate(info.bin_count,o.fragment_size,h,o.fpr);
        auto budgeted = estimate(info.bin_count,o.fragment_size,h,proposed);
        report << "{\"hashes\":" << h << ",\"same_p_bin_bytes\":" << current.first
               << ",\"same_p_bin_error\":" << json_string(current.second)
               << ",\"budgeted_p_bin\":" << proposed << ",\"budgeted_bytes\":" << budgeted.first
               << ",\"budgeted_error\":" << json_string(budgeted.second) << '}';
    }
    report << "],\n\"fragment_alternatives\":[";
    for (std::size_t scale = 1; scale <= 4; scale *= 2) {
        if (scale != 1) report << ',';
        if (o.fragment_size > std::numeric_limits<std::size_t>::max() / scale) {
            report << "{\"scale\":" << scale << ",\"fragment_size\":null,\"estimated_bytes\":null,\"error\":\"fragment size overflow\",\"applied\":false}";
            continue;
        }
        auto fragment = o.fragment_size * scale;
        std::size_t bins = 0;
        for (auto length : info.lengths) bins += length / fragment + (length % fragment != 0);
        auto alternative = estimate(bins,fragment,o.hash_functions,o.fpr);
        report << "{\"fragment_size\":" << fragment << ",\"bins\":" << bins
               << ",\"estimated_bytes\":" << alternative.first << ",\"error\":" << json_string(alternative.second)
               << ",\"applied\":false}";
    }
    report << "]\n}\n"; report.close();
    if (!report) throw std::runtime_error("Failed writing IBF tuning diagnostics");
    if (o.ibf_max_index_memory_gb > 0 && double(bytes) > o.ibf_max_index_memory_gb * 1e9)
        throw std::runtime_error("IBF estimated bytes " + std::to_string(bytes) + " exceed --ibf-max-index-memory-gb; increase budget or inspect hash sensitivity in " + o.ibf_tuning_report.string());
    if (unsafe && o.ibf_fail_on_unsafe_settings)
        throw std::runtime_error("Unsafe global_ibf settings; " + recommendation);
    if (!o.save_index.empty()) {
        auto parent = o.save_index.parent_path().empty() ? fs::path{"."} : o.save_index.parent_path();
        fs::create_directories(parent);
        if (fs::space(parent).available < bytes)
            throw std::runtime_error("Insufficient disk for saved IBF; choose a larger destination or smaller design");
    }
    production_log.event("ibf-diagnostics", "stage_completed", "tuning_report=" + o.ibf_tuning_report.string());
    return union_estimate;
}
