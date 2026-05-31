#include "config.hpp"
#include "fragmenter.hpp"
#include "logger.hpp"
#include "query_processor.hpp"
#include "reference_index.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>

namespace fs = std::filesystem;
using namespace seqan3::literals;

namespace
{

std::string to_string(seqan3::dna5_vector const & seq)
{
    std::string out;
    out.reserve(seq.size());

    for (auto const base : seq)
        out.push_back(seqan3::to_char(base));

    return out;
}

fs::path make_temp_dir(std::string const & name)
{
    auto dir = fs::temp_directory_path() / name;
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

void write_fasta(fs::path const & path,
                 std::string const & id,
                 std::string const & sequence)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Failed to create FASTA file: " + path.string());

    out << '>' << id << '\n' << sequence << '\n';
}

void write_text(fs::path const & path, std::string const & text)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Failed to create text file: " + path.string());

    out << text;
}

void expect(bool condition, std::string const & message)
{
    if (!condition)
        throw std::runtime_error(message);
}

std::string read_text(fs::path const & path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Failed to read file: " + path.string());

    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

void test_fragmenter_returns_overlapping_fragments()
{
    auto temp_dir = make_temp_dir("lncrna_mers_test_fragments");
    auto fasta = temp_dir / "ref.fa";
    write_fasta(fasta, "ref1", "ACGTACGTAA");

    Config cfg;
    cfg.fragment_size = 5;

    FragmenterDna5 fragmenter{cfg};
    auto fragments = fragmenter.fragment_reference(fasta, "ref1");

    expect(fragments.size() == 5, "expected 5 fragments");
    expect(to_string(fragments[0]) == "ACGTA", "unexpected fragment 0");
    expect(to_string(fragments[1]) == "GTACG", "unexpected fragment 1");
    expect(to_string(fragments[2]) == "ACGTA", "unexpected fragment 2");
    expect(to_string(fragments[3]) == "GTAA", "unexpected fragment 3");
    expect(to_string(fragments[4]) == "AA", "unexpected fragment 4");
}

void test_fragmenter_writes_fragment_fasta_when_enabled()
{
    auto temp_dir = make_temp_dir("lncrna_mers_test_fragment_output");
    auto fasta = temp_dir / "ref.fa";
    auto out_dir = temp_dir / "out";
    write_fasta(fasta, "ref2", "AACCGGTT");

    Config cfg;
    cfg.fragment_size = 4;
    cfg.store_fragments = true;
    cfg.output_dir = out_dir;

    FragmenterDna5 fragmenter{cfg};
    auto fragments = fragmenter.fragment_reference(fasta, "ref2");

    auto fragment_file = out_dir / "ref2_fragments.fasta";
    expect(fs::exists(fragment_file), "fragment FASTA was not created");
    expect(fragments.size() == 8, "expected 8 stored fragments");

    std::ifstream in(fragment_file);
    std::string contents((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());

    expect(contents.find(">ref2_frag0") != std::string::npos, "missing first fragment ID");
    expect(contents.find(">ref2_frag7") != std::string::npos, "missing last fragment ID");
    //expect(contents.find("> ref2_frag7") != std::string::npos, "missing last fragment ID");
}

void test_fragmenter_rejects_too_small_fragment_size()
{
    auto temp_dir = make_temp_dir("lncrna_mers_test_fragment_validation");
    auto fasta = temp_dir / "ref.fa";
    write_fasta(fasta, "ref3", "ACGT");

    Config cfg;
    cfg.fragment_size = 3;

    FragmenterDna5 fragmenter{cfg};

    bool threw = false;
    try
    {
        (void) fragmenter.fragment_reference(fasta, "ref3");
    }
    catch (std::runtime_error const & ex)
    {
        threw = std::string{ex.what()}.find("fragment_size must be >= 4") != std::string::npos;
    }

    expect(threw, "expected fragment_size validation error");
}

void test_reference_index_rejects_invalid_kmer_size()
{
    Config cfg;
    cfg.kmer_size = 6;
    cfg.ibf.hash_functions = 2;
    cfg.ibf.fpr = 0.01;

    std::vector<seqan3::dna5_vector> fragments{seqan3::dna5_vector{"ACGT"_dna5}};

    bool threw = false;
    try
    {
        (void) ReferenceIndexDna5{"ref4", fragments, cfg};
    }
    catch (std::runtime_error const & ex)
    {
        threw = std::string{ex.what()}.find("kmer_size must be > 0 and <=") != std::string::npos;
    }

    expect(threw, "expected kmer_size validation error");
}

void test_reference_index_reports_fragment_count_as_bin_count()
{
    auto temp_dir = make_temp_dir("lncrna_mers_test_ibf");

    Config cfg;
    cfg.kmer_size = 3;
    cfg.ibf.hash_functions = 2;
    cfg.ibf.fpr = 0.01;
    cfg.output_dir = temp_dir;

    std::vector<seqan3::dna5_vector> fragments{
        seqan3::dna5_vector{"ACGTA"_dna5},
        seqan3::dna5_vector{"CGTAC"_dna5},
        seqan3::dna5_vector{"GTACC"_dna5}
    };

    ReferenceIndexDna5 index{"ref5", fragments, cfg};
    expect(index.bin_count() == fragments.size(), "unexpected IBF bin count");
}

void test_config_loader_reads_sectioned_toml()
{
    auto temp_dir = make_temp_dir("lncrna_mers_test_config_ok");
    auto config_path = temp_dir / "config.toml";

    write_text(config_path,
               "[general]\n"
               "index_method = \"hibf\"\n"
               "ref_dir = \"/refs\"\n"
               "query_file = \"/queries.fa\"\n"
               "fragment_size = 8000\n"
               "kmer_size = 15\n"
               "hit_threshold = 13\n"
               "threads = 8\n"
               "output_dir = \"./out\"\n"
               "output_file = \"results.tsv\"\n"
               "log_file = \"run.log\"\n"
               "store_fragments = false\n"
               "store_index = true\n"
               "cleanup_index = true\n"
               "single_results_writer = true\n"
               "\n"
               "[ibf]\n"
               "hash_functions = 3\n"
               "fpr = 0.01\n"
               "\n"
               "[hibf]\n"
               "hash_functions = 2\n"
               "maximum_fpr = 0.05\n"
               "relaxed_fpr = 0.30\n"
               "threads = 4\n"
               "sketch_bits = 12\n");

    auto cfg = load_config_from_toml(config_path);

    expect(cfg.ref_dir == "/refs", "unexpected ref_dir");
    expect(cfg.query_file == "/queries.fa", "unexpected query_file");
    expect(cfg.index_method == IndexMethod::hibf, "unexpected index method");
    expect(cfg.fragment_size == 8000, "unexpected fragment_size");
    expect(cfg.kmer_size == 15, "unexpected kmer_size");
    expect(cfg.hit_threshold == 13, "unexpected hit_threshold");
    expect(cfg.threads == 8, "unexpected threads");
    expect(cfg.output_dir == "./out", "unexpected output_dir");
    expect(cfg.output_file == "results.tsv", "unexpected output_file");
    expect(cfg.log_file == "run.log", "unexpected log_file");
    expect(cfg.store_fragments == false, "unexpected store_fragments");
    expect(cfg.store_index == true, "unexpected store_index");
    expect(cfg.cleanup_index == true, "unexpected cleanup_index");
    expect(cfg.single_results_writer == true, "unexpected single_results_writer");
    expect(cfg.ibf.hash_functions == 3, "unexpected ibf.hash_functions");
    expect(cfg.ibf.fpr == 0.01, "unexpected ibf.fpr");
    expect(cfg.hibf.hash_functions == 2, "unexpected hibf.hash_functions");
    expect(cfg.hibf.maximum_fpr == 0.05, "unexpected hibf.maximum_fpr");
    expect(cfg.hibf.relaxed_fpr == 0.30, "unexpected hibf.relaxed_fpr");
    expect(cfg.hibf.threads == 4, "unexpected hibf.threads");
}

void test_config_loader_rejects_missing_required_key()
{
    auto temp_dir = make_temp_dir("lncrna_mers_test_config_missing");
    auto config_path = temp_dir / "config.toml";

    write_text(config_path,
               "query_file = \"/queries.fa\"\n"
               "fragment_size = 8000\n");

    bool threw = false;
    try
    {
        (void) load_config_from_toml(config_path);
    }
    catch (std::runtime_error const & ex)
    {
        threw = std::string{ex.what()}.find("ref_dir") != std::string::npos;
    }

    expect(threw, "expected missing required key error");
}

void test_query_processor_writes_matching_kmer_hit_file()
{
    auto temp_dir = make_temp_dir("lncrna_mers_test_query_processor");
    auto query_fasta = temp_dir / "queries.fa";
    auto results_path = temp_dir / "results_ref6.tsv";

    write_fasta(query_fasta, "query0", "ACGTAC");

    Config cfg;
    cfg.query_file = query_fasta;
    cfg.output_dir = temp_dir;
    cfg.kmer_size = 3;
    cfg.hit_threshold = 1;
    cfg.ibf.hash_functions = 2;
    cfg.ibf.fpr = 0.01;

    std::vector<seqan3::dna5_vector> fragments{
        seqan3::dna5_vector{"ACGTA"_dna5},
        seqan3::dna5_vector{"CGTAC"_dna5}
    };

    ReferenceIndexDna5 index{"ref6", fragments, cfg};
    QueryProcessorDna5 processor{cfg, index, "ref6"};
    processor.run_write_per_ibf(results_path);

    auto kmer_hits = read_text(temp_dir / "ref6_kmers.tsv");
    expect(kmer_hits == "query_index\tkmer_indices\tkmer_counts\n"
                        "0\t1/2/3/4\t1/2/2/1\n",
           "unexpected k-mer hit file contents");

    auto unique_hits = read_text(temp_dir / "unique_mers" / "ref6.tsv");
    expect(unique_hits == "query_index\tibf_unique_kmer\n"
                          "0\tACG\n"
                          "0\tTAC\n",
           "unexpected unique k-mer file contents");
}

} // namespace

int main()
{
    auto log_dir = make_temp_dir("lncrna_mers_test_logs");
    Logger::init((log_dir / "tests.log").string());

    try
    {
        test_fragmenter_returns_overlapping_fragments();
        test_fragmenter_writes_fragment_fasta_when_enabled();
        test_fragmenter_rejects_too_small_fragment_size();
        test_reference_index_rejects_invalid_kmer_size();
        test_reference_index_reports_fragment_count_as_bin_count();
        test_config_loader_reads_sectioned_toml();
        test_config_loader_rejects_missing_required_key();
        test_query_processor_writes_matching_kmer_hit_file();
    }
    catch (std::exception const & ex)
    {
        std::cerr << "Test failure: " << ex.what() << '\n';
        return 1;
    }

    std::cout << "All unit tests passed.\n";
    return 0;
}
