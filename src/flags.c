#include <stdint.h>

#include "flags.h"
#include "util.h"

struct opt_config *
load_compress_args()
{
    struct opt_config *conf;
    int32_t cpus;

    conf = opt_config_init();
    cpus = num_cpus();

    opt_flag_int(conf, 
        &compress_flags.gapped_window_size, "gapped-window-size", 100,
        "The size of the gapped match window.");
    opt_flag_int(conf,
        &compress_flags.ungapped_window_size, "ungapped-window-size", 10,
        "The size of the ungapped match window.");
    opt_flag_int(conf,
        &compress_flags.match_extend, "match-extend", 30,
        "The maximum number of residues to blindly extend a match without\n"
        "\tregard to sequence identity. This is to avoid small sequences\n"
        "\tin the coarse database.");
    opt_flag_int(conf,
        &compress_flags.match_kmer_size, "match-kmer-size", 10,
        "The size of the K-mer fragments to match in ungapped extension.");
    opt_flag_int(conf,
        &compress_flags.match_seq_id_threshold, "match-seq-id-threshold", 70,
        "The sequence identity threshold of an entire match.");
    opt_flag_int(conf,
        &compress_flags.procs, "procs", cpus,
        "The number of total CPUs to use to divide work.");
    opt_flag_int(conf,
        &compress_flags.map_seed_size, "map-seed-size", 10,
        "The size of a seed in the K-mer map. This size combined with "
        "'ext-seed-size' forms the total seed size.");
    opt_flag_int(conf,
        &compress_flags.ext_seed_size, "ext-seed-size", 0,
        "The additional residues to require for each seed match.");
    opt_flag_int(conf,
        &compress_flags.ext_seq_id_threshold, "ext-seq-id-threshold", 50,
        "The sequence identity threshold of [un]gapped extension.");

    opt_flag_int(conf,
        &compress_flags.min_match_len, "min-match-len", 300,
        "The minimum size of an extension needed for a match to be recorded.");
    opt_flag_int(conf,
        &compress_flags.max_kmer_freq, "max-kmer-freq", 500,
        "The maximum number of entries for a k-mer in the seeds table.");
    opt_flag_int(conf,
        &compress_flags.overlap, "overlap", 100,
        "The maximum number of entries for a k-mer in the seeds table.");
    opt_flag_int(conf,
        &compress_flags.max_chunk_size, "max-chunk-size", 10000,
        "The maximum number of bases that are checked before adding a "
        "chunk without a match.");
    opt_flag_int(conf,
        &compress_flags.min_progress, "min-progress", 50,
        "The minimum progress that needs to be made in each direction of "
        "attempt_ext to consider a seed as a potential match.");
    opt_flag_int(conf,
        &compress_flags.consec_match_clump_size, "consec-match-clump-size", 4,
        "The number of consecutive matches needed to check for a bad window in "
        "align_ungapped.");
    opt_flag_int(conf,
        &compress_flags.window_ident_thresh, "window-ident-thresh", 85,
        "The minimum number of matches that must be in the array of matches "
        "for the current window not to be considered a bad window.");
    opt_flag_int(conf,
        &compress_flags.btwn_match_min_dist_check, 
        "btwn-match-min-dist-check", 10,
        "The number of bases scanned since the last clump of matches in "
        "align_ungapped that results in checking for 50%% identity since the "
        "last clump.");
    opt_flag_double(conf,
        &compress_flags.btwn_match_ident_thresh,
        "btwn-match-ident_thresh", ((double)0.5),
        "The identity threshold for continuing align_ungapped if it has been "
        "more than 10 bases since the last clump of matches.");
    opt_flag_int(conf,
        &compress_flags.attempt_ext_len, "attempt-ext-len", 50,
        "The minimum total length of an extension after running attempt_ext in "
        "both directions needed to call extend_match.");

    return conf;
}


struct opt_config *
load_search_args()
{
    struct opt_config *conf;
    int32_t cpus;

    conf = opt_config_init();
    cpus = num_cpus();

    opt_flag_int(conf,
        &search_flags.map_seed_size, "map-seed-size", 10,
        "The size of a seed in the K-mer map. This size combined with "
        "'ext-seed-size' forms the total seed size.");
    opt_flag_string(conf,
        &search_flags.coarse_evalue, "coarse-evalue", "1e-20",
        "The e-value used during coarse search.  To set the e-value for fine "
        "search, add it as an argument in --blast-args, which is the list of "
        "arguments to pass into BLAST during fine search.");
    opt_flag_bool(conf,
        &search_flags.no_cleanup, "no-cleanup",
        "Set to true to keep the coarse search results XML file, the fine "
        "database FASTA file, and the last fine query FASTA file.");
    opt_flag_bool(conf,
        &search_flags.show_hit_info, "show-hit-info",
        "Set to true to output information on each fine BLAST hit.");

    return conf;
}

