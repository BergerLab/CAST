#ifndef __CABLAST_OPT_H__
#define __CABLAST_OPT_H__

#include <stdint.h>
#include <stdbool.h>

#include "clibs/include/opt.h"

struct compress_flags {
    int32_t gapped_window_size;
    int32_t min_match_len;
    int32_t procs;
    int32_t map_seed_size;
    int32_t ext_seq_id_threshold;

    int32_t max_kmer_freq;
    int32_t overlap;
    int32_t max_chunk_size;
    int32_t min_progress;
    int32_t consec_match_clump_size;
    int32_t window_ident_thresh;
    int32_t btwn_match_min_dist_check;
    int32_t notify_every;
    int32_t max_consec_mismatch;
} compress_flags;

struct search_flags {
    char    *coarse_evalue;
    bool    no_cleanup;
    bool    show_hit_info;
    bool    hide_progress;
    bool    load_coarse_residues;
    bool    load_coarse_links;
    bool    load_coarse_db;
    bool    load_compressed_db;
    bool    fine_blast_db;
    int32_t link_block_size;
} search_flags;

struct cablat_flags {
    bool    load_coarse_residues;
    bool    load_coarse_links;
    bool    load_coarse_db;
    bool    load_compressed_db;
    bool    fine_blat_db;
    int32_t link_block_size;
    char    *output_expanded_fasta;
} cablat_flags;


struct opt_config *load_compress_args();
struct opt_config *load_search_args();
struct opt_config *load_cablat_args();

#endif
