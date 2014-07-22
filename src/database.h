#ifndef __CABLAST_DATABASE_H__
#define __CABLAST_DATABASE_H__

#include <stdbool.h>

#include "coarse.h"
#include "compressed.h"

#define CABLAST_COARSE_FASTA "coarse.fasta"
#define CABLAST_COARSE_LINKS "coarse.links"
#define CABLAST_COARSE_LINKS_INDEX "coarse.links.index"
#define CABLAST_COARSE_LINKS_COUNT_INDEX "coarse.links.count.index"
#define CABLAST_COARSE_LINKS_BASE_INDEX "coarse.links.base.index"
#define CABLAST_COARSE_FASTA_INDEX "coarse.fasta.index"
#define CABLAST_COARSE_FASTA_BASE_INDEX "coarse.fasta.base.index"
#define CABLAST_COMPRESSED "compressed.cb"
#define CABLAST_COMPRESSED_INDEX "compressed.cb.index"
#define CABLAST_PARAMS "params"

struct cb_database {
    char *name;
    struct cb_coarse *coarse_db;
    struct cb_compressed *com_db;
};

struct cb_database_r {
    char *name;
    struct cb_coarse_r *coarse_db;
    struct cb_compressed *com_db;
};

struct cb_database *cb_database_init(char *dir, int32_t seed_size);

struct cb_database_r *
cb_database_r_init(char *dir, int32_t seed_size,
                   bool load_coarse_residues, bool load_coarse_links,
                   bool load_compressed_db, int32_t link_block_size);

void cb_database_free(struct cb_database *db);

void cb_database_r_free(struct cb_database_r *db);

#endif
