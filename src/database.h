#ifndef __CABLAST_DATABASE_H__
#define __CABLAST_DATABASE_H__

#include <stdbool.h>

#include "coarse.h"
#include "compressed.h"

#define CABLAST_COARSE_FASTA "coarse.fasta"
#define CABLAST_COARSE_LINKS "coarse.links"
#define CABLAST_COARSE_SEEDS "coarse.seeds"
#define CABLAST_COARSE_LINKS_INDEX "coarse.links.index"
#define CABLAST_COARSE_FASTA_INDEX "coarse.fasta.index"
#define CABLAST_COMPRESSED "compressed.cb"
#define CABLAST_COMPRESSED_INDEX "compressed.cb.index"
#define CABLAST_PARAMS "params"

struct cb_database {
    char *name;
    struct cb_coarse *coarse_db;
    struct cb_compressed *com_db;
};

struct cb_database *
cb_database_init(char *dir, int32_t seed_size, bool add);

struct cb_database *
cb_database_read(char *dir, int32_t seed_size, bool load_coarse_residues);

void cb_database_populate(struct cb_database *db, const char *pfasta,
                           const char *plinks);

void
cb_database_free(struct cb_database *db);

#endif
