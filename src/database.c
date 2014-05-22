#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "coarse.h"
#include "compressed.h"
#include "database.h"
#include "fasta.h"

static FILE * open_db_file(char *path, char *fopen_mode);

static char * path_join(char *a, char *b);

static char * basename(char *path);

struct cb_database *
cb_database_init(char *dir, int32_t seed_size, bool add)
{
    struct cb_database *db;
    struct stat buf;
    FILE *ffasta, *fseeds, *flinks, *fcompressed, *findex_coarse_links,
         *findex_coarse_fasta, *findex_compressed, *findex_params;
    char *pfasta, *pseeds, *plinks, *pcompressed, *pindex_coarse_links,
         *pindex_coarse_fasta, *pindex_compressed, *pindex_params;

    pfasta = path_join(dir, CABLAST_COARSE_FASTA);
    pseeds = path_join(dir, CABLAST_COARSE_SEEDS);
    plinks = path_join(dir, CABLAST_COARSE_LINKS);
    pindex_coarse_links = path_join(dir, CABLAST_COARSE_LINKS_INDEX);
    pindex_coarse_fasta = path_join(dir, CABLAST_COARSE_FASTA_INDEX);
    pcompressed = path_join(dir, CABLAST_COMPRESSED);
    pindex_compressed = path_join(dir, CABLAST_COMPRESSED_INDEX);
    pindex_params = path_join(dir, CABLAST_PARAMS);

    /* If we're not adding to a database, make sure `dir` does not exist. */
    if (!add && 0 == stat(dir, &buf)) {
        /* fprintf(stderr, */
            /* "The directory '%s' already exists. A new compressed " */
            /* "database cannot be created in the same directory as an " */
            /* "existing database. If you want to append to an existing " */
            /* "database, use the '--append' flag.\n", dir); */
        /* exit(1); */

        /* Just for testing purposes. */
        unlink(pfasta);
        unlink(pseeds);
        unlink(plinks);
        unlink(pindex_coarse_links);
        unlink(pindex_coarse_fasta);
        unlink(pcompressed);
        unlink(pindex_compressed);
        rmdir(dir);
    }
    /* Otherwise, check to make sure it *does* exist. */
    if (add && 0 != stat(dir, &buf)) {
        fprintf(stderr, "Could not open '%s' database for appending.", dir);
        exit(1);
    }

    if (0 != mkdir(dir, 0777)) {
        fprintf(stderr, "cb_database_init: 'mkdir %s' failed because: %s\n",
            dir, strerror(errno));
        exit(1);
    }

    db = malloc(sizeof(*db));
    assert(db);

    db->name = basename(dir);

    ffasta = open_db_file(pfasta, "r+");
    fseeds = open_db_file(pseeds, "r+");
    flinks = open_db_file(plinks, "r+");
    findex_coarse_links = open_db_file(pindex_coarse_links, "r+");
    findex_coarse_fasta = open_db_file(pindex_coarse_fasta, "r+");
    fcompressed = open_db_file(pcompressed, "r+");
    findex_compressed = open_db_file(pindex_compressed, "r+");
    findex_params = open_db_file(pindex_params, "r+");

    db->coarse_db = cb_coarse_init(seed_size, ffasta, fseeds, flinks,
                                    findex_coarse_links, findex_coarse_fasta,
                                    findex_params);
    db->com_db = cb_compressed_init(fcompressed, findex_compressed);

    free(pfasta);
    free(pseeds);
    free(plinks);
    free(pindex_coarse_links);
    free(pindex_coarse_fasta);
    free(pcompressed);
    free(pindex_compressed);
    free(pindex_params);

    return db;
}

struct cb_database *
cb_database_read(char *dir, int32_t seed_size)
{
    struct cb_database *db;
    struct stat buf;
    FILE *ffasta, *fseeds, *flinks, *fcompressed, *findex_coarse_links,
         *findex_coarse_fasta, *findex_compressed, *findex_params;
    char *pfasta              = path_join(dir, CABLAST_COARSE_FASTA);
    char *pseeds              = path_join(dir, CABLAST_COARSE_SEEDS);
    char *plinks              = path_join(dir, CABLAST_COARSE_LINKS);
    char *pindex_coarse_links = path_join(dir, CABLAST_COARSE_LINKS_INDEX);
    char *pindex_coarse_fasta = path_join(dir, CABLAST_COARSE_FASTA_INDEX);
    char *pcompressed         = path_join(dir, CABLAST_COMPRESSED);
    char *pindex_compressed   = path_join(dir, CABLAST_COMPRESSED_INDEX);
    char *pindex_params       = path_join(dir, CABLAST_PARAMS);

    /* Make sure the database directory exists. */
    if (0 != stat(dir, &buf)) {
        fprintf(stderr, "Could not open '%s' database for reading.", dir);
        exit(1);
    }

    db = malloc(sizeof(*db));
    assert(db);

    db->name = basename(dir);

    ffasta              = open_db_file(pfasta, "r");
    fseeds              = open_db_file(pseeds, "r");
    flinks              = open_db_file(plinks, "r");
    findex_coarse_links = open_db_file(pindex_coarse_links, "r");
    findex_coarse_fasta = open_db_file(pindex_coarse_fasta, "r");
    fcompressed         = open_db_file(pcompressed, "r");
    findex_compressed   = open_db_file(pindex_compressed, "r");
    findex_params       = open_db_file(pindex_params, "r");

    db->coarse_db = cb_coarse_init(seed_size, ffasta, fseeds, flinks,
                                    findex_coarse_links, findex_coarse_fasta,
                                    findex_params);
    db->com_db = cb_compressed_init(fcompressed, findex_compressed);
    /*cb_database_populate(db, pfasta, plinks);*/

    return db;
}

void cb_database_populate(struct cb_database *db, const char *pfasta,
                           const char *plinks){
    struct fasta_seq_gen *fsg = fasta_generator_start(pfasta, "", 100);
    FILE *flinks = fopen(plinks, "r");
    /*read_coarse(db->coarse_db, flinks, fsg);*/
    fclose(flinks);
}

void
cb_database_free(struct cb_database *db)
{
    /* All files opened in cb_database_init are close in subsequent frees. */
    cb_coarse_free(db->coarse_db);
    cb_compressed_free(db->com_db);
    free(db->name);
    free(db);
}

static FILE *
open_db_file(char *path, char *fopen_mode)
{
    struct stat buf;
    FILE *fp;

    if (0 != stat(path, &buf))
        if (-1 == creat(path, 0666)) {
            fprintf(stderr, "open_db_file: 'creat %s' failed: %s\n",
                path, strerror(errno));
            exit(1);
        }
    if (NULL == (fp = fopen(path, fopen_mode))) {
        fprintf(stderr, "open_db_file: 'fopen %s' failed: %s\n",
            path, strerror(errno));
        exit(1);
    }

    return fp;
}

static char *
path_join(char *a, char *b)
{
    char *joined;

    joined = malloc((1 + strlen(a) + 1 + strlen(b)) * sizeof(*joined));
    assert(joined);

    sprintf(joined, "%s/%s", a, b);
    return joined;
}

static char *
basename(char *path)
{
    char *base;
    int i;
    int len;

    len = strlen(path);
    for (i = len; i >= 0 && path[i] != '/'; i--);
    if (i > 0)
        i++;

    base = malloc((1 + len - i) * sizeof(*base));
    assert(base);

    strncpy(base, path + i, len - i);

    return base;
}
