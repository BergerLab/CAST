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

static FILE *open_db_file(char *path, char *fopen_mode);

static char *path_join(char *a, char *b);

static char *basename(char *path);

/*@param dir: The name of the directory to store the database in
 *@param seed_size: The size of the k-mers in the database
 *
 *Initializes a new CaBLAST database and its coarse and compressed databases.
 */
struct cb_database *cb_database_init(char *dir, int32_t seed_size){
    struct cb_database *db;
    struct stat buf;
    FILE *ffasta, *flinks, *fcompressed, *findex_coarse_links,
         *findex_coarse_links_base, *findex_coarse_links_count,
         *findex_coarse_fasta, *findex_coarse_fasta_base, *findex_compressed,
         *findex_params;

    //Get the file paths for the new database's file
    char *pfasta                    = path_join(dir, CABLAST_COARSE_FASTA),
         *plinks                    = path_join(dir, CABLAST_COARSE_LINKS),
         *pindex_coarse_links       =
           path_join(dir, CABLAST_COARSE_LINKS_INDEX),
         *pindex_coarse_links_base  =
           path_join(dir, CABLAST_COARSE_LINKS_BASE_INDEX),
         *pindex_coarse_links_count =
           path_join(dir, CABLAST_COARSE_LINKS_COUNT_INDEX),
         *pindex_coarse_fasta       =
           path_join(dir, CABLAST_COARSE_FASTA_INDEX),
         *pindex_coarse_fasta_base  =
           path_join(dir, CABLAST_COARSE_FASTA_BASE_INDEX),
         *pcompressed               = path_join(dir, CABLAST_COMPRESSED),
         *pindex_compressed         = path_join(dir, CABLAST_COMPRESSED_INDEX),
         *pindex_params             = path_join(dir, CABLAST_PARAMS);

    //Make sure `dir` does not exist.
    if (0 == stat(dir, &buf)) {
        fprintf(stderr,
                "The directory '%s' already exists. A new compressed "
                "database cannot be created in the same directory as an "
                "existing database.\n", dir);
        exit(1); 

        //Just for testing purposes.
        unlink(pfasta);
        unlink(plinks);
        unlink(pindex_coarse_links);
        unlink(pindex_coarse_links_base);
        unlink(pindex_coarse_links_count);
        unlink(pindex_coarse_fasta);
        unlink(pindex_coarse_fasta_base);
        unlink(pcompressed);
        unlink(pindex_compressed);
        unlink(pindex_params);
        rmdir(dir);
    }

    if (0 != mkdir(dir, 0777)) {
        fprintf(stderr, "cb_database_init: 'mkdir %s' failed because: %s\n",
            dir, strerror(errno));
        exit(1);
    }

    db = malloc(sizeof(*db));
    assert(db);

    db->name = basename(dir);

    /*Open the files for the database for binary output (except for the coarse
      FASTA file, which prints ASCII output).*/
    ffasta                    = open_db_file(pfasta, "r+");
    flinks                    = open_db_file(plinks, "rb+");
    findex_coarse_links       = open_db_file(pindex_coarse_links, "rb+");
    findex_coarse_links_base  = open_db_file(pindex_coarse_links_base, "rb+");
    findex_coarse_links_count = open_db_file(pindex_coarse_links_count, "rb+");
    findex_coarse_fasta       = open_db_file(pindex_coarse_fasta, "rb+");
    findex_coarse_fasta_base  = open_db_file(pindex_coarse_fasta_base, "rb+");
    fcompressed               = open_db_file(pcompressed, "rb+");
    findex_compressed         = open_db_file(pindex_compressed, "rb+");
    findex_params             = open_db_file(pindex_params, "rb+");

    //Initialize the coarse and compressed databases
    db->coarse_db = cb_coarse_init(seed_size, ffasta, flinks,
                                 findex_coarse_links, findex_coarse_links_base,
                                 findex_coarse_links_count, findex_coarse_fasta,
                                 findex_coarse_fasta_base, findex_params,
                                 false);
    db->com_db    = cb_compressed_init(fcompressed, findex_compressed, false);

    free(pfasta);
    free(plinks);
    free(pindex_coarse_links);
    free(pindex_coarse_links_base);
    free(pindex_coarse_links_count);
    free(pindex_coarse_fasta);
    free(pindex_coarse_fasta_base);
    free(pcompressed);
    free(pindex_compressed);
    free(pindex_params);

    return db;
}

/*@param dir: The name of the directory where the database is stored
 *@param seed_size: The size of the k-mers in the database
 *@param load_coarse_residues: Whether or not to load all of the residues in
 *  the database's coarse FASTA file into memory.
 *@param load_coarse_links: Whether or not to load all of the links in the
 *  database's coarse links file into memory.
 *@param load_compressed_db: Whether or not to load the compressed database
 *  into memory.
 *@param link_block_size: The number of residues of the coarse FASTA file
 *  each "block" of links represents.
 *
 *Loads a CaBLAST database for reading
 */
struct cb_database_r *
cb_database_r_init(char *dir, int32_t seed_size,
                   bool load_coarse_residues, bool load_coarse_links,
                   bool load_compressed_db, int32_t link_block_size){
    struct cb_database_r *db;
    struct stat buf;
    FILE *ffasta, *flinks, *fcompressed, *findex_coarse_links,
         *findex_coarse_links_base, *findex_coarse_links_count,
         *findex_coarse_fasta, *findex_coarse_fasta_base, *findex_compressed,
         *findex_params;

    char *pfasta                    = path_join(dir, CABLAST_COARSE_FASTA),
         *plinks                    = path_join(dir, CABLAST_COARSE_LINKS),
         *pindex_coarse_links       =
           path_join(dir, CABLAST_COARSE_LINKS_INDEX),
         *pindex_coarse_links_base  =
           path_join(dir, CABLAST_COARSE_LINKS_BASE_INDEX),
         *pindex_coarse_links_count =
           path_join(dir, CABLAST_COARSE_LINKS_COUNT_INDEX),
         *pindex_coarse_fasta       =
           path_join(dir, CABLAST_COARSE_FASTA_INDEX),
         *pindex_coarse_fasta_base  =
           path_join(dir, CABLAST_COARSE_FASTA_BASE_INDEX),
         *pcompressed               = path_join(dir, CABLAST_COMPRESSED),
         *pindex_compressed         = path_join(dir, CABLAST_COMPRESSED_INDEX),
         *pindex_params             = path_join(dir, CABLAST_PARAMS);

    //Make sure the database directory exists.
    if (0 != stat(dir, &buf)) {
        fprintf(stderr, "Could not open '%s' database for reading.", dir);
        exit(1);
    }

    db = malloc(sizeof(*db));
    assert(db);

    db->name = basename(dir);

    ffasta                    = open_db_file(pfasta, "r");
    flinks                    = open_db_file(plinks, "rb");
    findex_coarse_links       = open_db_file(pindex_coarse_links, "rb");
    findex_coarse_links_base  = open_db_file(pindex_coarse_links_base, "rb");
    findex_coarse_links_count = open_db_file(pindex_coarse_links_count, "rb");
    findex_coarse_fasta       = open_db_file(pindex_coarse_fasta, "rb");
    findex_coarse_fasta_base  = open_db_file(pindex_coarse_fasta_base, "rb");
    fcompressed               = open_db_file(pcompressed, "rb");
    findex_compressed         = open_db_file(pindex_compressed, "rb");
    findex_params             = open_db_file(pindex_params, "rb");

    db->coarse_db = cb_coarse_r_init(seed_size, ffasta, flinks,
                                     findex_coarse_links,
                                     findex_coarse_links_base,
                                     findex_coarse_links_count,
                                     findex_coarse_fasta,
                                     findex_coarse_fasta_base,
                                     findex_params, load_coarse_residues,
                                     load_coarse_links, link_block_size);
    db->com_db    = cb_compressed_init(fcompressed, findex_compressed,
                                       load_compressed_db);

    return db;
}

/*Freeing function for a cb_database.  Frees the database and its coarse and
  compressed databases and closes its files.*/
void cb_database_free(struct cb_database *db){
    //All files opened in cb_database_init are closed in subsequent frees.
    cb_coarse_free(db->coarse_db);
    cb_compressed_free(db->com_db);
    free(db->name);
    free(db);
}

/*Freeing function for a cb_database_r.  Frees the database and its coarse and
  compressed databases and closes its files.*/
void cb_database_r_free(struct cb_database_r *db){
    //All files opened in cb_database_init are closed in subsequent frees.
    cb_coarse_r_free(db->coarse_db);
    cb_compressed_free(db->com_db);
    free(db->name);
    free(db);
}

/*Takes in a file path and the mode for opening the file and either opens the
 *file at that path or creates and opens the file if it does not exist.  The
 *program exits if the file cannot be created or opened.
 */
static FILE *open_db_file(char *path, char *fopen_mode){
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

//Joins two strings to create a new file path string.
static char *path_join(char *a, char *b){
    char *joined;

    joined = malloc((1+strlen(a)+1+strlen(b))*sizeof(*joined));
    assert(joined);

    sprintf(joined, "%s/%s", a, b);
    return joined;
}

//Takes in a file path and returns a new string containing the basename.
static char *basename(char *path){
    char *base;
    int i;
    int len;

    len = strlen(path);
    for (i = len; i >= 0 && path[i] != '/'; i--);
    if (i > 0)
        i++;

    base = malloc((1+len-i)*sizeof(*base));
    assert(base);

    strncpy(base, path + i, len - i);

    return base;
}
