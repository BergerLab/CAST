#include <assert.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "clibs/include/ds.h"
#include "clibs/include/opt.h"

#include "bitpack.h"
#include "coarse.h"
#include "compressed.h"
#include "compression.h"
#include "database.h"
#include "decompression.h"
#include "fasta.h"
#include "flags.h"
#include "progress-bar.h"
#include "seq.h"
#include "util.h"
#include "xml.h"

//Joins two strings to create a new file path string.
static char *path_join(char *a, char *b){
    char *joined;

    joined = malloc((1+strlen(a)+1+strlen(b))*sizeof(*joined));
    assert(joined);

    sprintf(joined, "%s/%s", a, b);
    return joined;
}

char *get_blast_args(struct opt_args *args);

/*Runs BLAST on the coarse database and stores the results in a temporary
  XML file.*/
void blast_coarse(struct opt_args *args, uint64_t dbsize){
    char *blastn,
         *input_path = path_join(args->args[0], CABLAST_COARSE_FASTA),
         *blastn_command =
           "blastn -db  -outfmt 5 -query  -dbsize  -task blastn -evalue  "
           "> CaBLAST_temp_blast_results.xml";
    int command_length = strlen(blastn_command) + strlen(input_path) +
                                                  strlen(args->args[1]) + 31;

    blastn = malloc(command_length*sizeof(*blastn));
    assert(blastn);

    sprintf(blastn,
            "blastn -db %s -outfmt 5 -query %s -dbsize %lu -task blastn"
            " -evalue %s > CaBLAST_temp_blast_results.xml",
            input_path, args->args[1], dbsize, search_flags.coarse_evalue);

    if (!search_flags.hide_progress)
        fprintf(stderr, "%s\n", blastn);

    system(blastn);

    free(input_path);
    free(blastn);
}

//Runs BLAST on the fine FASTA file
void blast_fine(struct opt_args *args, uint64_t dbsize){
    char *blastn,
         *blast_args      = get_blast_args(args),
         *fine_blast_args = is_substring("-evalue", blast_args) ?
                              blast_args : "-evalue 1e-30";

    int command_length = 1024;

    blastn = malloc(command_length*sizeof(*blastn));
    assert(blastn);

    sprintf(blastn,
            "blastn %s CaBLAST_fine.fasta -query %s "
            "-dbsize %lu -task blastn -outfmt 5 %s > "
            "CaBLAST_results.xml",
            search_flags.fine_blast_db ? "-db" : "-subject", args->args[1],
            dbsize, fine_blast_args);

    if (!search_flags.hide_progress)
          fprintf(stderr, "\n%s\n", blastn);

    system(blastn); //Run fine BLAST

    free(blast_args);
    free(blastn);
}

/*A function for traversing a parsed XML tree.  Takes in the root node, a
 *void * accumulator, and a function that takes in an xmlNode and a void *
 *accumulator and traverses the tree, applying the function on each node.
 */
void traverse_blast_xml(xmlNode *root, void (*f)(xmlNode *,void *), void *acc){
    for (; root; root = root->next) {
        f(root, acc);
        traverse_blast_xml(root->children, f, acc);
    }
}

/*A wrapper function for traverse_blast_xml that treats the xmlNode passed into
  r as the root of the XML tree, skipping any sibling nodes r has*/
void traverse_blast_xml_r(xmlNode *r, void (*f)(xmlNode *, void *), void *acc){
    f(r, acc);
    traverse_blast_xml(r->children, f, acc);
}


/*A function for traversing a parsed XML tree.  Takes in the root node, a
 *void * accumulator, a function that takes in an xmlNode and a void *
 *accumulator, and the name of a node and traverses the tree, applying the
 *function on each node but ending the traversal if the name of the current node
 *is the name passed into stop.
 */
void traverse_blast_xml_until(xmlNode *root, void (*f)(xmlNode *,void *),
                              void *acc, char *stop){
    for (; root; root = root->next) {
        f(root, acc);
        if (strcmp((char *)root->name, stop) != 0)
            traverse_blast_xml_until(root->children, f, acc, stop);
    }
}

/*A wrapper function for traverse_blast_xml_until that treats the xmlNode
 *passed into r as the root of the XML tree, skipping any sibling nodes r has.
 *If the name of r is the same as the name passed into stop, traversal does not
 *happen past r.
 */
void traverse_blast_xml_until_r(xmlNode *r, void (*f)(xmlNode *, void *),
                                void *acc, char *stop){
    f(r, acc);
    if (strcmp((char *)r->name, stop) != 0)
        traverse_blast_xml_until(r->children, f, acc, stop);
}


/*Takes in the xmlNode representing a BLAST hsp and populates a struct hsp
  with its data.*/
struct hsp *populate_blast_hsp(xmlNode *node){
    struct hsp *h;
    bool hit_frame_fwd = true;

    h = malloc(sizeof(*h));
    assert(h);

    h->xml_name = "Hsp";
    for (; node; node = node->next){
        if (!strcmp((char *)node->name, "Hsp_num"))
            h->num = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_evalue"))
            h->evalue = atof((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_query-from"))
            h->query_from = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_query-to"))
            h->query_to = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_hit-from"))
            h->hit_from = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_hit-to"))
            h->hit_to = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_hit-frame"))
            hit_frame_fwd = atoi((char *)node->children->content) == 1;
    }
    if (!hit_frame_fwd) {
        int hit_to = h->hit_from;
        h->hit_from = h->hit_to;
        h->hit_to = hit_to;
    }
    return h;
}

/*Takes in the xmlNode representing a BLAST hit and populates a struct hit
  with its data.*/
struct hit *populate_blast_hit(xmlNode *node){
    struct hit *h = malloc(sizeof(*h));
    assert(h);

    h->xml_name = "Hit";
    for (; node; node = node->next) {
        if (!strcmp((char *)node->name, "Hit_num"))
            h->num = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hit_accession")) 
            h->accession = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hit_hsps")){
            struct DSVector *hsps = ds_vector_create();
            xmlNode *hsp_node = node->children;
            for (; hsp_node; hsp_node = hsp_node->next)
                if (!strcmp((char *)hsp_node->name, "Hsp"))
                    ds_vector_append(hsps,
                           (void *)populate_blast_hsp(hsp_node->children));
            h->hsps = hsps;
        }
    }
    return h;
}

/*Takes in an xmlNode and a vector of xmlNodes converted to a void *
  and adds the node to the vector if the node represents a BLAST hit.*/
void add_blast_hit(xmlNode *node, void *hits){
    if (!strcmp((char *)(node->name), "Hit"))
        ds_vector_append((struct DSVector *)hits,
                         (void *)populate_blast_hit(node->children));
}

/*Takes in the xmlNode for the root of a parsed BLAST XML tree and returns
  a vector of all of the nodes in the tree that represent BLAST hits.*/
struct DSVector *get_blast_hits(xmlNode *node){
    struct DSVector *hits = ds_vector_create();
    traverse_blast_xml_until_r(node, add_blast_hit, hits, "Hit");
    return hits;
}

/*Takes in an xmlNode and a vector of xmlNodes converted to a void *
  and adds the node to the vector if the node represents an iteration.*/
void add_blast_iteration(xmlNode *node, void *iterations){
    if (!strcmp((char *)(node->name), "Iteration"))
        ds_vector_append((struct DSVector *)iterations, (void *)node);
}


/*Takes in the xmlNode for the root of a parsed BLAST XML tree and returns
  a vector of all of the nodes in the tree that represent iterations.*/
struct DSVector *get_blast_iterations(xmlNode *node){
    struct DSVector *iterations = ds_vector_create();
    traverse_blast_xml_until_r(node, add_blast_iteration,
                               iterations, "Iteration");
    return iterations;
}

/*Takes in a vector of expansion structs for the expanded BLAST hits from a
  query and outputs them to the FASTA file CaBLAST_fine.fasta.*/
void write_fine_fasta(struct DSVector *oseqs){
    FILE *temp = fopen("CaBLAST_fine.fasta", "w");
    int i;

    if (!temp) {
        fprintf(stderr, "Could not open CaBLAST_fine.fasta for writing\n");
        return;
    }
    for (i = 0; i < oseqs->size; i++) {
        struct cb_seq *current_seq =
          ((struct cb_hit_expansion *)ds_vector_get(oseqs, i))->seq;
        fprintf(temp, "> %s\n%s\n", current_seq->name, current_seq->residues);
    }
    fclose(temp);

    if (search_flags.fine_blast_db)
        system("makeblastdb -dbtype nucl -in CaBLAST_fine.fasta -out "
               "CaBLAST_fine.fasta");
}

/*Takes in the arguments for the program and returns a string of all of the
 *arguments for fine BLAST (the arguments after --blast-args) concatenated and
 *separated by a space.
 */
char *get_blast_args(struct opt_args *args){
    int i = 0, index = -1, length = 1;
    char *blast_args = NULL;

    for (i = 0; i < args->nargs; i++)
        if (index >= 0)
            length += strlen(args->args[i]) + 1;
        else
            index = strcmp(args->args[i], "--blast_args") == 0 ? i : -1;

    blast_args = malloc(length*sizeof(*args));
    assert(blast_args);

    *blast_args = '\0';
    if (index != -1)
        for (i = index + 1; i < args->nargs; i++) {
            blast_args = strcat(blast_args, args->args[i]);
            if (i < args->nargs - 1)
                blast_args = strcat(blast_args, " ");
        }
    return blast_args;
}

/*Takes in as input the vector of iterations from coarse BLAST, an index into
 *the vector (representing the index of which query we want to get the results
 *from), and the database we are using for search and returns a vector of
 *every original sequence section re-created from the calls to
 *cb_coarse_expand for the hits in the iteration we are expanding.
 */
struct DSVector *expand_blast_hits(struct DSVector *iterations, int index,
                                   struct cb_database_r *db){
    struct DSVector *expanded_hits = ds_vector_create(),
                    *hits = get_blast_hits(
                              (xmlNode *)ds_vector_get(iterations, index));
    int i = 0, j = 0, k = 0;

    for (i = 0; i < hits->size; i++) {
        struct hit *current_hit = (struct hit *)ds_vector_get(hits, i);
        struct DSVector *hsps = current_hit->hsps;
        for (j = 0; j < hsps->size; j++) {
            struct DSVector *oseqs;
            struct hsp *h = (struct hsp *)ds_vector_get(hsps, j);
            int32_t coarse_start  = h->hit_from-1, coarse_end = h->hit_to-1,
                    coarse_seq_id = current_hit->accession;

            oseqs = cb_coarse_expand(db->coarse_db, db->com_db, coarse_seq_id,
                                     coarse_start, coarse_end, 50);
            for (k = 0; k < oseqs->size; k++)
                ds_vector_append(expanded_hits, ds_vector_get(oseqs, k));

            ds_vector_free_no_data(oseqs);
        }

        ds_vector_free(current_hit->hsps);
        free(current_hit);
    }
    ds_vector_free_no_data(hits);

    return expanded_hits;
}

int main(int argc, char **argv){
    FILE *query_file = NULL, *test_hits_file = NULL;
    struct cb_database_r *db = NULL;
    struct opt_config *conf;
    struct opt_args *args;
    struct DSVector *iterations = NULL, *expanded_hits = NULL, *queries = NULL,
                    *oseqs = ds_vector_create();
    struct fasta_seq *query = NULL;
    xmlDoc *doc = NULL;
    xmlNode *root = NULL;
    uint64_t dbsize = 0;
    int i = 0, j = 0;

    conf = load_search_args();
    args = opt_config_parse(conf, argc, argv);

    if (args->nargs < 2) {
        fprintf(stderr, 
                "Usage: %s [flags] database-dir fasta-file "
                "[ --blast_args BLASTN_ARGUMENTS ]\n", argv[0]);
        opt_config_print_usage(conf);
        exit(1);
    }

    system("rm CaBLAST_results.xml");

    if (!search_flags.hide_progress)
        fprintf(stderr, "Loading database data\n\n");


    db = cb_database_r_init(args->args[0],
                            (search_flags.load_coarse_db ||
                             search_flags.load_coarse_residues),
                            (search_flags.load_coarse_db ||
                             search_flags.load_coarse_links),
                            search_flags.load_compressed_db,
                            search_flags.link_block_size);
    dbsize = read_int_from_file(8, db->coarse_db->db->file_params);

    if (!search_flags.hide_progress)
        fprintf(stderr, "Running coarse BLAST\n\n");
    blast_coarse(args, dbsize);

    query_file = fopen(args->args[1], "r");
    queries = ds_vector_create();
    query = fasta_read_next(query_file, "");
    while (query) {
        ds_vector_append(queries, (void *)query);
        query = fasta_read_next(query_file, "");
    }

    fclose(query_file);

    if (!search_flags.hide_progress)
        fprintf(stderr, "Processing coarse BLAST hits for fine BLAST\n\n");
    if (search_flags.show_hit_info)
        test_hits_file = fopen("CaBLAST_hits.txt", "w");

    //Parse the XML file generated from coarse BLAST and get its iterations.
    doc = xmlReadFile("CaBLAST_temp_blast_results.xml", NULL, 0);
    if (doc == NULL) {
        fprintf(stderr, "Could not parse CaBLAST_temp_blast_results.xml\n");
        return 0;
    }
    root = xmlDocGetRootElement(doc);
    iterations = get_blast_iterations(root);

    for (i = 0; i < iterations->size; i++) {
        if (!search_flags.hide_progress) {
            int32_t digits_full = floor(log10((double)iterations->size)),
                    digits_i    = floor(log10((double)i)), spaces;
            char *bar = progress_bar(i, iterations->size);
            spaces = digits_full - digits_i;
            fprintf(stderr, "\r");
            fprintf(stderr, "iteration: %d/%d", i+1, iterations->size);
            for (j = 0; j < spaces; j++)
                putc(' ', stderr);
            fprintf(stderr, " %s ", bar);
            free(bar);
        }

        /*Expand any BLAST hits we got from the current query sequence during
          coarse BLAST.*/
        expanded_hits = expand_blast_hits(iterations, i, db);

        for (j = 0; j < expanded_hits->size; j++)
            ds_vector_append(oseqs, ds_vector_get(expanded_hits, j));
        ds_vector_free_no_data(expanded_hits);
    }

    write_fine_fasta(oseqs);

    for (i = 0; i < oseqs->size; i++)
        cb_hit_expansion_free(
          (struct cb_hit_expansion *)ds_vector_get(oseqs, i));
    ds_vector_free_no_data(oseqs);

    blast_fine(args, dbsize);

    if (!search_flags.hide_progress)
        fprintf(stderr, "\n"); //Make a newline after the progress bar

    for (i = 0; i < queries->size; i++)
        fasta_free_seq((struct fasta_seq *)ds_vector_get(queries, i));
    ds_vector_free_no_data(queries);

    if (search_flags.show_hit_info)
        fclose(test_hits_file);

    //Free the XML data and expanded hits
    for (i = 0; i < iterations->size; i++) {
        struct DSVector *iteration =
          (struct DSVector *)ds_vector_get(iterations, i);
        for (j = 0; j < iteration->size; j++) {
            struct hit *h = (struct hit *)ds_vector_get(iteration, j);
            ds_vector_free(h->hsps);
            free(h);
        }
    }
    ds_vector_free_no_data(iterations);

    cb_database_r_free(db);
    xmlFreeDoc(doc);

    /*Free the coarse BLAST results file if the --no-cleanup flag is not being
      used.*/
    if (!search_flags.no_cleanup)
        system("rm CaBLAST_temp_blast_results.xml");

    opt_args_free(args);
    opt_config_free(conf);

    return 0;
}
