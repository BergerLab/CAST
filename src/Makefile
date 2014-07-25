CC=gcc
CFLAGS=-g -O3 -std=c99 -Wall -Wextra -pedantic -lpthread -I. \
       -I/usr/include/libxml2
LDLIBS=clibs/libds/libds.a -lpthread clibs/libopt/libopt.a -lxml2

COMPRESS_OBJS=align.o bitpack.o coarse.o compressed.o \
              compression.o database.o DNAutils.o edit_scripts.o fasta.o \
              flags.o link_to_coarse.o link_to_compressed.o progress-bar.o \
              seeds.o seq.o util.o

COMPRESS_HEADERS=align.o coarse.h compressed.h compression.h database.h \
                 fasta.h flags.h link_to_coarse.h \
                 link_to_compressed.h progress-bar.h seeds.h seq.h util.h

all: cablast-compress cablast-decompress cablast-search cablat-search

compress: cablast-compress

cablast-compress: cablast-compress.o $(COMPRESS_HEADERS) $(COMPRESS_OBJS)

cablast-compress.o: cablast-compress.c database.h fasta.h



DECOMPRESS_OBJS=align.o bitpack.o coarse.o compressed.o compression.o \
                database.o decompression.o DNAutils.o \
                edit_scripts.o fasta.o flags.o link_to_coarse.o \
                link_to_compressed.o progress-bar.o seeds.o seq.o util.o

DECOMPRESS_HEADERS=align.o coarse.h compressed.h compression.h database.h \
                   decompression.h fasta.h flags.h \
                   link_to_coarse.h link_to_compressed.h progress-bar.h \
                   seeds.h seq.h util.h

decompress: cablast-decompress

cablast-decompress: cablast-decompress.o $(DECOMPRESS_HEADERS) \
                                         $(DECOMPRESS_OBJS)

cablast-decompress.o: cablast-decompress.c database.h fasta.h



SEARCH_OBJS=align.o bitpack.o coarse.o compressed.o compression.o database.o \
            decompression.o DNAutils.o edit_scripts.o fasta.o \
            flags.o link_to_coarse.o link_to_compressed.o progress-bar.o \
            seeds.o seq.o util.o

SEARCH_HEADERS=align.o coarse.h compressed.h compression.h database.h \
               decompression.h fasta.h flags.h link_to_coarse.h \
               link_to_compressed.h progress-bar.h seeds.h seq.h util.h xml.h

search: cablast-search

cablast-search: cablast-search.o $(SEARCH_HEADERS) $(SEARCH_OBJS)

cablast-search.o: cablast-search.c database.h fasta.h



CABLAT_OBJS=align.o bitpack.o coarse.o compressed.o compression.o database.o \
            decompression.o DNAutils.o edit_scripts.o fasta.o \
            flags.o link_to_coarse.o link_to_compressed.o progress-bar.o \
            seeds.o seq.o util.o

CABLAT_HEADERS=align.o coarse.h compressed.h compression.h database.h \
               decompression.h fasta.h flags.h link_to_coarse.h \
               link_to_compressed.h progress-bar.h seeds.h seq.h util.h xml.h

cablat-search: cablat

cablat: cablat.o $(CABLAT_HEADERS) $(CABLAT_OBJS)

cablat.o: cablat.c database.h fasta.h


align.o: DNAutils.c DNAutils.h align.c align.h
bitpack.o: bitpack.c bitpack.h
coarse.o: coarse.c coarse.h link_to_compressed.h seq.h
compressed.o: compressed.c compressed.h link_to_coarse.h
compression.o: compression.c compression.h align.h coarse.h compressed.h seq.h
database.o: database.c database.h coarse.h compressed.h \
            link_to_coarse.h link_to_compressed.h
decompression.o: decompression.c decompression.h
DNAutils.o: DNAutils.c DNAutils.h
edit_scripts.o: link_to_coarse.h edit_scripts.c edit_scripts.h
fasta.o: fasta.c fasta.h util.h
flags.o: flags.c flags.h util.h
link_to_coarse.o: link_to_coarse.c link_to_coarse.h
link_to_compressed.o: link_to_compressed.c link_to_compressed.h
progress-bar.o: progress-bar.c progress-bar.h
seeds.o: seeds.c seeds.h
seq.o: seq.c seq.h
uitl.o: util.c util.h

loc: find ./ -name '*.[ch]' -print0 | xargs -0 wc -l
