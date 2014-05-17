CC=gcc
CFLAGS=-g -O3 -ansi -Wall -Wextra -pedantic -lpthread -I. -I/usr/include/libxml2
LDLIBS=-lds -lpthread -lopt -lxml2

COMPRESS_OBJS=align.o DNAalphabet.o DNAmatrix.o \
							bitpack.o coarse.o compressed.o compression.o database.o DNAutils.o edit_scripts.o fasta.o flags.o \
							progress-bar.o seeds.o seq.o util.o
COMPRESS_HEADERS=align.o coarse.h compressed.h compression.h \
							database.h DNAalphabet.h fasta.h flags.h link_to_coarse.h link_to_compressed.h progress-bar.h seeds.h seq.h util.h

all: cablast-compress cablast-decompress cablast-search

compress: DNAalphabet.h cablast-compress

cablast-compress: cablast-compress.o $(COMPRESS_HEADERS) $(COMPRESS_OBJS)

cablast-compress.o: cablast-compress.c database.h DNAalphabet.h fasta.h



DECOMPRESS_OBJS=align.o \
							bitpack.o coarse.o compressed.o compression.o database.o decompression.o DNAalphabet.o DNAmatrix.o DNAutils.o edit_scripts.o fasta.o flags.o \
							progress-bar.o seeds.o seq.o util.o
DECOMPRESS_HEADERS=align.o coarse.h compressed.h compression.h \
							database.h decompression.h DNAalphabet.h fasta.h flags.h link_to_coarse.h link_to_compressed.h progress-bar.h seeds.h seq.h util.h

decompress: DNAalphabet.h cablast-decompress

cablast-decompress: cablast-decompress.o $(DECOMPRESS_HEADERS) $(DECOMPRESS_OBJS)

cablast-decompress.o: cablast-decompress.c database.h DNAalphabet.h fasta.h



SEARCH_OBJS=align.o bitpack.o coarse.o compressed.o compression.o database.o decompression.o DNAalphabet.o DNAmatrix.o DNAutils.o edit_scripts.o fasta.o flags.o \
							progress-bar.o seeds.o seq.o util.o
SEARCH_HEADERS=align.o coarse.h compressed.h compression.h \
							database.h decompression.h DNAalphabet.h fasta.h flags.h link_to_coarse.h link_to_compressed.h progress-bar.h seeds.h seq.h util.h xml.h

search: DNAalphabet.h cablast-search

cablast-search: cablast-search.o $(SEARCH_HEADERS) $(SEARCH_OBJS)

cablast-search.o: cablast-search.c database.h DNAalphabet.h fasta.h


align.o: align.c align.h DNAalphabet.h
bitpack.o: bitpack.c bitpack.h
coarse.o: coarse.c coarse.h link_to_compressed.h seq.h
compressed.o: compressed.c compressed.h link_to_coarse.h
compression.o: compression.c compression.h align.h coarse.h compressed.h seq.h
database.o: database.c database.h coarse.h compressed.h link_to_coarse.h link_to_compressed.h
decompression.o: decompression.c decompression.h
DNAalphabet.o: DNAalphabet.c DNAalphabet.h
DNAmatrix.o: DNAmatrix.c DNAalphabet.h
DNAutils.o: DNAutils.c DNAutils.h
edit_scripts.o: link_to_coarse.h edit_scripts.c edit_scripts.h
fasta.o: fasta.c fasta.h util.h
flags.o: flags.c flags.h util.h
progress-bar.o: progress-bar.c progress-bar.h
seeds.o: seeds.c seeds.h
seq.o: seq.c seq.h
uitl.o: util.c util.h

#blosum62_matrix.c: ../scripts/mkBlosum
#	../scripts/mkBlosum > blosum62_matrix.c

tests: test-extension test-nw test-ungapped

# test-extension: tests/test-extension.o align.o blosum62.o blosum62_matrix.o \ 
								# compression.o 
	# $(CC) $(LDFLAGS) \ 
		# tests/test-extension.o align.o blosum62.o blosum62_matrix.o compression.o \ 
		# $(LDLIBS) \ 
		# -o test-extension 

test-extension: tests/test-extension.o $(COMPRESS_OBJS)
	$(CC) $(LDFLAGS) \
		tests/test-extension.o $(COMPRESS_OBJS) \
		$(LDLIBS) \
		-o test-extension

test-nw: tests/test-nw.o align.o DNAalphabet.o DNAmatrix.o
	$(CC) $(LDFLAGS) \
		tests/test-nw.o align.o DNAalphabet.o DNAmatrix.o \
		$(LDLIBS) \
		-o test-nw

test-ungapped: tests/test-ungapped.o align.o DNAalphabet.o DNAmatrix.o
	$(CC) $(LDFLAGS) \
		tests/test-ungapped.o align.o DNAalphabet.o DNAmatrix.o \
		$(LDLIBS) \
		-o test-ungapped

clean:
	rm -f *.o tests/*.o
	rm -f cablast-compress
	rm -f test-*

loc:
	find ./ -name '*.[ch]' -print0 | xargs -0 wc -l

