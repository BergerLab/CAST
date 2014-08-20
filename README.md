This is an implementation of the Compressively Accelerated Search Tools (CAST)
in C99 built from the C implementation of CaBLASTP, which can be found on
https://github.com/BurntSushi/c_cablastp.  It is made to work like Po-Ru Loh's
version of the code, which can be found in a prototype at
 http://cast.csail.mit.edu/.

The compressed databases created by this code are not in any way compatible 
with the compressed databases produced by the C++ implementation. Namely, the
data itself should be the same (or similar), but its representation on disk is
completely different. There may be other, subtler differences as well.

**Multithreaded compression introduces some variance into the size of the
  resulting compressed database.**

The official Go implementation of CaBLASTP can be found here:
https://github.com/BurntSushi/cablastp

**Notes about CaBLAT:**
* Because of the smaller database size in CaBLAT, CaBLAT may produce hits
  that do not occur when running BLAT on a full uncompressed database.
* For CaBLAT, it is recommended to pass in the flags --number-queries and
  --number-targets to prevent errors in fine BLAT related to duplicate query and
  target sequence names.
* Because of the differences between how fine BLAT is implemented in this
  version of CaBLAT and in Po-Ru's prototype version, there are some
  differences in the results of this version and the prototype.

Installation
============
CAST depends on three libraries: `opt`, `ds` and `pthread`. `pthread` 
should be installed via your system's package manager. `opt` and `ds` can be 
found in Andrew Gallant's [clibs respository]
(https://github.com/BurntSushi/clibs), which is included in the clibs directory
in src but is modified to produce .a files instead of .so files.

Briefly, the following commands should get CAST into a working state:

```bash
git clone git://github.com/BergerLab/c_cablast
cd c_cablast/src/clibs
make
cd ../..
make
./cablast-compress --help
```

The usage of `cablast-compress` is:

```bash
cablast-compress [flags] database-directory fasta-file [fasta-file ...]
```

The usage of `cablast-search` is:
```bash
./cablast-search [flags] database-dir fasta-file [ --blast_args BLASTN_ARGUMENTS ]
```

The usage of `cablat` is:
```bash
./cablat [flags] database-dir fasta-file results-file [ --blat_args BLAT_ARGUMENTS ]
```
