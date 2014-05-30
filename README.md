This is an implementation of CaBLAST in ANSI C built from the C implementation
of CaBLASTP found on https://github.com/BurntSushi/c_cablastp.  It is made to
work like Po-Ru Loh's version of the code, which can be found in a prototype at
http://cast.csail.mit.edu/.

The compressed databases created by this code are not in any way compatible 
with the compressed databases produced by the official implementation. Namely, 
the data itself should be the same (or similar), but its representation on disk 
is completely different. There may be other, subtler differences as well.

**This code is currently under development.**

There are no benchmarks.

The official Go implementation of CaBLASTP can be found here:
https://github.com/BurntSushi/cablastp


Installation
============
c_cablast depends on three libraries: `opt`, `ds` and `pthread`. `pthread` 
should be installed via your system's package manager. `opt` and `ds` can be 
found in Andrew Gallant's [clibs respository](https://github.com/BurntSushi/clibs).

Briefly, the following commands should get c_cablast into a working state:

```bash
mkdir c_cablast
git clone git://github.com/BurntSushi/clibs
git clone git://github.com/BergerLab/c_cablast
cd clibs
make
export C_INCLUDE_PATH=$(pwd)/include
export LIBRARY_PATH=$(pwd)/lib
cd ../c_cablast
make
./cablast-compress --help
```

The usage of `cablast-compress` is:

```bash
cablast-compress [flags] database-directory fasta-file [fasta-file ...]
```


Current progress:  Can compress and decompress two same-direction matches, two reverse-complement matches,
one same-direction match and one reverse-complement match, and two matches that are at least one chunk
apart.  On real test data, can currently compress and decompress all ten Brucella sequences from the test
file, producing as output coarse FASTA files, links tables, and seeds tables that perfectly match these data
in Po-Ru Loh's C++ version of the code, and the order of sequences does not affect whether or not compression
and decompression work correctly.
