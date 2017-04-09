# slaMEM

*slaMEM* is a tool used to efficiently retrieve _MEMs_ (`Maximal Exact Matches`)
between a reference genome sequence and one or more query sequences, similarly to
these software tools:
* [E-MEM](http://www.csd.uwo.ca/~ilie/E-MEM/)
* [MUMmer](http://mummer.sourceforge.net/)
* [essaMEM](https://github.com/readmapping/essaMEM)
* [sparseMEM](https://github.com/zia1138/sparseMEM)
* [backwardMEM](http://www.uni-ulm.de/in/theo/research/seqana.html#c102393)

*slaMEM* relies on an [FM-Index][1] together with a new data structure called
*SSI[LCP][2]* (`Sampled Search Intervals from Longest Common Prefixes`) to store
information about _parent intervals_ in a time- and space-efficient way.

*slaMEM* also includes an useful feature to display the locations of the found
MEMs, generating images like the one below.

[![MEMs of 57 E.coli strains](https://github.com/fjdf/slaMEM/blob/master/images/all-ecoli-strains-mems_small.jpg)](https://github.com/fjdf/slaMEM/blob/master/images/all-ecoli-strains-mems.bmp?raw=true)

### Reference
If you use *slaMEM*, please cite:
[Fernandes, Francisco and Ana T. Freitas. **slaMEM: efficient retrieval of maximal exact matches using a sampled LCP array**. Bioinformatics 30.4 (2014): 464-471.](http://dx.doi.org/10.1093/bioinformatics/btt706)

[1]: https://en.wikipedia.org/wiki/FM-index
[2]: https://en.wikipedia.org/wiki/LCP_array

### Manual

#### Install
```bash
make
```
#### Usage
```bash
./slaMEM (<options>) <reference_file> <query_file(s)>
```
##### Options:
- `mem` : find MEMs: any number of occurrences in both ref and query (default)
- `mam` : find MAMs: unique in ref but any number in query
- `l`   : minimum match length (default=20)
- `o`   : output file name (default="*-mems.txt")
- `b`   : process both forward and reverse strands
- `n`   : discard 'N' characters in the sequences
- `m`   : minimum sequence size (e.g. to ignore small scaffolds)
- `r`   : load only the reference(s) whose name(s) contain(s) this string
##### Extra:
- `v` : generate MEMs map image from this MEMs file
##### Example:
```bash
./slaMEM -b -l 10 ./ref.fna ./query.fna
./slaMEM -v ./ref-mems.txt ./ref.fna ./query.fna
```
