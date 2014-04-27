# slaMEM

*slaMEM* is a tool used to efficiently retrieve _MEMs_ (`Maximal Exact Matches`)
between a reference sequence and one or more query sequences, similarly to
these software tools:
* [MUMmer](http://mummer.sourceforge.net/)
* [essaMEM] (https://github.ugent.be/ComputationalBiology/essaMEM)
* [sparseMEM] (http://compbio.cs.princeton.edu/mems/)
* [backwardMEM](http://www.uni-ulm.de/in/theo/research/seqana.html#c102393)

*slaMEM* relies on an [FM-Index][1] together with a new data structure called
*SSILCP* (`Sampled Search Intervals from Longest Common Prefixes`) to store
information about _parent intervals_ in a space and time efficient way.

*slaMEM* also includes an useful feature to display the locations of the found
MEMs, generating images like the one below.

![MEMs of 57 E.coli strains](https://raw.github.com/fjdf/slaMEM/master/images/all-ecoli-strains-mems_small.jpg "E.coli MEMs")

### Install
Simply compile it with `make`.

### Usage
The basic command line is:
```
./slaMEM <reference_file> <query_files>
```
For a full list of options run `./slaMEM` with no arguments.

### Reference
To reference *slaMEM* please cite:

[Fernandes,F. and Freitas,A.T. (2013) slaMEM: Efficient retrieval of Maximal Exact Matches using a Sampled LCP Array. Bioinformatics.](http://dx.doi.org/10.1093/bioinformatics/btt706)


[1]: http://dl.acm.org/citation.cfm?id=796543


---

> Copyright (C) 2013 INESC-ID KDBio (fjdf@kdbio.inesc-id.pt)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
