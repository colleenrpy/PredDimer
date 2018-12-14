# PredDimer
PredDimer is an algorithm to predict heterodimeric protein complexes with pairwise kernel function and Support Vector Machines (SVM).
Please read paper [Improving prediction of heterodimeric protein complexes using combination with pairwise kernel](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2017-5) for details.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

What things you need to install the software and how to install them

#### Environment

* Library: libsvm-3.0.2
* Programming language: C++

#### Data

* cyc2008.data - downloaded from [cyc2008](http://wodaklab.org/cyc2008/)
* pro200600448_3_s.csv - downloaded from [WI-PHI](http://www.wiley-vch.de/contents/jc_2120/2007/pro200600448_s.html)
* uniprot_sprot_fungi.dat - downloaded from [UniProtKB](https://www.uniprot.org/uniprot/)

### Usage

```
g++ pairwise-kernel.cpp -o dimer
./dimer 1 2 3 4 5 8 9 a filename1.txt filename2.txt
./svm-train -q -t 4 -v 10 -w1 b -w-1 c filename2.txt
```
a: 
b: weight of positive examples
c: weight of negative examples
