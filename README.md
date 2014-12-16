Codon Optimization
==================

Generally, codon optimization for expression of a gene in a non-native expression system is achieved by inserting the most-frequently used codons for the corresponding amino-acids. However, recent studies suggest the existance of encoded information in how rare and frequenty used codons are distributed in the gene [1,2]. This encoded information is lost during such "blind" optimization procedures. Here, I have developed a new methodology to optimize a given DNA sequence such that its codon usage pattern in the expression system/organism closely matches that of the source organism.

[1] Clarke T.F., Clark P.L. Rare codons cluster. PLoS ONE. 2008;3:e3412. 

[2] Chartier M, Gaudreault F, Najmanovich R. Large-scale analysis of conserved rare codon clusters suggests an involvement in co-translational molecular recognition events. Bioinformatics 2012;28(11):1438-1445

Usage
=====
Type
`python codon-optimize.py -h` 
to get help about usage

Package Requirements
====================
1. Python 2.7
2. Python external library dependencies (install using pip or easy_install)
 
 a) matplotlib
 ```Shell
 pip install matplotlib
 ```
 b) numpy (numeric python)
  ```Shell
  pip install numpy
  ```

Author
======
Venkatramanan Krishnamani (venky.krishna@me.com)
