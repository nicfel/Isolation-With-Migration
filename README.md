# Inference of species histories in the presence of gene flow

When populations become isolated, individuals of the different populations can diverge genetically over time. If this isolation is incomplete, genes can still be exchanged between these otherwise isolated populations. Current methods to infer past speciation processes over time either assume that the species tree and the order of speciation events is known *a priori* or ignore gene flow completely. Ignoring gene flow can lead to biased inference of species history if gene flow is present. 
We here present a new approach based on a recently introduced structured coalescent method that allows to co-infer the species tree alongside evolutionary parameters of interest while accounting for possible exchange of alleles after speciation. Using simulations, we show that our newly introduced approach is able to reliably infer species trees, speciation times, migration rates and effective population sizes unbiased from genetic sequence data. This is in contrast to the multi-species coalescent which shows biases in inference of species tree topologies as well as speciation times. We then apply our newly introduced method to infer the species history of six great apes species gene flow after population isolation. We find evidence of gene flow between bonobos and common chimpanzees.


## Source code

The source code of the approximate isolation-with-migration model in [Beast2](https://github.com/CompEvol/beast2) is part of the package [StarBeast2](https://github.com/genomescale/starbeast2)

## Jar

A compiled jar file **AIM.jar** is given in the Software folder and can be run as a normal BEAST2 file using the same command line input
