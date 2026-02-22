# evop-phylogenomics

This hands-on course will introduce you to the awesomeness of phylogenomics.

You will learn the most fundamental steps in the phylogenomics pipeline, allowing you to go from a bunch of sequences from different species to a phylogenetic tree that represents how these species are related to each other in evolutionary terms. 

The phylogenomics pipeline can become very complex, many additional steps might be included (particularly at the stage of dataset assembly!) and some analyses can take weeks to complete. Pipelines are modular, meaning they can (and should) be improved and modified to fit the particular data or question at hand. Here, we will cover the basics of a phylogenomic pipeline: identification of homologs, multiple sequence alignment, alignment trimming, and phylogenetic inference with concatenation and coalescent methods.



## Objective and data


We will use a dataset extracted [from this paper](https://royalsocietypublishing.org/rspb/article/288/1963/20212168/79242/Unexpected-cryptic-species-among-streptophyte). The starting point is a subset of proteins obtained from genomes/transcriptomes for algae and plants. We aim to reconstruct these species' phylogeny using data concatenation and coalescent approaches. In practice, we are using a small subset of the original data to to speed up computations.

Let's start by cloning this repository. The starting data can be found in the `proteomes` folder.

```
# clone repository, data, and software
git clone https://github.com/iirisarri/evop-phylogenomics.git
```



## Inferring ortholog groups

To identify homologs among all proteins, we will use [OrthoFinder](https://github.com/davidemms/OrthoFinder). Try to figure out the command by yourself by looking at the [online documentation](https://davidemms.github.io/menu/tutorials.html) or by the help provided when typing `orthofinder -h`. 


<details>
  <summary>Need help?</summary>
  To run OrtoFinder with default parameters, just provide the path to the folder containing the input proteomes.
  
  
```
orthofinder -f proteomes
```
</details>

Look for Orthofinder's results inside your `proteomes` folder. A bunch of interesting information is contained there. The Orthogroups can be found in `Orthogroup_Sequences`. Each file corresponds to one orthogroup ("gene") containing one sequence per species.

Let's check what the orthogroups look like. How many sequences does each orthogroup contain? Do you see anything unexpected? In case, remove any orthogroup if it contains more sequences than the total number of taxa.

<details>
  <summary>Need help?</summary>
  
```
less -S OG0000006.fa

grep -c ">" *.fa
```
</details>

You will see that sequences are named `Genus@sequence_name`. To concatenate genes at a later stage, sequences of the same taxa need to have uniform names. Thus, we can homogenize them already now. Can you use bash to remove the bits that differ between different sequences of the same species?

<details>
  <summary>Need help?</summary>
  
```
for f in *fa; do sed -E '/>/ s/@.+//g' $f > out; mv out $f; done
# Explanation: in headers (lines starting with ">"), remove "@" and everything afterwards
```
</details>

Don't forget to check the output: is your command doing what you want?


**NOTE ABOUT ORTHOLOGY**: Ensuring orthology is a difficult issue and often using a tool like Orthofinder might not be enough. Paralogy is a tricky business! Research has shown (e.g. [here](https://www.nature.com/articles/s41559-017-0126) [here](https://academic.oup.com/sysbio/article/71/1/105/6275704?login=false) or [here](https://academic.oup.com/mbe/article/36/6/1344/5418531)) that including paralogs can bias phylogenetic relationships and molecular clock estimates, particularly when the phylogenetic signal is weak. Paralogs should always be removed before phylogenetic inference. But identifying them can be difficult and time-consuming. One could build single-gene trees and look for sequences producing extremely long branches or those that cluster outside of the remaining sequences. [Automatic pipelines](https://github.com/fethalen/phylopypruner) also exist.



## Pre-alignment and quality filtering (optional)


Often, transcriptomes and genomes have stretches of erroneous, non-homologous amino acids or nucleotides, produced by sequencing errors, assembly errors, or errors in genome annotation. Until recently, these types of errors had been mostly ignored because no automatic tool could deal with them.

We will use [PREQUAL](https://academic.oup.com/bioinformatics/article/34/22/3929/5026659?login=true), a software that takes sets of (homologous) unaligned sequences and identifies sequence stretches (amino acids or codons) sharing no evidence of (residue) homology, which are then masked in the output. Note that homology can be invoked at the level of sequences as well as of residues (amino acids or nucleotides). 

Running PREQUAL for each set orthogroup is easy if we use a for loop in Bash or the command `parallel` (you might need to find the path to the program `prequal`):
```sh
parallel 'prequal {}' ::: *fa

for f in *fa; do prequal $f ; done
```
The filtered (masked) alignments are in .filtered, whereas .prequal contains relevant information such as the number of residues filtered.



## Multiple sequence alignment


The next step is to infer multiple sequence alignments from orthogroups. Multiple sequence alignments allow us to *propose* which amino acids/ nucleotides are homologous. A simple yet accurate tool is [MAFFT](https://mafft.cbrc.jp/alignment/server/).

We will align gene files separately using a for loop:

```sh
for f in *filtered; do mafft $f > $f.mafft; done
```



## Alignment trimming


Some gene regions (e.g., fast-evolving) are difficult to align and thus positional homology can be uncertain. It is unclear (i.e., problem-specific) whether trimming suspicious regions [improves](https://academic.oup.com/sysbio/article/56/4/564/1682121) or [worsens](https://academic.oup.com/sysbio/article/64/5/778/1685763) tree inference. However, gently trimming very incomplete positions (e.g. with >95% gaps) will speed up computation in the next steps without a significant loss of phylogenetic information.

To trim alignment positions we can use [ClipKIT]([https://github.com/JLSteenwyk/ClipKIT]) but several other software are also available.


```sh
for f in *mafft; do clipkit -g 0.95 $f -o $f.g95; done
```

While diving into phylogenomic pipelines, it is always advisable to check a few intermediate results to ensure we are doing what we should. Multiple sequence alignments can be visualized in [SeaView](http://doua.prabi.fr/software/seaview) or [AliView](https://github.com/AliView/AliView) on your local machine. Also, one could have a quick look at alignments using command line tools (`less -S`).



## Concatenate alignment


To infer our phylogenomic tree we need to concatenate the trimmed single-gene alignments we generated. There are many tools that you can use for this step (e.g [concat_fasta.pl](https://github.com/santiagosnchez/concat_fasta) or [catsequences](https://github.com/ChrisCreevey/catsequences)). Here, we will use [PhyKIT](https://jlsteenwyk.com/PhyKIT/usage/index.html), which will read (option `-a`) a text file as input that needs to contain the file names of alignments to be concatenated and `-p` will indicate the prefix of the output file. Three output files are generated: a concatenated alignment, a partition file, and an occupancy file. Try to build the command by yourself!


<details>
  <summary>Need help?</summary>
  
  
```
phykit create_concat -a concat_files -p concatenation
```
</details>



Is your concatenated file what you expected? It should contain 21 taxa and 108 genes. If that is not the case, you might have forgotten someting on the way. Looking good? Then your concatenated dataset is ready to rock!!



## Concatenation: Maximum likelihood


One of the most common approaches in phylogenomics is gene concatenation: the signal from multiple genes is "pooled" together with the aim of increasing resolution power. This method is best when among-gene discordance is low.

We will use [IQTREE](http://www.iqtree.org/), an efficient and accurate software for maximum likelihood analysis. Another great alternative is [RAxML](https://github.com/stamatak/standard-RAxML). The most simple analysis is to treat the concatenated dataset as a single homogeneous entity. We need to provide the number of threads to use (`-nt 1`) input alignment (`-s`), tell IQTREE to select the best-fit evolutionary model with BIC (`-m TEST -merit BIC -msub nuclear`) and ask for branch support measures such as non-parametric bootstrapping and approximate likelihood ratio test (`-bb 1000 -alrt 1000 -bnni`):

```sh
iqtree3 -s concatenation.fa -m TEST -msub nuclear -bb 1000 -alrt 1000 -nt AUTO -bnni -pre unpartitioned
```

A more sophisticated approach would be to perform a partitioned maximum likelihood analysis, where different genes (or other data partitions) are allowed to have different evolutionary models. This should provide a better fit to the data but will increase the number of parameters too. To launch this analysis we need to provide a file containing the coordinates of the partitions (`-p`) and we can ask IQTREE to select the best-fit models for each partition, in this case, according to AICc (more suitable for shorter alignments).

```sh
iqtree3 -s concatenation.fa -p FcC_supermatrix_partition.txt -m TEST -msub nuclear -merit AICc -bb 1000 -alrt 1000 -nt AUTO -bnni -pre partitioned
```

Congratulations!! If everything went well, you should get your maximum likelihood estimation of the plant phylogeny (`.treefile`)! Looking into the file you will see a tree in parenthetical (newick) format. See below how to create a graphical representation of your tree.



## Coalescence analysis


An alternative to concatenation is to use a multispecies coalescent approach. Unlike maximum likelihood, coalescent methods account for incomplete lineage sorting (ILS; an expected outcome of evolving populations). These methods are particularly useful when we expect high levels of ILS, e.g. when speciation events are rapid and leave little time for allele coalescence.

We will use [ASTRAL](https://github.com/smirarab/ASTRAL), a widely used tool that scales up well to phylogenomic datasets. It takes a set of gene trees as input and will generate the coalescent "species tree". ASTRAL assumes that gene trees are estimated without error.

Thus, before running ASTRAL, we will need to estimate individual gene trees. This can be easily done by calling IQTREE in a for loop:

```sh
parallel 'iqtree3  -s {} -m TEST -msub nuclear -merit AICc' ::: *g95

for f in *g95; do iqtree3  -s $f -m TEST -msub nuclear -merit AICc -nt AUTO; done
```

After all gene trees are inferred, we should put them all into a single file:

```sh
cat *g95.treefile > my_gene_trees.tre
```

Now running ASTRAL is trivial, providing the input file with the gene trees and the desired output file name (you might need to find the path to the program `astral.5.7.8.jar`):

```sh
java -jar astral.5.7.8.jar -i my_gene_trees.tre -o species_tree_ASTRAL.tre 2> out.log
```

Congratulations!! You just got your coalescent species tree!! Is it different from the concatenated maximum likelihood trees? 



## Tree visualization


Trees are just text files representing relationships with parentheses; did you see that already? But it is more practical to plot them as a graph, for which we can use tools such as [iTOL](https://itol.embl.de), [FigTree](https://github.com/rambaut/figtree/releases), [TreeViewer](https://treeviewer.org/), [iroki](https://www.iroki.net/), or R (e.g. [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html), [phytools](http://www.phytools.org/); see provided R code).

Upload your trees to iTOL. Trees need to be rooted with an outgroup, in this case, at the branch that separates Chlorophyta from Streptophyta. In iToL "Tree Structure/Reroot the tree here". Branch support values can be shown under the "Advanced" menu. The tree can be modified in many other ways, and finally, a graphical tree can be exported. In Figtree, a rooted tree can be saved using File/Export Trees.../"Save as currently displayed".

[Well done!](https://media0.giphy.com/media/v1.Y2lkPTc5MGI3NjExdG9mZjBqNXU5Y2xxdnpyMWQ0d3RrM2I2aDZwNmptOGNscHd3NnFtNSZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9dg/qjTLsLjcu5PpscKwge/giphy.gif)



## Ancestral character state reconstruction

What can we do with the trees we obtained? So many things!!! The very basics: to investigate how species are related to each other in the green lineage (Chloroplastida), and interpret their evolutionary history.

Something more advanced is to study the evolution of a given trait in the group of interest. Here, we will reconstruct the evolution of multicellularity in the green lineage using [phytools](https://peerj.com/articles/16505/]). In particular, we can use the function `ace` to reconstruct ancestral character states of multicellularity using (a) trait information from extant species at the tips of the phylogeny and (b) the evolutionary history of the group, i.e. the phylogeny, that reflects the evolutionary relationships among the species (topology) and their relative divergences (branch lengths).

Thus, the input files are (a) a text file indicating whether a given species is unicellular (0) or multicellular (1) (`multicellularity.txt`) and (b) one of the phylogenies that you inferred above. Very important: the phylogeny should be properly rooted, in this case at the branch that separates Chlorophyta from Streptophyta. The tree might be rerooted in FigTree and exported in Newick (File/Export Trees.../"Save as currently displayed").

A phytools tutorial for ancestral character state reconstruction under different models is available [here](http://www.phytools.org/Cordoba2017/ex/8/Anc-states-discrete.html). Try to follow the tutorial until you can identify the first multicellular ancestor. If you find difficulties, there might be help somewhere in this repository :-)

Now, a more advanced question: will different trees reconstruct the same evolutionary history for the evolution of multicellularity? Look at Figure 3 of [this paper](https://www.cell.com/current-biology/fulltext/S0960-9822(23)01770-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0960982223017700%3Fshowall%3Dtrue). What do you notice?




## Software links

* Orthofinder (https://github.com/davidemms/OrthoFinder)
* PREQUAL (https://github.com/simonwhelan/prequal)
* MAFFT (https://mafft.cbrc.jp/alignment/software/source.html)
* MUSCLE v5 (https://github.com/rcedgar/muscle)
* ClipKIT (https://github.com/JLSteenwyk/ClipKIT)
* PhyKIT (https://jlsteenwyk.com/PhyKIT/usage/index.html)
* IQTREE (http://www.iqtree.org/)
* ASTRAL (https://github.com/smirarab/ASTRAL)
* FIGTree V1.4.4 (https://github.com/rambaut/figtree/releases) 
* TreeViewer (https://treeviewer.org/)
* iTOL (https://itol.embl.de/)
* SeaView (https://doua.prabi.fr/software/seaview_data/seaview5-64.tgz)
  
