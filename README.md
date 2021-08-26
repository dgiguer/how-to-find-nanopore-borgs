# How to find borgs in nanopore data 

Daniel Giguere 

[![DOI](https://zenodo.org/badge/391498265.svg)](https://zenodo.org/badge/latestdoi/391498265)

**tl;dr at [bottom](https://github.com/dgiguer/how-to-find-nanopore-borgs#summary-tldr)**

### Background 

borgs are double stranded linear DNA elements recently proposed by [@themicrobeguy and others](https://twitter.com/themicrobeguy/status/1414202238537449473). they were identified using short-read metagenomic data - **so how can we find them using long-read metagenomic nanopore data?** tl;rdr: do an **all-vs-all alignment** for all long reads, **remove alignments where less than 90% of the read actually aligns**, then **create a network graph** of target and query read names. each cluster of reads represents the end of a linear dsDNA element. the number of reads in the cluster approximates the read depth for that genomic element. if this "end" of a linear element has reasonably high coverage and was already assembled into a long linear contig during your metagenome assembly, **congratulations! you found a linear dsDNA element, that could be a borg**. now you just need to make sure that it's not a eukaryotic chromosome, mega-phage, or other dsDNA element. 


### Filter alignments by query coverage

if you filter alignments by % query coverage (i.e., how much of that read actually aligns against the target), then circular dsDNA elements will have coverage that tails off at the beginning and end of the fasta file because a read that truly aligns to both sides will be broken into two alignments, resulting in two reads with around 50% query coverage (more or less). this is not the case for linear dsDNA elements because the DNA fragments actually all start at roughly the same spot. if you filter your read alignments by retaining only those with > 90% query coverage, the coverage depth plot will look like D) for circular elements in blue and E) for linear elements in orange.

![figure1](figs/coverage_figure.png)

### Create a network graph 


notice how in B) and C) that long reads all end at the same spot? that means they will align to all other reads that start at the "end" of a dsDNA linear with near 100% query coverage. these are retained by our previous filtering step. this should only occur throughout the rest of the genome by random chance. so if you have ~ 100X sequencing coverage of a linear dsDNA element, you should have a group of about 100 reads that all align to each other with near 100% query coverage. we know this because we have done this for other similar applications in a [eukaryotic genome assembly](https://www.biorxiv.org/content/10.1101/2021.05.04.442596v1.abstract), and have [completed circular bacterial genomes](https://www.biorxiv.org/content/10.1101/2020.04.08.032540v1) directly from metagenomes. 

in the eukaryotic genome assembly, we found **the number of clusters approximated the true number of chromosomes in the microalgae**, because it allowed us to estimate the number of linear dsDNA "ends" (i.e., telomeres in that case). we found about 100 clusters, which conveniently equals 25 chromosomes * 2 telomeres per chromosome * 2 haplotypes (which is what our assembly proposes). these clusters represent the "end" of any linear dsDNA element, including small linear dsDNA viral genomes, eukaryotic chromosomes, and most importantly in this context, borgs. this should allow you to identify borg candidate contigs in a sequence- and annotation-independent manner relatively quickly.

a visualization of the network graph is shown below in the context of overlapping read coverage, where the ends of a linear dsDNA element have a very strong signal of reads that all map 100% to each other (either query or target coverage), whereas this occurs infrequently in the middle of the chromosome. if you had a borg with 100X long-read coverage, you would therefore expect two clusters with about 100 reads in each. 

![figure2](figs/network_graph.png)

### Summary, tl;dr 

you can tell the diference between circular and linear dsDNA elements by visualizing coverage depth after filtering alignment by query coverage because of how minimap2 reports alignments. you can also identify all "ends" of linear dsDNA elements by creating a network graph of all-vs-all alignment output (filtered by high query coverage), where each cluster is actually a group of reads that aligns to all other reads in the group. combining both allows you to find all ends of linear dsDNA elements with high enough coverage - these are your borg candidates. you now have your borg candidates, and can verify that they are indeed borgs using the methods described in [@themicrobeguy's paper](https://twitter.com/themicrobeguy/status/1414202238537449473).

[questions? find borgs in nanopore data? click here to email me.](mailto:dgiguer@uwo.ca)

the citeable github repository (along with pseudo code) is available [here](https://github.com/dgiguer/how-to-find-nanopore-borgs).

### The work flow 

This is section is intended as guidance rather than actually providing a script. 

Use minimap2 to do all-vs-all read alignment. 

```
# output to paf for network graph of reads. use longer reads only (minimum 10 kb)
minimap2 -x ava-ont -t 40 reads_10kb.fastq.gz reads_10kb.fastq.gz > read_alignments.paf
```

You can filter reads by query coverage using a tool Alec and I developed called [gerenuq[(https://github.com/abahcheli/gerenuq). 

```
# retain only read alignments with 90% query coverage. this also reduces file size for R
gerenuq -i read_alignments.paf -m 0.9 -o read_alignments_filtered.paf
```

Now bring this in to R for the network graph using the iGraph package. 

```
R

d <- read.table("read_alignments_filtered.paf", sep = "\t", header= FALSE, row.names=NULL)

# 

library(igraph)

d.sub <- data.frame(from=d$V1, to=d$V6)


g <- graph_from_data_frame(d.sub)
clu <- components(g)

# show individual groups of reads 
groups(clu)

# plot a histogram of reads to see if high coverage linear dsDNA exist.  
final <- vector()

for (i in seq(length(groups(clu)))) {
    
    final[i] <- length(groups(clu)[[i]])
    
}

hist(final, breaks = 99)
```

If high coverage groups exist, these will represent ends of linear dsDNA elements. Find the contig in your assembly by mapping a read or two to it, then verify that it's a borg using the criteria described in [@themicrobeguy](https://twitter.com/themicrobeguy)'s pre-print!
 

