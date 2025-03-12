# Lab 9 Analyzing the SARS-CoV-2 genome

Let's first set up our lab environment. Log into the cluster.

```bash
mkdir lab_9
cd lab_9
```
We also need to set up a few packages, as we will be using python! Yay!

```bash
module load python
pip install biopython matplotlib jupyterlab
```

A lot of things will print to screen if it is working. It may take a few minutes but we will get there.

# Part 1: Taking a smart approach and validating things in a well-known genome.

Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) is a strain of coronavirus that causes COVID-19. Throughout these labs, we will analyze this virus and try to understand its origins and inner workings. We will implement the necessary bioinformatics tools and apply them to further our understanding of this pesky little virus.

First, we will analyze the main component of every organism - its genetic material. Our focus will be on the genes, parts of the genetic material that code for proteins. Proteins are the main macromolecular actors in every organism.

Why are we so interested in genes?

Genes dictate the behavior of an organism, such as replication, viral assembly, and even innate immune evasion. If we compare the genes from this new virus with genes from other known viruses, we can get a good idea of how this virus works and maybe even how to stop it. When these genes are translated into proteins, they start acting out their function. Some proteins can attach to human cells and allow viruses to enter them. If we can figure out which genes these are, they will make good candidates for drug targets.

We can find potential genes in a genome by looking for common patterns shared across all genes. However, validating that these potential candidates are, in fact, real genes requires experimental confirmation.

In this lab, we combine tools that we already know and begin to explore new tools. Our first goal is to find potential genes in the SARS-CoV-2 genome and take a first stab at figuring out what their corresponding proteins do. When looking for genes in an unknown genome, we consider all the Open Reading Frames (ORFs) as potential gene candidates. However, as we remember from our braker and BUSCO labs, many of our ORFs will not be true genes, and since each ORF must be experimentally validated in the lab, we will try to reduce the number of ORFs only to include the most likely gene candidates. One approach for filtering ORFs is to perform a permutation test to determine a threshold indicating the minimum length of ORFs. This will reduce the number of false positives that we generate. So, our approach might be to first find the ORFs, remove ORFs that are likely too short, and then determine what each ORF might do.

However, before we jump right into the SARS-CoV-2 genome, we first want to convince ourselves that this is, in fact, a reasonable approach that produces good results. Since we are going to pretend that SARS-CoV-2 is an unknown virus, we can't check our results to see if they are correct. So instead, we will validate our approach on one of the most well-understood organisms in existence: E. Coli. Here, we will be able to check how many of our ORFs are true genes and how well our threshold selection method actually works.

Once we've convinced ourselves that this is, in fact, a good approach, we will start exploring the SARS-CoV-2 genome. Finally, once we have found and filtered our coronavirus ORFs, we will implement a simple classification technique to determine the function of their corresponding proteins. We will explore computational techniques for determining protein function in the next homework.


First we need to initiate python. 
```bash
python
```


```python
from Bio import Entrez

Entrez.email = "your_email@example.com"  # Always provide your email when using Entrez

# Fetch the E. coli K-12 MG1655 genome (NC_000913.3)
handle = Entrez.efetch(db="nucleotide", id="NC_000913.3", rettype="gbwithparts", retmode="text")

# Read and save the genome data
with open("ecoli_genome.gb", "w") as output_file:
    output_file.write(handle.read())

handle.close()
```



[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
