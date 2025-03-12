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


## Escherichia coli
Escherichia coli (E. Coli) is a bacteria commonly found in the human intestine. Most strains are harmless to humans, at worse, causing food poisoning and diarrhea. They can survive outside a host for only a short period of time, making it a potential indicator of fecal contamination.

Over the years, E. Coli has been intensely studied and is probably one of the most well-understood organisms in existence. We've learned how to grow them in an optimal environment where they can reproduce up to once every 20 minutes. Due to their rapid growth and easy manipulation, biologists often use them to produce recombinant proteins.
![image](https://github.com/user-attachments/assets/d884f29d-2c41-4da7-a9be-4615edee1196)
Image source: https://doi.org/10.1371/journal.pbio.0030045

Recombinant proteins are proteins that wouldn't naturally appear in that organism. For instance, we can insert genes that code for fluorescence into plants, making them glow in the dark. Or, perhaps more usefully, we can take the human gene that codes for insulin and convince E. Coli to produce insulin instead. Insulin that we can then use to treat diabetic patients. Or proteins used in cancer treatment. Or, more recently, we can insert fragments of the SARS-CoV-2 virus into E. Coli and use that to produce COVID-19 vaccines. E. coli is a very well-studied organism, therefore, well-annotated. This means we can quickly check our work for any analysis we might perform because we have the ground truth, which biologists have spent decades meticulously gathering for us. We will examine the DNA sequence of E. coli and implement an algorithm for finding potential gene candidates. Because the ground truth is readily available, we can check how many of our ORF candidates are actual genes and how many candidates are false positives.

Every organism in NCBI has an associated unique identifier, which we can use to download various genomes and other kinds of data. The NCBI ID for E. coli is NC_000913.

First we need to initiate python. 
```bash
python
```
We are now in a python scripting environemnt. Not bash! Your bash commands will no longer work! It should looks something like this:
<img width="744" alt="Screenshot 2025-03-11 at 8 52 23 PM" src="https://github.com/user-attachments/assets/3ac28adb-4065-42d5-84e2-6fa10221dc60" />



Fetch the E. coli genome from NCBI. To check you have the right one, the genome should be 4,641,652 base pairs long. It might take up-to a minute to download.

We are using the Entrez.efetch (http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec149) and SeqIO.read (http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec37) functions from biopython. We are using rettype="gbwithparts" to fetch all required features.

```python
from Bio import Entrez

Entrez.email = "lyohe1@charlotte.edu"  # Always provide your email when using Entrez

# Fetch the E. coli K-12 MG1655 genome (NC_000913)
handle = Entrez.efetch(db="nucleotide", id="NC_000913", rettype="gbwithparts", retmode="text")

# Read and save the genome data
with open("ecoli_genome.gb", "w") as output_file:
    output_file.write(handle.read())

handle.close()
```

```python
from Bio import Entrez, SeqIO

# Always provide your email when using Entrez
Entrez.email = "lyohe1@charlotte.edu"

# Specify the genome you want to download (e.g., E. coli K-12)
genome_id = "NC_000913"

# Fetch the genome
handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")

# Save the FASTA file
with open(f"{genome_id}.fasta", "w") as output_file:
    output_file.write(handle.read())

handle.close()
```



[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
