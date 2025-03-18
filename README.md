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

SARS-CoV-2, the coronavirus causing COVID-19, will be analyzed in this lab series through bioinformatics approaches to understand its genetic mechanisms and potential therapeutic targets. We focus on identifying protein-coding genes within its genome, as these determine critical functions like host cell entry, replication, and immune evasion. By comparing SARS-CoV-2 genes with those of known viruses, we aim to infer functional roles and prioritize drug targets. Our methodology involves detecting open reading frames (ORFs) as candidate genes, filtering them using permutation tests to eliminate short non-functional sequences, and validating predictions through protein function classification algorithms. To ensure methodological rigor, we first test this pipeline on the well-characterized E. coli genome before applying it to SARS-CoV-2. This stepwise approach balances computational efficiency with biological relevance, streamlining the transition from gene prediction to functional annotation while minimizing false positives requiring experimental validation.

## Escherichia coli
Escherichia coli (E. Coli) is a bacteria commonly found in the human intestine. Most strains are harmless to humans, at worse, causing food poisoning and diarrhea. They can survive outside a host for only a short period of time, making it a potential indicator of fecal contamination. Over the years, E. Coli has been intensely studied and is probably one of the most well-understood organisms in existence. We've learned how to grow them in an optimal environment where they can reproduce up to once every 20 minutes. Due to their rapid growth and easy manipulation, biologists often use them to produce recombinant proteins.

Recombinant proteins are proteins that wouldn't naturally appear in that organism. For instance, we can insert genes that code for fluorescence into plants, making them glow in the dark. Or, perhaps more usefully, we can take the human gene that codes for insulin and convince E. Coli to produce insulin instead. Insulin that we can then use to treat diabetic patients. Or proteins used in cancer treatment. Or, more recently, we can insert fragments of the SARS-CoV-2 virus into E. Coli and use that to produce COVID-19 vaccines. E. coli is a very well-studied organism, therefore, well-annotated. This means we can quickly check our work for any analysis we might perform because we have the ground truth, which biologists have spent decades meticulously gathering for us. We will examine the DNA sequence of E. coli and implement an algorithm for finding potential gene candidates. Because the ground truth is readily available, we can check how many of our ORF candidates are actual genes and how many candidates are false positives.
![DFE-Ecoli-Kinds-2019-Pathogen-Illustration](https://github.com/user-attachments/assets/2c00a1e4-6f49-4b20-b2e3-4a48904a1510)




Every organism in NCBI has an associated unique identifier, which we can use to download various genomes and other kinds of data. The NCBI ID for E. coli is NC_000913.

First we need to initiate python. 
```bash
python
```
We are now in a python scripting environemnt. Not bash! Your bash commands will no longer work! It should looks something like this:
<img width="744" alt="Screenshot 2025-03-11 at 8 52 23 PM" src="https://github.com/user-attachments/assets/3ac28adb-4065-42d5-84e2-6fa10221dc60" />



Fetch the E. coli genome from NCBI. We already learned how to do this in bash. Now we are going to learn with biopython.  It might take up-to a minute to download.

We are using the Entrez.efetch (http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec149) and SeqIO.read (http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec37) functions from biopython.

Here are two ways to 
```python
from Bio import Entrez

Entrez.email = "your.email@charlotte.edu"  # Always provide your email when using Entrez

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
Entrez.email = "your.email@charlotte.edu"

# Specify the genome you want to download (e.g., E. coli K-12)
genome_id = "NC_000913"

# Fetch the genome
handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")

# Save the FASTA file
with open(f"{genome_id}.fasta", "w") as output_file:
    output_file.write(handle.read())

handle.close()
```

To exit python and return to bash, type the folloiwng command:

```python
exit()
```

## LQ 1.a: Reading python script

What are two differences in the script above?

## LQ 1.b Understanding the output

Inspect both outputs using head. How do they differ?

## LQ 1.c 
 Determine how many base pairs are in the E.coli genome.

## LQ 1.d
Paste the command you used to determin this or explain how you got your answer.

Okay, great!!! Moving right along. We need to get set up to do some more elegant analyses like looking for open reading frames. Make some functions accessible here.

```bash
#make sure you are in your lab_9 directory
python
```

I am not going to have you WRITE any python code yourself today, but we should begin to look at what it is doing. Here, I am defining three fucntions we will use:
```bash
find_orfs()
find_all_orfs()
```
You can copy and paste all of this code directly below into the python prompt. It will create a set of functions we can access in this session.

```python
from Bio import SeqIO
from Bio.Seq import Seq
from typing import List, Dict

def find_orfs(seq_record: SeqIO.SeqRecord, min_length: int = 30) -> List[Dict]:
    """
    Find all open reading frames (ORFs) in a sequence record
    
    Args:
        seq_record: Biopython SeqRecord object containing the sequence
        min_length: Minimum ORF length in nucleotides (default: 30)
        
    Returns:
        List of dictionaries with ORF information
    """
    orfs = []
    sequence = seq_record.seq
    seq_length = len(sequence)
    
    # Define genetic code parameters
    start_codons = ['ATG', 'GTG', 'TTG']
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    # Check all six reading frames (3 forward, 3 reverse)
    for frame in range(6):
        # Set frame parameters
        is_reverse = frame >= 3
        frame_offset = frame % 3
        
        # Get appropriate sequence for frame
        if is_reverse:
            working_seq = sequence.reverse_complement()
        else:
            working_seq = sequence
            
        # Adjust for frame offset
        working_seq = working_seq[frame_offset:]
        
        # Scan through the sequence
        pos = 0
        while pos + 3 <= len(working_seq):
            codon = str(working_seq[pos:pos+3])
            
            if codon in start_codons:
                start = pos
                end = pos + 3
                has_stop = False
                
                # Look for stop codon
                while end + 3 <= len(working_seq):
                    codon = str(working_seq[end:end+3])
                    if codon in stop_codons:
                        has_stop = True
                        end += 3
                        break
                    end += 3
                
                if has_stop:
                    orf_length = end - start
                    if orf_length >= min_length:
                        # Convert positions to original sequence coordinates
                        if is_reverse:
                            original_start = seq_length - (frame_offset + end)
                            original_end = seq_length - (frame_offset + start)
                        else:
                            original_start = frame_offset + start
                            original_end = frame_offset + end
                            
                        orfs.append({
                            'sequence_id': seq_record.id,
                            'start': original_start,
                            'end': original_end,
                            'length': orf_length,
                            'strand': '-' if is_reverse else '+',
                            'frame': frame + 1,
                            'protein': working_seq[start:end].translate()
                        })
                pos = end  # Skip past this ORF
            else:
                pos += 3  # Move to next codon
                
    return orfs

def find_all_orfs(genome_file: str, min_length: int = 30) -> List[Dict]:
    """
    Find all ORFs in a genome file (FASTA/GenBank format)
    
    Args:
        genome_file: Path to input file
        min_length: Minimum ORF length in nucleotides
        
    Returns:
        List of ORF dictionaries from all sequences in file
    """
    all_orfs = []
    
    # Auto-detect file format
    file_format = 'fasta' if genome_file.endswith(('.fasta', '.fa')) else 'genbank'
    
    with open(genome_file, 'r') as fh:
        for record in SeqIO.parse(fh, file_format):
            orfs = find_orfs(record, min_length)
            all_orfs.extend(orfs)
    
    return all_orfs

```

## LQ 9.2a

After analyzing the code, cosider the following as TRUE or FALSE.

The function find_all_orfs() implements find_orfs().

## LQ 9.2b

What are the differences between the two functions find_all_orfs() and find_orfs()?


That was fun, wasn't it? Okay, do not close your python terminal. Python is an object-oriented programming language. We are going to create first "object" here that will inherit the output of the functions we just created. 

```python
orfs = find_all_orfs("NC_000913.fasta")
```
It may take a second to run.
## LQ 9.2c

We made an object that inherited the output of find_all_orfs(). What did we name this object?

Let's take some time to inspect this object. If we look at the documentation, the output of find_all_orfs() should give us a "List of ORF dictionaries from all sequences in file." This is super exciting becauase we essentially are able to carry all of our wanted output in a single variable.

If you just type 
```python
orfs
```
a bunch of stuff prints to screen. 
```python
type(orfs)
```

## LQ 9.2d

What type of object is our object?

Now we are going to make an object of the first item in our object! Notice how just typing what you want and typing your object name give you the same answer.
```python
[orfs[0]]
first = [orfs[0]]
first
```
The object first inherits all of the output of our command, saving us lots of typing and enhancing our ability to scale up. Cool.

## LQ 9.2e
What does the second line in our orfs[] list look like? Find its length and submit the answer.


Okay, now we are going to clean things up a little bit. You can see the first and second ORF in our orfs[] list are very different lengths. Let's only select the ORFs longer than 300bp.

```python
#here is a cute mini for loop function for you!!! :) 
long_orfs = [orf for orf in orfs if orf['length'] >= 300]
print(f"Found {len(long_orfs)} ORFs longer than 300bp")
```

## LQ 9.3a
How many orfs are longer than 300bp?

Alter the print function slightly! 

## LQ 9.3b
How many orfs did we originally have before we trimmed?

# PART 2: Working with a novel genome
Okay, now time for things to get interesting. We are ready to work with the SARS-CoV-2 genome!
![sars_cov_2_shutterstock_1889563108_crop_0](https://github.com/user-attachments/assets/e19e8c07-e6de-4040-a235-bb8dddb0aff1)
I am sure most of you are well aware of what the disease does, so we'll skip any long-winded introduction. Instead, we'll jump straight into the genome of this pesky little virus.



Exit out of python.
```python
exit()
```

I have kindly provided you with the genome in the shared lab folder. Bring it over to your lab_9 directory.
```
cp /projects/class/binf3101_001/sars_cov_2.fasta ~/lab_9
```


Let's see how long the virus of SARS-CoV-2 is.
```bash
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'   sars_cov_2.fasta
```

# LQ9.4

The genome of SARS-CoV-2 is ________ bp and it is ________ bp smaller than the ecoli genome. Tiny! (hint: use substraction from questsion LQ9.1c)


Okay, now we are going to do the same thing to the SARS-CoV-2 genome as we did with E.coli. We need to find all of the ORFs and filter out the ones less than 300bp. Don't be scared, you are on your way to being a python whiz.

Get your python envinrment set up in case you are starting from a new log in.

```bash
module load python
python
```
We will have to redefine the functions again, but there is one key difference.
The start codon and stop codons of SARS-CoV-2 differ sligthly than E.Coli. Change the two lines of code from the functions we made from E.coli before copying and pasting.
```
start_codons = ["ATG"]
stop_codons = ["TGA", "TAA", "TAG"]
```

# LQ 9.5a
SARS-CoV-2 has _______ total ORFs and ____ ORFs less than 300 bp.

# LQ 9.5b
The length of the first ORF in the SARS-CoV-2 genome from the TRIMMED ORFs object (long_orfs) is ____ bp.

Quite different than our good friend E.Coli! Nothing but bare-bones in most viruses. :) Wild how they can wreak havoc with just small levels of genomic machinery.



# PART 3: Determining gene functionality from amino acid sequence
Okay, so we've found and filtered down our gene candidates, but now what? The next step is to figure out what these genes might do. So, what do these genes do? Well, nothing by themselves. If DNA is the cookbook, then the genes are the recipes. They tell us how to prepare each dish. Recipes by themselves are nothing but information. It is the actual dishes we care about: the proteins. It is the proteins that actually do things inside our cells.

An organism can read codons in ORFs (hence open-READING-frame) and translate the language of codons into amino acids. Codons in a DNA sequence are translated into amino acids consecutively, forming a long chain. The amino acid chain folds into a macromolecule called protein (as shown in the picture below). Their 3D structure is crucial for their function. We can predict protein characteristics from the amino acid sequence, such as their location in a cell. However, for a more concrete prediction of the functionality, we would have to use more sophisticated approaches, such as BLAST or AlphaFold.


## Intro to hydrophobicity
The hydrophobicity of an amino acid is its tendency to repel water molecules. These amino acids tend to stick together in a non-polar environment. One well-known example of hydrophobicity is the lipid membrane of cells. This hydrophobic cell membrane ensures that the environment outside the cell (water and anything inside that water) stays outside, and all the important bits of the cell (e.g., the mitochondria, nucleus, ribosomes) stay inside the cell.

Hydrophilic proteins cannot pass through the hydrophobic membrane since its hydrophobic interior will repel them, but hydrophobic proteins can sit happily embedded inside the membrane. These hydrophobic proteins are called transmembrane proteins. Transmembrane proteins typically function as transport channels between the interior and exterior of the cell and enable cell signaling and absorption. In the case of viruses, e.g., SARS-CoV-2, we expect at least some proteins to be embedded in the viral membrane. For instance, SARS-CoV-2 viruses can only enter human cells when the spike protein that sticks out on the exterior of the virus binds to a specific receptor in human cells.

Transmembrane proteins usually contain more hydrophobic amino acids than other proteins because their environment is a strictly hydrophobic lipid membrane. Therefore, one naive approach to finding transmembrane proteins might be to look for proteins with a high degree of hydrophobicity.

The less-biological explanation of all this might be that the little, yellow hydrophobic tails of the lipids in the membrane attract each other because they are all hydrophobic. They can also attract proteins, provided that the amino acids inside that protein are also hydrophobic. The little hydrophobic tails also repel anything hydrophilic, so hydrophilic proteins can't pass through the hydrophobic membrane. So, simply put, hydrophobic molecules attract other hydrophobic molecules and repel any hydrophilic molecules. Finally, hydrophobic molecules repel water, while hydrophilic molecules are drawn toward water.

Yada, yada, yada.

First, some housekeeping. In your bash terminal, still in the lab_9 directory, we need to install a few more modules and copy over a table with the hydrophobicity values of each amino acid.
```bash
pip install seaborn
cp /projects/class/binf3101_001/aminoacid_properties.csv ~/lab_9
```

We now are going to translate our open reading frames and export them into a fasta file of amino acid sequences using this beautiful code I am providing you.

```python
from typing import List, Dict

def export_trimmed_proteins(orfs: List[Dict], output_filename: str) -> None:
    """
    Export ORF protein sequences in FASTA format with stop codons removed
    
    Args:
        orfs: List of ORF dictionaries (from find_orfs/find_all_orfs)
        output_filename: Name for output FASTA file
    
    Returns:
        None (writes to file)
    """
    with open(output_filename, 'w') as fasta_file:
        for i, orf in enumerate(orfs, 1):
            # Trim trailing stop codon and convert to string
            protein_seq = str(orf['protein']).rstrip('*')
            
            # Create informative FASTA header
            header = (
                f">{orf['sequence_id']}_ORF{i}|"
                f"start={orf['start']}|end={orf['end']}|"
                f"strand={orf['strand']}|length={orf['length']}aa"
            )
            
            # Write formatted sequence (60 chars per line)
            fasta_file.write(header + "\n")
            for j in range(0, len(protein_seq), 60):
                fasta_file.write(protein_seq[j:j+60] + "\n")
```

You can run the command on your long_orfs you saved from SARS-CoV-2. 

```python
export_trimmed_proteins(long_orfs, "sars_cov_2_proteome.fasta")
```

I know it may be frustrating for you to move from python to bash. If you want to inspect your file names in your directory, you can use the following commands. Just note it will return them as a list, but it will still allow you to see what is there.


```python
import os
os.listdir()
```
Now, here is where the magic happens. 
Here is a script to 
```python
import csv
from Bio import SeqIO
from typing import List, Dict
import re

def load_hydrophobicity_table(csv_path: str) -> Dict[str, float]:
    """
    Load amino acid hydrophobicity values from CSV file
    
    Args:
        csv_path: Path to CSV file with columns 'amino_acid' and 'hydrophobicity'
        
    Returns:
        Dictionary mapping single-letter amino acid codes to hydrophobicity values
    """
    hydro_dict = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            aa = row['aminoacid'].upper()
            hydro_dict[aa] = float(row['hydrophobicity'])
    return hydro_dict

def calculate_hydrophobicities(fasta_file: str, hydro_table: Dict[str, float]) -> List[Dict]:
    """
    Calculate average hydrophobicity for proteins in FASTA file
    
    Args:
        fasta_file: Path to FASTA file from export_trimmed_proteins()
        hydro_table: Hydrophobicity dictionary from load_hydrophobicity_table()
        
    Returns:
        List of dictionaries with protein metadata and hydrophobicity scores
    """
    results = []
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        # Parse metadata from FASTA header
        header_parts = re.split(r'\||=', record.description)
        
        # Extract sequence information
        seq_id_orf = header_parts[0].lstrip('>')
        seq_parts = seq_id_orf.split('_ORF')
        sequence_id = seq_parts[0]
        orf_number = int(seq_parts[1]) if len(seq_parts) > 1 else 1
        
        metadata = {
            'sequence_id': sequence_id,
            'orf_number': orf_number,
            'start': int(header_parts[2]),
            'end': int(header_parts[4]),
            'strand': header_parts[6],
            'aa_length': int(header_parts[8].replace('aa', '')),
            'protein_id': record.id
        }
        
        # Calculate hydrophobicity
        total = 0.0
        valid_aas = 0
        sequence = str(record.seq).upper()
        
        for aa in sequence:
            if aa in hydro_table:
                total += hydro_table[aa]
                valid_aas += 1
            else:
                print(f"Warning: No hydrophobicity value for '{aa}' in {record.id}")
        
        avg_hydro = total / valid_aas if valid_aas > 0 else float('nan')
        
        results.append({
            **metadata,
            'average_hydrophobicity': avg_hydro,
            'valid_aa_count': valid_aas,
            'total_aa_count': len(sequence)
        })
    
    return results

```

Run the command

```python
 hydro_table = load_hydrophobicity_table("aminoacid_properties.csv")
```

[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
