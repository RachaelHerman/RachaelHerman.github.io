progressiveMauve --threads=90 \
--log=progressiveMauve.log \
--output=synteny_alignment_pygoscelis.xmfa \
--backbone-output=synteny_backbone.txt \
/gpfs/scratch/rherman/pacbio/SYNTENY/adelie_pacbio_ref/ADPE_BROW_1.asm.bp.p_ctg.fasta \
/gpfs/scratch/rherman/pacbio/SYNTENY/gentoo_pacbio_ref/GEPE-SCRY-2.asm.bp.p_ctg.fasta


nucmer --maxmatch /gpfs/scratch/rherman/pacbio/SYNTENY/adelie_pacbio_ref/ADPE_BROW_1.asm.bp.p_ctg.fasta /gpfs/scratch/rherman/pacbio/SYNTENY/gentoo_pacbio_ref/GEPE-SCRY-2.asm.bp.p_ctg.fasta \ delta-filter -r -q out.delta > out.filtered.delta
mummerplot out.filtered.delta




[rherman@dg-mem SYNTENY]$ bash progressive_mauve.sh
Storing raw sequence at /tmp/rawseq3905387.000
Sequence loaded successfully.
/gpfs/scratch/rherman/pacbio/SYNTENY/Megadyptes_genome/GCA_010078485.1_data/GCA_010078485.1/GCA_010078485.1_BGI_Mant.V1_genomic.fna 1314943623 base pairs.                                           
Storing raw sequence at /tmp/rawseq3905387.001
Sequence loaded successfully.
/gpfs/scratch/rherman/pacbio/SYNTENY/adelie_pacbio_ref/ADPE_BROW_1.asm.bp.p_ctg.fasta 1430289379 base pairs.                                                                                         
Storing raw sequence at /tmp/rawseq3905387.002
Sequence loaded successfully.
/gpfs/scratch/rherman/pacbio/SYNTENY/gentoo_pacbio_ref/GEPE-SCRY-2.asm.bp.p_ctg.fasta 1433789152 base pairs.                                                                                         
Using weight 21 mers for initial seeds
Sorted mer list loaded successfully
Sorted mer list loaded successfully
Sorted mer list loaded successfully
0%..1%..2%..3%..4%..5%..6%..7%..8%..9%..10%..
11%..12%..13%..14%..15%..16%..17%..18%..19%..20%..
21%..22%..23%..24%..25%..26%..27%..28%..29%..30%..
31%..32%..33%..34%..35%..36%..37%..38%..39%..40%..
41%..42%..43%..44%..45%..46%..47%..48%..49%..50%..
51%..52%..53%..54%..55%..56%..57%..58%..59%..60%..
61%..62%..63%..64%..65%..66%..67%..68%..69%..70%..
71%..72%..73%..74%..75%..76%..77%..78%..79%..80%..
81%..82%..83%..84%..85%..86%..87%..88%..89%..90%..
91%..92%..93%..94%..95%..96%..97%..98%..99%..100%..
done.
using default bp penalty: 212629
using default bp estimate min score: 637887
Starting with 49318877 multi-matches
Computing genome content distance matrix...


Genome conservation distance matrix:
0       0.305006        0.29555
0.305006        0       0.236544
0.29555 0.236544        0

Writing guide tree to /tmp/guide_tree3905387.000
reading tree...
initializing alignment tree...
Constructing seed occurrence lists for repeat detection
Calculating pairwise breakpoint distances
Pair 0, 1 has 2162176 initial LCBs
Using scaled bp penalty: 148823
Pair (0,1) has 8592 well-supported breakpoints
Pair 0, 2 has 2168190 initial LCBs
Using scaled bp penalty: 157088
Pair (0,2) has 9273 well-supported breakpoints

10%..11%..12%..13%..14%..15%..16%..17%..18%..19%..
20%..21%..22%..23%..24%..25%..26%..27%..28%..29%..
30%..31%..32%..33%..34%..35%..36%..37%..38%..39%..
40%..41%..42%..43%..44%..45%..done
Arrived at 1138 intervals
Adding unaligned intervals
addUnalignedIntervals yields 3373 intervals
Merging unaligned intervals
Marbling gaps
Propagating descendant breakpoints
descendant 0(2) has 1 intervals
descendant 1(4) has 1 intervals
propagateDescendantBreakpoints yields 1138 intervals
Creating ancestral ordering
Previous anchoring score: -1.79769e+308, new anchor score: 9.29869e+10
Backing up alignment tree...
propagating ancestral breakpoints
recursive anchor search
0,0 have 48521 new matches outside LCBs
0,0 has an additional 263715 matches
Restoring backed up alignment tree...
1,2 has 263715 pairwise matches
Performing Sum-of-pairs Greedy Breakpoint Elimination
construct LCB tracking matches
There are 11588055 tracking matches
There are 23176110 / 57940275 components used
init tracking match LCB tracking
pairwise score tracking matches
get pairwise LCBs
there are 77285 pairwise LCBs
scaling bp penalty by conservation weight:
0.236544


scaling bp penalty by bp weight:
0.0261881

Greedy BPE
Scoring with scaled breakpoint penalty: 68503
1%..2%..3%..4%..5%..6%..7%..8%..9%..
10%..11%..12%..13%..14%..15%..16%..17%..18%..19%..
20%..21%..22%..23%..24%..25%..26%..27%..28%..29%..
30%..31%..32%..33%..34%..35%..36%..37%..38%..39%..
40%..41%..42%..43%..44%..45%..46%..47%..48%..49%..
50%..51%..52%..53%..54%..55%..56%..57%..58%..59%..
60%..61%..62%..63%..64%..65%..66%..67%..68%..69%..
70%..71%..72%..73%..74%..75%..76%..77%..78%..79%..
80%..81%..82%..83%..84%..85%..86%..87%..88%..89%..
90%..done
Arrived at 1044 intervals
Adding unaligned intervals
addUnalignedIntervals yields 3027 intervals
Merging unal

Creating ancestral ordering
Previous anchoring score: 9.29869e+10, new anchor score: 8.95594e+10
propagating ancestral breakpoints
performing a gapped alignment


#!/usr/bin/env python

import random

# Define a list of colors (Circos-compatible)
colors = ["red", "blue", "green", "purple", "orange", "cyan", "magenta", "yellow", "black"]

# Read the input file
with open("circos_input.txt", "r") as infile, open("links.txt", "w") as outfile:
    for line in infile:
        if line.strip():
            color = random.choice(colors)  # Pick a random color
            outfile.write(f"{line.strip()} color={color}\n")

print("Formatted Circos input saved to links.txt")

#! /bin/bash

# Assuming you have genome FASTA files with .fai index files already generated
for genome in GCA_010078485.1_BGI_Mant.V1_genomic.fna ADPE_BROW_1.asm.bp.p_ctg.fasta GEPE-SCRY-2.asm.bp.p_ctg.fasta; do
    # Extract chromosome names and lengths from the .fai file
    while read chr len; do
        # Format the output for Circos (adjust color if needed)
        echo "chr - ${genome}_${chr} ${chr} 0 ${len} ${RANDOM}"
    done < ${genome}.fai
done > karyotype.txt


#!/bin/bash

# Define a list of colors for the three genomes
colors=("blue" "green" "purple")

# Define the genome names (replace these with your actual genome names)
genomes=("GCA_010078485.1_BGI_Mant.V1_genomic.fna" "ADPE_BROW_1.asm.bp.p_ctg.fasta" "GEPE-SCRY-2.asm.bp.p_ctg.fasta")

# Initialize a counter for genomes
genome_index=0

# Loop through each genome
for genome in "${genomes[@]}"; do
    # Extract chromosome names and lengths from the .fai file
    while read chr len; do
        # Get the color for the current genome
        color=${colors[$genome_index]}

        # Format the output for Circos with the color
        echo "chr - ${genome}_${chr} ${chr} 0 ${len} ${color}"
    done < ${genome}.fai
    
    # Increment the genome index and move to the next color
    genome_index=$((genome_index + 1))
done > karyotype.txt


awk '{ $7=""; $8=""; $9=""; print $1, $2, $3, $4, $5, $6, $10 }' karyotype.txt > karyotype_cleaned.txt



import sys

def convert_backbone_to_circos(backbone_file, output_file):
    with open(backbone_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            cols = line.strip().split()
            
            if len(cols) != 6:
                continue  # Skip malformed lines
            
            g1_start, g1_end = int(cols[0]), int(cols[1])
            g2_start, g2_end = int(cols[2]), int(cols[3])
            g3_start, g3_end = int(cols[4]), int(cols[5])
            
            # Assign chromosome names
            chr1, chr2, chr3 = "chr1", "chr2", "chr3"
            
            # Generate links, handling reverse strand cases
            if g1_start < 0:
                g1_start, g1_end = abs(g1_end), abs(g1_start)
            if g2_start < 0:
                g2_start, g2_end = abs(g2_end), abs(g2_start)
            if g3_start < 0:
                g3_start, g3_end = abs(g3_end), abs(g3_start)

            # Link Genome 1 ↔ Genome 2
            outfile.write(f"{chr1} {g1_start} {g1_end} {chr2} {g2_start} {g2_end} color=blue\n")
            
            # Link Genome 1 ↔ Genome 3 (if valid)
            if g3_start != 0 and g3_end != 0:
                outfile.write(f"{chr1} {g1_start} {g1_end} {chr3} {g3_start} {g3_end} color=green\n")
            
            # Link Genome 2 ↔ Genome 3 (if valid)
            if g3_start != 0 and g3_end != 0:
                outfile.write(f"{chr2} {g2_start} {g2_end} {chr3} {g3_start} {g3_end} color=red\n")

if __name__ == "__main__":
    backbone_file = "alignment.xmfa.backbone"  # Change this to your actual file
    output_file = "circos_input.txt"
    convert_backbone_to_circos(backbone_file, output_file)
    print(f"Circos link file saved as {output_file}")
