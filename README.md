# StORF-Reporter

## Overview
This project focuses on enhancing genome annotation by identifying and characterising stop-to-stop open reading frames (StORFs) within bacterial genomes.
Traditional genome annotation methods often overlook genes due to stringent criteria such as specific start codons or minimum length thresholds. 
StORF-Reporter addresses these limitations by extracting and analysing previously unannotated genomic regions (URs) to detect potential genes, 
including those with alternative start codons,small genes, and overlapping genes.


## Goals

- Identify novel StORFs and conserved small open reading frames (con-stORFs) from bacterial genomes.

- Assess and refine criteria for accurately reporting the longest and most biologically relevant StORFs.

- Improve genome annotation quality by integrating previously missed genetic elements.


## Aims

- Analyse the results to determine if the longest StORF should always be selected in an unannotated region.

- Analyse the results to determine if a reported con-stORF overlaps with a reported StORF.

- Determine better rules to filter StORFs so that the sequence most likely to be a gene is reported.


## Indicative Tasks

- Run StORF Reporter on a bacterial genome to identify StORFs and con-stORFs.

- Evaluate if the longest StORF consistently represents the biologically most relevant gene.

- Assess overlap between reported con-stORFs and StORFs to identify potential annotation issues.
- Develop and test improved criteria for filtering StORFs.



## Features and Tools

#### UR-Extractor

- Extracts unannotated regions from genomic FASTA files using corresponding GFF annotations.
- Generates separate FASTA and GFF files for unannotated regions, facilitating targeted analysis.

#### StORF-Finder

- Scans unannotated regions to identify potential StORFs based on stop-to-stop codon sequences.
- Produces annotated FASTA and
 GFF files containing candidate StORFs.

#### StORF-Reporter

- Automates the complete workflow from UR extraction to StORF identification and validation.
- Int
egrates external tools such as Pyrodigal for comprehensive gene prediction.

#### StORF-Extractor

- Allows filtering and extraction of specific subsets of StORFs based on user-defined criteria (e.g., length, genomic location).

#### StORF-Remover

- Refines StORF datasets by removing sequences that do not meet specified validation criteria.



## Analysis Workflow
1. Data Preparation: Extract URs from annotated genomes.

2. StORF Identification: Use StORF-Finder to scan extracted URs.
3. Validation: Conduct BLAST analyses to verify the biological relevance of detected StORFs and con-stORFs.
4. Clustering and Annotation: Cluster identified StORFs using tools like CD-HIT, annotate novel genetic elements, and integrate them into existing genome datasets.
5. Result Interpretation: Evaluate whether the longest identified StORF is consistently the most conserved and biologically relevant.


## Usage

#### UR Extraction
UR-Extractor -f genome.fna -gff annotations.gff -odir output_directory -gz True

#### StORF Finding
StORF-Finder -f genomic_UR.fasta -gff genomic_UR.gff -odir output_directory -aa True -gff True -gz True

#### Automated Pipeline
StORF-Reporter -anno Pyrodigal -p genomic_UR.fasta -odir output_directory -aa True -sout True

#### StORF Extraction
StORF-Extractor -f storfs.fasta -min_len 30 -max_len 500 -odir output_directory

#### Filtering Unwanted StORFs
StORF-Remover -f storfs.fasta -min_len 50 -odir output_directory


## Dependencies

- Python 3.x
- Pandas, NumPy, Matplotlib, Seaborn
- BLAST and DIAMOND for sequence alignment
- CD-HIT for clustering analysis

