# pbb-drylab
dry lab methods used by [@PuckerLab](https://www.tu-braunschweig.de/en/ifp/pbb)

WARNING: Some of this content is not publicly available. If you are not a member of the PuckerLab, you might not be able to use it.

## Introduction 

All scripts are available in /grp/pbb/scripts. Additional tools can be found in /grp/pbb/tools. Each member is encouraged to have a folder in /grp/pbb/members to make scripts available to the group. Please make sure that no large files are shared through this group volumen, because the size is rather small. Dedicated project volumes should be used for the storage and exchange of larger files.

Please do not use normal terminals for heavy computing. Only terminals on dedicated compute nodes should be used for this purpose.


## Scripts (developed by group)

### automatic_SRA_download.py

### BWA_MEM_wrapper.py

### collect_best_BLAST_hits.py

### construct_DESeq2_input.py (requires modification!)

### DGE_analysis.R

### extract_myb_domain.py

### filter_RNAseq_samples.py

### get_DEGs_from_DESeq2_output.py

### kallisto_pipeline.py

### mask_monophyly_bp.py

### merge_kallisto_output.py

### reads2counts2.py

### rename_FASTA_seqs.py

### tree.py

### tree2.py




## Standard tools (developed by others)

A large number of tools is installed globally thus you just have to type the name of your tool into the terminal to use it.

### FastQC / MultiQC

### BWA-MEM

### FastTree

### IGV

### Picard-tools

### RAxML

### Samtools

### SRAtoolkit (fastq-dump)

### STAR

### Trimmomatic

### Trinity


### Single cell RNA-seq (scRNA-seq)
CellRanger for quantification of reads. 

There is an excellent book about the following analysis: [Orchestrating Single Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/)



