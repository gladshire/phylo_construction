# Transcriptome assembly pipeline

***NOTE: This pipeline is not finished. As of now, it can perform all assembly steps up to and including chimera detection and removal. Downstream steps, such as clustering, tree construction, and tree analysis are not yet implemented***

This pipeline, heavily adapted from that created by Yang Lab [[1]](#1), assembles *de novo* transcriptomes from raw SRA data retrievable from the NCBI database. As an extended version of Yang's, this pipeline enables automatic retrieval and assembly of RNA-seq data, and is more adapted to performing bulk-assemblies in fewer steps. Through the wrapping and use of several third-party packages, all of which are listed below, this pipeline performs all pre-processing, filtration, and QC steps necessary for the bulk assembly of ready-to-use transcriptomes.

## Dependencies
**NOTE: for each package listed, the associated version used in the pipeline is given in parentheses. I am not sure if other / newer versions of the packages will work.**

- BLAST+ (v2.13.0)
- cd-hit (v4.8.1)
- Corset (v1.09)
- FastTree (v2.1.11)
- FastQC (v0.11.9)
- Gblocks (v0.91b)
- mafft (v7.490)
- MCL (v14-137)
- Pasta (v1.9.0)
- Phyx (v1.3)
- Prank (v0.170427)
- RAxML (v8.2.12)
- Rcorrector (v1.0.5)
- Salmon (v0.9.1)
- SRA Toolkit (v3.0.0)
- TransDecoder (v5.3.0)
- TreeShrink (v1.3.9)
- Trimmomatic (v0.39)
- Trinity (v2.14.0)

## Workflow

[1. Retrieval and read processing](#read_processing)  
[2. Assembly with Trinity](#trinity)  
[3. Transcript filtering and translation](#filt_trans)  
[4. Clustering](#clustering)  
[5. Build homolog trees](#homo_trees)  
[6. Paralogy pruning](#para_prune)  
[7. Construct supermatrix](#supermatrix)  


### 1. Retrieval and read processing <a name="read_processing"></a>

The following steps are performed during read processing:

1. Data retrieval with prefetch and fasterq-dump
2. Error correction with Rcorrector
3. Removal of unfixable errors found in (2)
4. Removal of adaptors with Trimmomatic
5. Filtration of foreign reads with Kraken2
6. Quality analysis with FastQC
7. Removal of over-represented reads found in (6)

Prior to execution of the read processing script, locate **sras.txt** in the scripts folder. This is where the SRA numbers for RNA-seq data retrieval are specified. In the text file, below the line which reads: 
```
Write SRA codes for retrieval below this line
```
Specify the SRA numbers for each run to be retrieved and processed by the pipeline, line-wise and without whitespace. Once all desired SRA numbers are entered, exit and run the **process_reads.py** python script, specifying the desired number of threads to dedicate:
```
python3 process_reads.py num_threads
```
This will perform all the above-listed processing steps on all SRA runs specified in the **sras.txt** file. Once completed, the processed transcript files can be found in the **"05-filter_over_represented"** folder in the current working directory.

### 2. Assembly with Trinity <a name="trinity"></a>

Following initial processing of the reads, Trinity will be used to generate *de novo* transcriptomes for each SRA run. This process is handled by the **assemble_reads.py** script, which takes as its arguments the number of threads and the maximum amount of RAM (in gigabytes) that may be dedicated towards the Trinity assemblies. To perform this step, simply run the script:
```
python3 assemble_reads.py num_threads max_memory_GB
```
This will initiate Trinity assembly for all processed transcripts in series. Note that this step can be quite time consuming, depending on the size of the dataset in question, and will likely be the longest step in the pipeline's total runtime.

Once complete, all output assemblies can be found in the **"06-trinity_assembly"** folder in the current working directory.

### 3. Transcript filtering and translation <a name="filt_trans"></a>

### 4. Clustering <a name="clustering"></a>

### 5. Build homolog trees <a name="homo_trees"></a>

### 6. Paralogy pruning <a name="para_prune"></a>

### 7. Construct supermatrix <a name="supermatrix"></a>

## References
<a id= "1">[1]</a>  Yang Lab (2021) Phylogenomic Dataset Construction. https://bitbucket.org/yanglab/phylogenomic_dataset_construction.git


