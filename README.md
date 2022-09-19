# Transcriptome assembly pipeline

This pipeline, heavily adapted from that created by Yang Lab [[1]](#1), assembles *de novo* transcriptomes from raw SRA data retrievable from the NCBI database. As an extended version of Yang's, this pipeline enables automatic retrieval and assembly of RNA-seq data, and is more adapted to performing bulk-assemblies in fewer steps. Through the wrapping and use of several third-party packages, all of which are listed below, this pipeline performs all pre-processing, filtration, and QC steps necessary for the bulk assembly of ready-to-use transcriptomes.

## Dependencies
**NOTE: for each package listed, the associated version used in the pipeline is given in parentheses. I am not sure if other / newer versions of the packages will work.**

- BLAST+ (v2.13.0)
- Corset (v1.09)
- FastQC (v0.11.9)
- Rcorrector (v1.0.5)
- Salmon (v1.9.0)
- SRA Toolkit (v3.0.0)
- TransDecoder (v5.3.0)
- Trimmomatic (v0.39)
- Trinity (v2.14.0)

## Workflow

[1. Retrieval and read processing](#read_processing)  
[2. Assembly with Trinity](#trinity)  
[3. Transcript filtering and translation](#filt_trans)  


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
By default, file outputs from intermediary processing steps will not be removed. To remove them automatically during processing, call the script with the **rem-inter** flag:
```
python3 process_reads num_threads rem-inter
```
This will perform all the above-listed processing steps on all SRA runs specified in the **sras.txt** file. Once completed, the processed transcript files can be found in the **"05-filter_over_represented"** folder in the current working directory.

### 2. Assembly with Trinity <a name="trinity"></a>

Following initial processing of the reads, Trinity will be used to generate *de novo* transcriptomes for each SRA run. This process is handled by the **assemble_reads.py** script, which takes as its arguments the number of threads and the maximum amount of RAM (in gigabytes) that may be dedicated towards the Trinity assemblies. To perform this step, simply run the script:

To assemble from a single SRA run:
```
python3 assemble_reads.py num_threads max_memory_GB
```
To assemble from several SRA runs (such as from multiple differing tissue samples):
```
python3 assembly_reads.py num_threads max_memory_GB mult_samples
```
This will initiate Trinity assembly for all processed transcripts in series. Note that this step can be quite time consuming, depending on the size and complexity of the dataset in question.

Once complete, all output assemblies can be found in the **"06-trinity_assembly"** folder in the current working directory.

### 3. Transcript filtering and translation <a name="filt_trans"></a>

The following steps are performed during transcript filtering and translation:

1. Chimera detection and removal with BLASTX
2. Transcript clustering with Corset
3. Filtering of Corset clusters
4. Translation with Transdecoder

With our newly-assembled Trinity transcripts, the next stage is to filter and translate them. All above steps are handled by the **post_assembly.py** script. The first step here utilizes BLAST along with a reference proteome to identify and remove chimeras from the Trinity transcripts. For the reference proteome, choose several species that are closely related to those being fed through the pipeline. For instance, when passing in several plants from the Caryophyllales family, one may assemble the reference proteome from spinach, beets, and arabidopsis proteomes.

For simple construction of the reference proteome, a script called **concat_fasta.py** has been provided, and can be executed thus:
```
python3 concat_files.py [proteome_1.fasta, proteome_2.fasta, ...] output_directory output_file.fasta
```
With a reference proteome ready, simply run the **post_assembly.py** script:

To process assemblies from single SRA runs:
```
python3 post_assembly.py proteome_reference.fasta num_threads
```
To process assemblies from several SRA runs:
```
python3 post_assembly.py proteome_reference.fasta num_threads
```
This command will handle all the post-processing steps listed above and generate the final coding sequences output for each assembly, all of which can be found in the **"09-translate"** directory upon completion.

## References
<a id= "1">[1]</a>  Yang Lab (2021) Phylogenomic Dataset Construction. https://bitbucket.org/yanglab/phylogenomic_dataset_construction.git


