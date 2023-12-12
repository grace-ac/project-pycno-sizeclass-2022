# Summer 2022 _Pycnopodia helianthoides_ RNAseq Analyses Pipeline

# RNAseq Sample Info:

32 libraries: [/grace-ac/project-pycno-sizeclass-2022/data/summer2022_samples_sequenced.csv](https://github.com/grace-ac/project-pycno-sizeclass-2022/blob/main/data/summer2022_samples_sequenced.csv)

More information on the project in main repository README.md: [/project-pycno-sizeclass-2022/README.md](https://github.com/grace-ac/project-pycno-sizeclass-2022/blob/main/README.md).

# RNAseq Workflow:     

1. Data Management
2. QC and Trimming
3. `HISAT2`: Alignment to _Pycnopodia helianthoides_ genome
4. Count Matrices with `kallisto` (or `HTseq`? something else? `StringTie`?)
5. `DESeq2`: Differential Gene Expression Analyses

# 1. Data Management
Received data from Azenta (FKA Genewiz).

Downloaded to genefish using: [cyberduck](https://cyberduck.io/download/).     
Azenta provided log-in credentials and information of how to access their servers. Azenta pdf of different ways to download the data from their servers: [here](https://f.hubspotusercontent00.net/hubfs/3478602/Sell%20Sheet%20Collateral%20Library/NGS/NGS%20User%20Guides/NGS_sFTP-Data-Download-Guide_Option%201_Nov03_2020.pdf).

All downloaded files were put into a folder in Documents on genefish called: "P_helianthoides_RNAseq", with a subfolder for Summer 2022 data is in folder: "30-833328413".

## Move data to OWL     
### A. `rsync` data to `owl/nightingales/P_helianthoides`
On genefish, in commandline in the directory where the data lives:
```
rsync --archive --progress --verbose PSC*.fastq.gz <owl username>@128.95.149.83:/volume1/web/nightingales/P_helianthoides/
```
Replace <owl_username_> with whatever username you use to login to owl (even replace the < and the >).

### B. `rsync` checksums to OWL on commandline:     
#### B.1. `ssh` into OWL:     
```
ssh username@owl.fish.washington.edu
```
Will be prompted for password.    

#### B.2. Navigate to folder you want to work in, in my case: `owl/nightingales/P_helianthoides`

#### B.3. Put below code in and hit ENTER:     
```
for fastq in PSC*.fastq.gz
do
  md5sum "${fastq}" >> checksums.md5
  echo "Generated checksum for ${file}."
  echo ""
done
```

All RNAseq data and checksums are now in [owl/nightingales/P_Helianthoides](https://owl.fish.washington.edu/nightingales/P_helianthoides/).

# 2. QC and Trimming
## 1. Untrimmed Data QC Part I: `FASTQC`  
### A. Get FastQC if you want to run on your laptop:    
[https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)     

### B. Get .fastq.gz files from OWL onto RAVEN   
#### i. `ssh` into Raven using credentials in command line.    
#### ii. Make a directory for all PSC.fastq.gz files (made one called `pycnornaseq2022` in my `/home/shared/8TB_HDD_02/graceac9/pycnornaseq2022`)  

### C. Then move files from OWL to pycnornaseq directory in Raven:

#### a. Have Husky OnNet App (BIG-IP Edge Client in Applications folder after downloaded)
#### b. Log in with UW credentials
#### c. Put RStudio IP into browser: http://172.25.149.12:8787
#### d. Log in using Raven Credentials
#### e. cd into pycnornaseq/ and run:
```
wget -r --no-directories --no-parent  -A "PSC-0*" https://owl.fish.washington.edu/nightingales/P_helianthoides
```

Then the directory will look like this:     
<img width="573" alt="PSC_files_on_RAVEN" src="https://github.com/grace-ac/project-pycno-sizeclass-2022/blob/main/protocols/images/screenshot-1.png">

### D. Get into Rstudio on Raven to run FASTQC:  
Follow steps a-d in Step C listed above.

Then, follow the code outline in this script: [project-pycno-sizeclass-2022/code/01-FastQC_pre-trim.Rmd](https://github.com/grace-ac/project-pycno-sizeclass-2022/blob/main/code/01-FastQC_pre-trim.Rmd)

The FASTQC files are saved on Raven: `/home/shared/8TB_HDD_02/graceac9/analyses/pycno2022`.

## 2. Untrimmed Data QC Part II: `MultiQC`   
In the terminal of the same Rstudio project used in Part D, run:    

```
eval "$(/opt/anaconda/anaconda3/bin/conda shell.bash hook)"
conda activate
```

Then, navigate into the directory where the FASTQC output lives, in this case: `/home/shared/8TB_HDD_02/graceac9/analyses/pycno2022`, then run:   
```
multiqc .
```

The report will generate in seconds to minutes.       
To view the report, transfer the .html report to Gannet or Owl, then you can view the .html report on your own browser.

I moved the untrimmed MultiQC report to Owl: navigate in terminal to directory where the .html report lives, then rsync file to where I want it on owl.

```
graceac9@raven:~/analyses/pycno2022$ rsync --archive --progress --verbose multiqc_report.html grace@owl.fish.washington.edu:/volume1/web/scaphapoda/grace/pycno_2021/multiqc
grace@owl.fish.washington.edu's password:
sending incremental file list
multiqc_report.html
      1,853,263 100%   72.34MB/s    0:00:00 (xfr#1, to-chk=0/1)

sent 1,853,829 bytes  received 34 bytes  195,143.47 bytes/sec
total size is 1,853,263  speedup is 1.00
graceac9@raven:~/analyses/pycno2022$
```

Untrimmed RNAseq data `MultiQC report`: [owl/scaphapoda/grace/pycno_2022/multiqc/multiqc_report.htm](http://owl.fish.washington.edu/scaphapoda/grace/pycno_2022/multiqc/multiqc_report.html).

## 3. Trim RNAseq data: Run `fastp` on Mox, then `MultiQC`
### A. `rsync` RNAseq data (fastq.gz) from `nightingales` to `/gscratch/srlab/graceac9/data/pycno/RNAseq/summer2022`.

#### i. Navigate into `/nightingales/P_helianthoides` in command line on Owl.
#### ii. Copy in code:     
```
rsync —archive —progress —verbose PSC-0*.fastq.gz graceac9@mox.hyak.uw.edu:/gscratch/srlab/graceac9/data/pycno/RNAseq/summer2022
```

#### iii. You'll be prompted for Mox password and 2FA, then you'll be good to go.

Takes ~2 hours for 32 libraries.

### B. Add .sh script to Mox.
`/gscratch/srlab/graceac9/jobs/20231206_pycno2022_fastp.sh`

Also, saved script to respository: [project-pycno-sizeclass-2022/code/02-20231206_pycno2022_fastp.sh](https://raw.githubusercontent.com/grace-ac/project-pycno-sizeclass-2022/main/code/02-20231206_pycno2022_fastp.sh)

### C. Run the job.
Navigate into `/gscratch/srlab/graceac9/jobs`, then run:
```
sbatch 20231206_pycno2022_fastp
```

and it will put output into `/gscratch/srlab/graceac9/analyses/pycno/20231206_PSC2022_trimming`


Check job status by running:
```
squeue | grep "srlab"
```

20231206 - started 14:26. Changed wall time to 5 days, because my 10 day wall time conflicted with the monthly Tuesday Mox maintenance. So I'll have to pause the job before maintenance and re-start after.   

20231206 - ended 16:06.

The output will be that each library has 4 files:
PSC-0###_R1_001.fastq.gz.fastp-trim.20231206.fq.gz
PSC-0###_R1_001.fastq.gz.fastp-trim.20231206.report.html
PSC-0###_R1_001.fastq.gz.fastp-trim.20231206.report.json
PSC-0###_R2_001.fastq.gz.fastp-trim.20231206.fq.gz

Where ## is the library number. Since they're paired end, the reports (report.html and report.json) contain info for both sets of reads (note from Sam White).

### D. Move the multiqc report to OWL
In same area where other one lives, but add the date and that it's for trimmed data.

Navigate in terminal to directory on Mox where the multiqc report lives

Then:     
```
[graceac9@mox2 20231206_PSC2022_trimming]$ rsync --archive --progress --verbose multiqc_report.html grace@owl.fish.washington.edu:/volume1/web/scaphapoda/grace/pycno_2022/multiqc/trimmed
grace@owl.fish.washington.edu's password:
sending incremental file list
multiqc_report.html
      1,590,859 100%  185.74MB/s    0:00:00 (xfr#1, to-chk=0/1)

sent 1,591,363 bytes  received 34 bytes  28,673.82 bytes/sec
total size is 1,590,859  speedup is 1.00
```

Trimmed MultiQC report: [owl.fish.washington.edu/scaphapoda/grace/pycno_2022/multiqc/trimmed/multiqc_report.html](http://owl.fish.washington.edu/scaphapoda/grace/pycno_2022/multiqc/trimmed/multiqc_report.html)

# 3. Align trimmed RNAseq data to published _Pycnopodia helianthoides_ genome using `HISAT2`

## 1. Move genome and annotation to Mox:
I did this as part of my Summer 2021 work:

### A. Download .zip genome onto laptop
Navigate to downloads in command line and run:  

```
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_032158295.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCA_032158295.1.zip" -H "Accept: application/zip
```

Code got from NCBI genbank page: [https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032158295.1/](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032158295.1/)

### B. `rsync` .zip file to Mox directory:
Be in Downloads folder in command line on laptop, then run:
```
rsync --archive --progress --verbose /path/to/file username@mox_IP:/path/to/mox/directory
```

The .zip file is now on Mox:
`/gscratch/srlab/graceac9/data/pycno/genome/GCA_032158295.1.zip`

### C. Get genome annotation file
interproscan.gff3 from [https://datadryad.org/stash/dataset/doi:10.5061/dryad.51c59zwfd](https://datadryad.org/stash/dataset/doi:10.5061/dryad.51c59zwfd)

In the same Mox directory where the .zip genome lives (`/gscratch/srlab/graceac9/data/pycno/genome`), run:
```
wget https://datadryad.org/stash/downloads/file_stream/2634383
```

Then in the `/gscratch/srlab/graceac9/data/pycno/genome` directory, you'll have:
```
[graceac9@mox2 genome]$ pwd
/gscratch/srlab/graceac9/data/pycno/genome
```

### D. Genome is too large to be unzipped on regular Mox, so unzipping will happen in `/gscratch/scrubbed`

#### i. `rsync` .zip genome from laptop downloads into scrubbed directory:
```
rsync --archive --progress --verbose GCA_032158295.1.zip graceac9@mox.hyak.uw.edu:/gscratch/scrubbed/graceac9/ncbi_dataset/data
```

It's now on Mox:
```
[graceac9@mox2 data]$ ls
GCA_032158295.1.zip
```
#### ii. unzip the file
```
unzip GCA_032158295.1.zip
```
Things got a bit messy, but it's now unzipped in a series of new directories:

`/gscratch/scrubbed/graceac9/ncbi_dataset/data/ncbi_dataset/data/GCA_032158295.1`

And the genome itself is: `/gscratch/scrubbed/graceac9/ncbi_dataset/data/ncbi_dataset/data/GCA_032158295.1/GCA_032158295.1_ASM3215829v1_genomic.fna`

The fasta cannot be moved out of `/gscratch/scrubbed/` because it is far too big!

## 2. Run `HISAT2` with genome and trimmed RNAseq reads
Make a directory for 2022 `HISAT2` output:     
`/gscratch/srlab/graceac9/analyses/20231207-hisat2-2022data`

### A. Create .sh for `HISAT2`
in `/gscratch/scrubbed/graceac9/jobs`, create a .sh script.            20231207_hisat2_2022pycno_align.sh

Be sure to save a copy to GitHub repo, because things get deleted from `/gscratch/scrubbed` every 21 days!       
Script: [project-pycno-sizeclass-2022/code/02-20231207_hisat2_2022pycno_align.sh](https://raw.githubusercontent.com/grace-ac/project-pycno-sizeclass-2022/main/code/02-20231207_hisat2_2022pycno_align.sh)

Additionally --> back up all of Mox to gannet. Log into gannet.fish.washington.edu on command line, navigate to `/volume2/web/gcrandall/bu-mox` and `rsync` contents:  
```
rsync -avz --progress graceac9@mox.hyak.uw.edu:/path/to/what/you/want/to/move /volume2/web/gcrandall/bu-mox
```

### B. Run `HISAT2` script:
Be in directory that script lives in: `/gscratch/scrubbed/graceac9/jobs`
```
sbatch 20231207_hisat2_2022pycno_align.sh
```

Began 20231207 16:01    
Ended 20231207 18:59

Output files are apparently here: `/gscratch/scrubbed/graceac9/ncbi_dataset/data/ncbi_dataset/data/GCA_032158295.1`, but I thought they were supposed to be .bam files... and I'm not sure why there are 8.

```
[graceac9@mox1 GCA_032158295.1]$ pwd
/gscratch/scrubbed/graceac9/ncbi_dataset/data/ncbi_dataset/data/GCA_032158295.1
[graceac9@mox1 GCA_032158295.1]$ ls
GCA_032158295.1_ASM3215829v1_genomic.fna  Phelianthoides_ref.4.ht2  Phelianthoides_ref.8.ht2
Phelianthoides_ref.1.ht2		  Phelianthoides_ref.5.ht2  sequence_report.jsonl
Phelianthoides_ref.2.ht2		  Phelianthoides_ref.6.ht2
Phelianthoides_ref.3.ht2		  Phelianthoides_ref.7.ht2
```
