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

## 1. Data Management
Received data from Azenta (FKA Genewiz).

Downloaded to genefish using: [cyberduck](https://cyberduck.io/download/).     
Azenta provided log-in credentials and information of how to access their servers. Azenta pdf of different ways to download the data from their servers: [here](https://f.hubspotusercontent00.net/hubfs/3478602/Sell%20Sheet%20Collateral%20Library/NGS/NGS%20User%20Guides/NGS_sFTP-Data-Download-Guide_Option%201_Nov03_2020.pdf).

All downloaded files were put into a folder in Documents on genefish called: "P_helianthoides_RNAseq", with a subfolder for Summer 2022 data is in folder: "30-833328413".

### Move data to OWL     
#### A. `rsync` data to `owl/nightingales/P_helianthoides`
On genefish, in commandline in the directory where the data lives:
```
rsync --archive --progress --verbose PSC*.fastq.gz <owl username>@128.95.149.83:/volume1/web/nightingales/P_helianthoides/
```
Replace <owl_username_> with whatever username you use to login to owl (even replace the < and the >).

#### B. `rsync` checksums to OWL on commandline:     
##### B.1. `ssh` into OWL:     
```
ssh username@owl.fish.washington.edu
```
Will be prompted for password.    

##### B.2. Navigate to folder you want to work in, in my case: `owl/nightingales/P_helianthoides`

##### B.3. Put below code in and hit ENTER:     
```
for fastq in PSC*.fastq.gz
do
  md5sum "${fastq}" >> checksums.md5
  echo "Generated checksum for ${file}."
  echo ""
done
```

All RNAseq data and checksums are now in [owl/nightingales/P_Helianthoides](https://owl.fish.washington.edu/nightingales/P_helianthoides/).

## 2. QC and Trimming
### 1. Untrimmed Data QC Part I: `FASTQC`  
#### A. Get FastQC if you want to run on your laptop:    
[https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)     

#### B. Get .fastq.gz files from OWL onto RAVEN   
##### i. `ssh` into Raven using credentials in command line.    
##### ii. Make a directory for all PSC.fastq.gz files (made one called `pycnornaseq2022` in my `/home/shared/8TB_HDD_02/graceac9/pycnornaseq2022`)  

#### C. Then move files from OWL to pycnornaseq directory in Raven:

##### a. Have Husky OnNet App (BIG-IP Edge Client in Applications folder after downloaded)
##### b. Log in with UW credentials
##### c. Put RStudio IP into browser: http://172.25.149.12:8787
##### d. Log in using Raven Credentials
##### e. cd into pycnornaseq/ and run:
```
wget -r --no-directories --no-parent  -A "PSC-0*" https://owl.fish.washington.edu/nightingales/P_helianthoides
```

Then the directory will look like this:     
<img width="573" alt="PSC_files_on_RAVEN" src="https://github.com/grace-ac/project-pycno-sizeclass-2022/blob/main/protocols/images/screenshot-1.png">

#### D. Get into Rstudio on Raven to run FASTQC:  
Follow steps a-d in Step C listed above.

Then, follow the code outline in this script: [project-pycno-sizeclass-2022/code/01-FastQC_pre-trim.Rmd](https://github.com/grace-ac/project-pycno-sizeclass-2022/blob/main/code/01-FastQC_pre-trim.Rmd)

The FASTQC files are saved on Raven: `/home/shared/8TB_HDD_02/graceac9/analyses/pycno2022`.

### 2. Untrimmed Data QC Part II: `MultiQC`   
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
