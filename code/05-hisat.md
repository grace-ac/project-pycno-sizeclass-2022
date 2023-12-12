05-Hisat
================
2023-12-12

# Trimmed Reads

On Gannet
<https://gannet.fish.washington.edu/gcrandall/bu-mox/20231206_PSC2022_trimming/>

``` bash
wget -r \
--no-directories --no-parent \
-P ../data \
-A "*.fq.gz" https://gannet.fish.washington.edu/gcrandall/bu-mox/20231206_PSC2022_trimming/
```

# HiSat Indexing of Genome

# genome

``` bash
cd ../data
/home/shared/datasets download genome accession GCA_032158295.1 --include gff3,rna,cds,protein,genome,seq-report
```

\#genome prep

Just fasta… 01

``` bash
/home/shared/hisat2-2.2.1/hisat2-build \
../data/GCA_032158295.1_ASM3215829v1_genomic.fna \
../analyses/05-hisat/GCA_032158295.1_01 \
-p 40 \
> ../analyses/05-hisat/GCA_032158295.1_01.txt
```

# genome feature files

<https://datadryad.org/stash/dataset/doi:10.5061/dryad.51c59zwfd>

``` bash
head ../data/aug*
ls ../data/aug*
```

# extract exons and splice sites

``` bash
/home/shared/hisat2-2.2.1/hisat2_extract_exons.py \
../data/augustus.hints.gtf \
> ../analyses/05-hisat/exon.tab
```

``` bash
head ../analyses/05-hisat/exon.tab
```

    ## pycn_heli.0001   46404   46597   -
    ## pycn_heli.0001   47866   47893   -
    ## pycn_heli.0001   53468   53545   -
    ## pycn_heli.0001   60935   60988   -
    ## pycn_heli.0001   63915   64003   -
    ## pycn_heli.0001   65837   65926   -
    ## pycn_heli.0001   67817   67907   -
    ## pycn_heli.0001   70745   70890   +
    ## pycn_heli.0001   73706   73787   +
    ## pycn_heli.0001   75928   76088   +

``` bash
/home/shared/hisat2-2.2.1/hisat2_extract_splice_sites.py \
../data/augustus.hints.gtf \
> ../analyses/05-hisat/splice_sites.tab
```

``` bash
head ../analyses/05-hisat/splice_sites.tab
```

    ## pycn_heli.0001   46597   47866   -
    ## pycn_heli.0001   47893   53468   -
    ## pycn_heli.0001   60988   63915   -
    ## pycn_heli.0001   64003   65837   -
    ## pycn_heli.0001   65926   67817   -
    ## pycn_heli.0001   70890   73706   +
    ## pycn_heli.0001   73787   75928   +
    ## pycn_heli.0001   76088   76784   +
    ## pycn_heli.0001   76922   78411   +
    ## pycn_heli.0001   78512   80347   +

# indexing with splice sites

2)  

``` bash
/home/shared/hisat2-2.2.1/hisat2-build \
../data/GCA_032158295.1_ASM3215829v1_genomic.fna \
../analyses/05-hisat/GCA_032158295.1_02 \
--exon ../analyses/05-hisat/exon.tab \
--ss ../analyses/05-hisat/splice_sites.tab \
-p 40 \
../data/augustus.hints.gtf \
2> ../analyses/05-hisat/GCA_032158295.1_02.txt
```

# Alignment

    Hisat2 alignments
    "${programs_array[hisat2]}" \
    -x "${genome_index_name}" \
    -1 "${fastq_list_R1}" \
    -2 "${fastq_list_R2}" \
    -S "${sample_name}".sam \
    2> "${sample_name}"-hisat2_stats.txt

test run

``` bash
/home/shared/hisat2-2.2.1/hisat2 \
-x ../analyses/05-hisat/GCA_032158295.1_02 \
-p 48 \
-1 ../data/PSC-0235_R1_001.fastq.gz.fastp-trim.20231206.fq.gz \
-2 ../data/PSC-0235_R2_001.fastq.gz.fastp-trim.20231206.fq.gz \
-S ../analyses/05-hisat/PSC-0235.sam \
> ../analyses/05-hisat/PSC-0235_hisat2_stats.txt
```

``` bash

# Loop through all  files
for file in ../data/*_R1*fq.gz; do
    # Remove the  part to get the base name
    base=$(basename "$file" _R1_001.fastq.gz.fastp-trim.20231206.fq.gz)
    # Construct the names of the pair of files
    file1=${base}_R1_001.fastq.gz.fastp-trim.20231206.fq.gz
    file2=${base}_R2_001.fastq.gz.fastp-trim.20231206.fq.gz
    # Run the hisat2 command
    /home/shared/hisat2-2.2.1/hisat2 \
    -x ../analyses/05-hisat/GCA_032158295.1_02 \
    -p 20 \
    -1 ../data/$file1 \
    -2 ../data/$file2 \
    -S ../analyses/05-hisat/${base}.sam \
    2> ../analyses/05-hisat/${base}_hisat2_stats.txt
done
```

``` bash
for samfile in ../analyses/05-hisat/*.sam; do
  bamfile="${samfile%.sam}.bam"
  sorted_bamfile="${samfile%.sam}.sorted.bam"
  
  # Convert SAM to BAM
  /home/shared/samtools-1.12/samtools view -bS -@ 20 "$samfile" > "$bamfile"
  
  # Sort BAM
  /home/shared/samtools-1.12/samtools sort -@ 20 "$bamfile" -o "$sorted_bamfile"
  
  # Index sorted BAM
  /home/shared/samtools-1.12/samtools index -@ 20 "$sorted_bamfile"
done
```
