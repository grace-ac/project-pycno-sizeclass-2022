05-Hisat
================
2023-12-12

Code for RNAseq of 2022 SeaStar experiment

# Sample info

``` r
mdata <- read.delim("../data/summer2022_samples_sequenced.csv", sep = ",")
head(mdata)
```

    ##   sample_date star_ID star_size treatment_grp experiment_day disease_sign
    ## 1    20220621       6     large       exposed              0      healthy
    ## 2    20220621      31     large       control              0      healthy
    ## 3    20220621      27     large       control              0      healthy
    ## 4    20220621      13     large       exposed              0      healthy
    ## 5    20220621      10     large       exposed              0      healthy
    ## 6    20220621       4     large       control              0      healthy
    ##   coelom_ID extraction_date                         kit          elution_liquid
    ## 1  PSC 0085      2023-02-27 Zymo-quickDNA/RNA-microprep RT DNAse-RNAse free H20
    ## 2  PSC 0087      2023-03-09 Zymo-quickDNA/RNA-microprep RT DNAse-RNAse free H20
    ## 3  PSC 0090      2023-03-09 Zymo-quickDNA/RNA-microprep RT DNAse-RNAse free H20
    ## 4  PSC 0093      2023-03-09 Zymo-quickDNA/RNA-microprep RT DNAse-RNAse free H20
    ## 5  PSC 0095      2023-03-09 Zymo-quickDNA/RNA-microprep RT DNAse-RNAse free H20
    ## 6  PSC 0096      2023-03-29 Zymo-quickDNA/RNA-microprep RT DNAse-RNAse free H20
    ##   elution_vol_ul qubit_run_date ul_run_on_qubit qubit_std_1 qubit_std_2
    ## 1             15     2023-03-02               1       82.15     1717.37
    ## 2             15     2023-03-09               1       82.14     1565.32
    ## 3             15     2023-03-09               1       82.14     1565.32
    ## 4             15     2023-03-09               1       82.14     1565.32
    ## 5             15     2023-03-09               1       82.14     1565.32
    ## 6             15     2023-03-29               1          NA          NA
    ##   RNA_readout units Nanodrop_run_date Nanodrop_ul.run X260.280 X260.230
    ## 1        9.92 ng/ul                NA              NA       NA       NA
    ## 2        20.4 ng/ul                NA              NA       NA       NA
    ## 3        9.62 ng/ul                NA              NA       NA       NA
    ## 4        84.8 ng/ul                NA              NA       NA       NA
    ## 5        70.2 ng/ul                NA              NA       NA       NA
    ## 6        <NA> ng/ul                NA              NA       NA       NA
    ##   RNA_ng.ul qubit_run_date_2 ul_run_on_qubit_2 qubit_std1_2 qubit_std2_2 RNA_2
    ## 1        NA                                 NA           NA           NA    NA
    ## 2        NA                                 NA           NA           NA    NA
    ## 3        NA                                 NA           NA           NA    NA
    ## 4        NA                                 NA           NA           NA    NA
    ## 5        NA                                 NA           NA           NA    NA
    ## 6        NA       2023-03-30                 1        96.77       1806.8  78.6
    ##       X units.1 final_vol_ul total.RNA_ng morethan_250ng_RNA_yes_no
    ## 1  9.92                   14       138.88                     FALSE
    ## 2 20.40                   14       285.60                      TRUE
    ## 3  9.62                   14       134.68                     FALSE
    ## 4 84.80                   14      1187.20                      TRUE
    ## 5 70.20                   14       982.80                      TRUE
    ## 6 78.60   ng/ul           13      1021.80                      TRUE
    ##   morethan_100ng_total_RNA_yes_no conc_original_qubt_10ngul_t_o_f
    ## 1                            TRUE                           FALSE
    ## 2                            TRUE                            TRUE
    ## 3                            TRUE                           FALSE
    ## 4                            TRUE                            TRUE
    ## 5                            TRUE                            TRUE
    ## 6                            TRUE                            TRUE
    ##   send_for_seq_yes.or.no
    ## 1                     NA
    ## 2                     NA
    ## 3                     NA
    ## 4                     NA
    ## 5                     NA
    ## 6                     NA
    ##                                                                    notes
    ## 1                                                                       
    ## 2                     treated cells like blood and did modified protocol
    ## 3                     treated cells like blood and did modified protocol
    ## 4                     treated cells like blood and did modified protocol
    ## 5                     treated cells like blood and did modified protocol
    ## 6 removed 1ul for qubit, but messed up math, so will have to re-do qubit
    ##   RNA_.080C_yes_no notes_on_star_at_sampling_time
    ## 1               NA                               
    ## 2               NA                               
    ## 3               NA                               
    ## 4               NA                               
    ## 5               NA                               
    ## 6               NA               arms twisted - 7

``` r
#simplified metadata
smdata <- select(mdata, coelom_ID, star_ID, star_size, treatment_grp, experiment_day, disease_sign)

smdata
```

    ##    coelom_ID star_ID star_size treatment_grp experiment_day
    ## 1   PSC 0085       6     large       exposed              0
    ## 2   PSC 0087      31     large       control              0
    ## 3   PSC 0090      27     large       control              0
    ## 4   PSC 0093      13     large       exposed              0
    ## 5   PSC 0095      10     large       exposed              0
    ## 6   PSC 0096       4     large       control              0
    ## 7   PSC 0099      28     large       exposed              0
    ## 8   PSC 0105      23     large       control              0
    ## 9   PSC 0118      60     small       control              0
    ## 10  PSC 0119      63     small       control              0
    ## 11  PSC 0132      55     small       control              0
    ## 12  PSC 0141      44     small       exposed              0
    ## 13  PSC 0149      28     large       exposed              8
    ## 14  PSC 0150       6     large       exposed              8
    ## 15  PSC 0156      28     large       exposed              9
    ## 16  PSC 0174      30     large       exposed             10
    ## 17  PSC 0177      34     small       control             11
    ## 18  PSC 0186      15     large       control             11
    ## 19  PSC 0187      57     small       exposed             11
    ## 20  PSC 0188      38     small       exposed             11
    ## 21  PSC 0190      10     large       exposed             11
    ## 22  PSC 0198      38     small       exposed             12
    ## 23  PSC 0202       4     large       control             12
    ## 24  PSC 0203      31     large       control             12
    ## 25  PSC 0206      30     large       exposed             12
    ## 26  PSC 0209       9     large       control             13
    ## 27  PSC 0217       6     large       exposed             13
    ## 28  PSC 0219      63     small       control             14
    ## 29  PSC 0228      61     small       exposed             14
    ## 30  PSC 0230      65     small       control             15
    ## 31  PSC 0231      13     large       exposed             15
    ## 32  PSC 0235      61     small       exposed             15
    ##                           disease_sign
    ## 1                              healthy
    ## 2                              healthy
    ## 3                              healthy
    ## 4                              healthy
    ## 5                              healthy
    ## 6                              healthy
    ## 7                              healthy
    ## 8                              healthy
    ## 9                              healthy
    ## 10                             healthy
    ## 11                             healthy
    ## 12                             healthy
    ## 13    FIRST arm dropped, 1 arm twisted
    ## 14   FIRST arm dropped, 4 arms twisted
    ## 15                  SECOND arm dropped
    ## 16   FIRST arm dropped, 6 arms twisted
    ## 17                             healthy
    ## 18                             healthy
    ## 19 1 arm dropped, twisting, stretching
    ## 20           1 arm dropped, stretching
    ## 21                   FIRST arm dropped
    ## 22       1 arm dropped, 4 arms twisted
    ## 23                             healthy
    ## 24                             healthy
    ## 25                  SECOND arm dropped
    ## 26                             healthy
    ## 27                      9 arms dropped
    ## 28                             healthy
    ## 29                       1 arm dropped
    ## 30                             healthy
    ## 31   FIRST arm dropped, 5 arms twisted
    ## 32   lost last 3 arms; just disc; dead

``` r
datatable(smdata, options = list(scrollX = TRUE, scrollY = "400px", scrollCollapse = TRUE, paging = FALSE))
```

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

then moved around

\#genome prep

Just fasta… 01 - will not initially use as does not have splice site
info. Will keep if ever want to do comparison.

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

    ## ==> ../data/augustus.hints.aa <==
    ## >g3450.t1
    ## MNYRTDDEIYEDEVDEETKLHHRIDERGVNGTDSRDEPTRVPPAKQLKIQAPVLASRLQL
    ## EAVRRPHPPPLHPPQIHLYLEPHLEALQQRCDSYKGIVG*
    ## >g3451.t1
    ## MKINKEKSKVMHLSRCKTPTDEHLPHYTLQETQSKLFCRLKGLRDYNFRSFILAFLQYNP
    ## SDLISVSDAFLHPFLARTSCSKVPCPLLPKMKIRRRQLGITHHRKTT*
    ## >g3452.t1
    ## MSWTLLWDSLMTVALTMSVLPVCAFVILVLYIAYEHRRFSHIPGPPRKEMENYGTVFSLF
    ## MFQKPVVVCLEPSFVKKLLYSTTHIKTPAEVKYFWRVFGQRYLHHGLLTECDVPKNLKRR
    ## ALFEPAFHRKYLKTLMCTFNESVSKMIERLTLKADGKTEVFMLDELNKLTLDVIAKTAFG
    ## 
    ## ==> ../data/augustus.hints.codingseq <==
    ## >g3450.t1
    ## ATGAACTACCGGACTGATGATGAAATATATGAAGATGAGGTGGACGAGGAGACAAAACTG
    ## CATCACAGAATTGATGAGCGTGGAGTGAACGGAACTGATTCTCGAGATGAACCAACAAGA
    ## GTTCCTCCAGCCAAGCAGCTGAAGATTCAAGCACCAGTATTGGCGTCTCGGCTCCAGCTG
    ## GAAGCAGTGAGGCGGCCCCACCCCCCTCCTCTTCATCCTCCTCAGATCCACCTGTACCTA
    ## GAACCACACCTAGAGGCATTACAACAGCGTTGTGATAGCTACAAAGGCATTGTTGGATAA
    ## >g3451.t1
    ## ATGAAGATCAACAAAGAAAAATCCAAAGTGATGCATTTGTCTCGCTGCAAAACACCCACG
    ## GACGAGCATCTTCCACATTACACACTGCAAGAGACACAGAGCAAACTGTTCTGCAGACTG
    ## AAGGGGCTGAGAGATTACAACTTTCGCAGTTTCATTCTGGCCTTTCTGCAGTACAATCCA
    ## 
    ## ==> ../data/augustus.hints.gtf <==
    ## pycn_heli.0392   AUGUSTUS    gene    1   4781    0.95    +   .   g1
    ## pycn_heli.0392   AUGUSTUS    transcript  1   4781    0.95    +   .   g1.t1
    ## pycn_heli.0392   AUGUSTUS    intron  1   4506    0.95    +   .   transcript_id "g1.t1"; gene_id "g1";
    ## pycn_heli.0392   AUGUSTUS    CDS 4507    4781    0.95    +   2   transcript_id "g1.t1"; gene_id "g1";
    ## pycn_heli.0392   AUGUSTUS    exon    4507    4781    .   +   .   transcript_id "g1.t1"; gene_id "g1";
    ## pycn_heli.0392   AUGUSTUS    stop_codon  4779    4781    .   +   0   transcript_id "g1.t1"; gene_id "g1";
    ## pycn_heli.0392   AUGUSTUS    gene    16379   16612   0.78    -   .   g2
    ## pycn_heli.0392   AUGUSTUS    transcript  16379   16612   0.78    -   .   g2.t1
    ## pycn_heli.0392   AUGUSTUS    stop_codon  16379   16381   .   -   0   transcript_id "g2.t1"; gene_id "g2";
    ## pycn_heli.0392   AUGUSTUS    CDS 16379   16612   0.78    -   0   transcript_id "g2.t1"; gene_id "g2";
    ## ../data/augustus.hints.aa
    ## ../data/augustus.hints.codingseq
    ## ../data/augustus.hints.gtf

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

# Indexing with splice sites

Genome indexing- `02` suffix; has exon and splice site tbas
incorporated.

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

# Converting sam to sorted bams

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

# Stringtie

``` bash
cd ../data
curl -O https://gannet.fish.washington.edu/seashell/snaps/interproscan.gff3
```

``` bash
head ../data/interproscan.gff3 
```

    ## ##gff-version 3
    ## ##feature-ontology http://song.cvs.sourceforge.net/viewvc/song/ontology/sofa.obo?revision=1.269
    ## ##interproscan-version 5.61-93.0
    ## ##sequence-region g3531.t1 1 635
    ## g3531.t1 .   polypeptide 1   635 .   +   .   ID=g3531.t1;md5=e585d160537c159881dfb8cf27e10150
    ## g3531.t1 CDD protein_match   288 537 2.58694E-146    +   .   date=19-04-2023;Target=g3531.t1 288 537;ID=match$1_288_537;signature_desc=ABC1_ADCK3;Name=cd13970;status=T;Dbxref="InterPro:IPR034646"
    ## g3531.t1 PANTHER protein_match   6   633 8.3E-243    +   .   date=19-04-2023;Target=g3531.t1 6 633;ID=match$2_6_633;Name=PTHR43851;status=T
    ## g3531.t1 MobiDBLite  protein_match   151 176 .   +   .   date=19-04-2023;Target=g3531.t1 151 176;ID=match$3_151_176;signature_desc=consensus disorder prediction;Name=mobidb-lite;status=T
    ## g3531.t1 MobiDBLite  protein_match   74  200 .   +   .   date=19-04-2023;Target=g3531.t1 74 200;ID=match$4_74_200;signature_desc=consensus disorder prediction;Name=mobidb-lite;status=T
    ## g3531.t1 MobiDBLite  protein_match   83  149 .   +   .   date=19-04-2023;Target=g3531.t1 83 149;ID=match$5_83_149;signature_desc=consensus disorder prediction;Name=mobidb-lite;status=T

``` bash
find ../analyses/05-hisat/*sorted.bam \
| xargs basename -s .sorted.bam | xargs -I{} \
/home/shared/stringtie-2.2.1.Linux_x86_64/stringtie \
-p 20 \
-o ../analyses/05-hisat/{}.gtf \
../analyses/05-hisat/{}.sorted.bam
```

``` bash
head ../analyses/05-hisat/*.gtf | head -25
```

    ## ==> ../analyses/05-hisat/PSC-0085.gtf <==
    ## # /home/shared/stringtie-2.2.1.Linux_x86_64/stringtie -p 20 -o ../analyses/05-hisat/PSC-0085.gtf ../analyses/05-hisat/PSC-0085.sorted.bam
    ## # StringTie version 2.2.1
    ## CM063243.1   StringTie   transcript  16862   21750   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "7.316682"; FPKM "2.272532"; TPM "6.570474";
    ## CM063243.1   StringTie   exon    16862   17039   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "8.224985";
    ## CM063243.1   StringTie   exon    18567   18598   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "5.343992";
    ## CM063243.1   StringTie   exon    21627   21750   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "3"; cov "6.521910";
    ## CM063243.1   StringTie   transcript  16862   19857   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.2"; cov "7.869413"; FPKM "2.444208"; TPM "7.066833";
    ## CM063243.1   StringTie   exon    16862   17039   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.2"; exon_number "1"; cov "15.679487";
    ## CM063243.1   StringTie   exon    18567   19857   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.2"; exon_number "2"; cov "6.792579";
    ## CM063243.1   StringTie   transcript  39611   40853   1000    -   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; cov "10.447793"; FPKM "3.245042"; TPM "9.382252";
    ## 
    ## ==> ../analyses/05-hisat/PSC-0087.gtf <==
    ## # /home/shared/stringtie-2.2.1.Linux_x86_64/stringtie -p 20 -o ../analyses/05-hisat/PSC-0087.gtf ../analyses/05-hisat/PSC-0087.sorted.bam
    ## # StringTie version 2.2.1
    ## CM063243.1   StringTie   transcript  16861   21728   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "13.773430"; FPKM "4.588499"; TPM "15.322474";
    ## CM063243.1   StringTie   exon    16861   17039   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "14.391528";
    ## CM063243.1   StringTie   exon    18567   18598   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "19.468750";
    ## CM063243.1   StringTie   exon    21627   21728   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "3"; cov "10.901960";
    ## CM063243.1   StringTie   transcript  16861   36000   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.2"; cov "3.747169"; FPKM "1.248337"; TPM "4.168599";
    ## CM063243.1   StringTie   exon    16861   17039   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.2"; exon_number "1"; cov "3.333803";
    ## CM063243.1   StringTie   exon    27014   27112   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.2"; exon_number "2"; cov "5.232324";
    ## CM063243.1   StringTie   exon    35926   36000   1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.2"; exon_number "3"; cov "2.773334";
    ## 
    ## ==> ../analyses/05-hisat/PSC-0090.gtf <==
