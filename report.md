# AAT Variant Calling Report

## Animal Sampling
Samples were collected in Shimba Hills National Park, Kwale county as shown in the map in Figure 1. 
The selection of the study area was determined by the high prevalence rate of trypanosomes in and around Shimba Hills National Park seen from previous trypanosome 
surveys and studies carried out.

![shimba_hills_map](https://user-images.githubusercontent.com/68329457/125247526-b669fc80-e2fb-11eb-8067-f735279d5e7d.jpg)

Figure 1: Geographical map of the traps uded to collect tsetse flies (Trap 1-20) and the two villages (Msulwa and Mlafyeni) targeted for blood collection in
livestock in Shimba Hills National Park, Kwale County, Kenya.

The traps were located across a transect from Msulwa and Mlafyeni villages with the aim of tracking the transmission of AAT from livestock outside the park and wildlife inside the park and vice versa. 
This will allow us to get an insight of the role wildlife have in the transmission of AAT in Shimba Hills National Park.

In this case, tsetse flies and livestock blood samples were collected and processed to explore the nature of the data acquired with the first round of the field work.
This is important as it enables us to optimize protocols used in the laboratory processes as well as bioinformatics analysis.

164 livestock blood samples were collected in EDTA coated vacutainers from both villages while ~1800 tsetse flies were collected from traps set in the park and preserved in falcon tubes containing 100% ethanol.
These were transported back to Nairobi for laboratory analysis.

### Laboratory Analysis

Tsetse flies and blood samples were extracted using a manual extraction method as well as a kit (Bioline) to assess the DNA quality and yield from both methods. The blood samples had a much lower DNA yield as compared to tsetse flies when using the kit, while tsetse flies had a lower DNA quality when using the manual extraction method. 
This resulted in using the Bioline kit in tsetse samples and a manual extraction method for blood samples for DNA extraction. 
All the blood samples and 327 tsetse flies were processed. The manual method used in DNA extraction from blood is as follows:

i. Add 600µl of Low Salt Buffer (10mM Tris-HCl pH 7.6, 10mM KCl, 10mM MgCl2, 2mM EDTA) and 35µl of triton X to a micro centrifuge tube containing 200µl of blood sample.

ii. Incubate at 56ºC for 5 minutes

iii. Centrifuge at 10000 rpm for 3 minutes and discard the supernatant

iv. Repeat step 3 three times with decreasing amount of 1X triton X until the RBC lysis is complete and white pellet is obtained

v. Add 300µl of High Salt Buffer (10mM Tris-HCl pH 7.6, 10mM KCl, 10mM MgCl2, 2mM EDTA, 0.4M NaCl) and 40µl of 10% SDS to the cell pellet and mix thoroughly
 
vi. Incubate at 56ºC for 5 minutes

vii. Add 100µl of 6M NaCl and vortex thoroughly.

viii. Centrifuge at 10000rpm for 5 minutes

ix. Transfer the supernatant into a new tube containing 300µl of isopropanol

x. Invert the tubes gently for 3 minutes

xi. Centrifuge at 10000 rpm for 10 minutes to pellet DNA

xii. Discard supernatant and add 500µl of 70% ethanol and mix slowly to remove excess salts

xiii. Centrifuge the tubes at 10000rpm for 5 minutes and discard the supernatant

xiv. Invert the tubes on a paper towel to dry for 1 hour

xv. Elute with 50µl of Tris-EDTA (TE) buffer


The extracted DNA samples were screened for Trypanosoma species using ITS1 primers. 
The table below (Table 1) summarizes the results of the PCR runs in both blood and tsetse flies.


| Trypanosoma Species | Blood samples | Tsetse flies |
|---------------------|---------------|--------------|
| Trypanosoma congolense | 10 | 3 |
| Trypanosoma vivax | 3 | 1 |
| Trypanosoma brucei | 3 | 2 |
| Mixed infections | 1 | 1 |
| Total infections | 17 | 7 |

Table 1 : Summary of the trypanosome infections in blood and tsetse flies from samples collected in and outside Shimba Hills Park

For Whole Genome Sequencing (WGS), the samples are required to have a DNA concentration of 500-1000ng; 260/280nm ratio of 1.8-2.0 and 260/230nm ratio of 2.0-2.2.
Four samples were selected for WGS and have been summarized in table 2

| Sample No.| Sample type | Location | Extraction method | DNA Concentration(ng/uL) | A260/A280 | A260/A230 | Diagnosis | Blood fed |
|--------|-------------|-------------------|-----------------|--------------------|----------|----------|---------|---------|
| T8 | Tsetse | Trap 12 | Bioline Kit | 104.9 | 2.10 | 2.26 | T.V | Yes |
| T70 | Tsetse | Trap 18 | Bioline Kit | 129.0 | 2.00 | 1.59 | T.C | Yes |
| B87 | Blood | Msulwa | Manual | 557.7 | 1.69 | 1.11 | T.C | - |
| B151 | Blood | Msulwa | Manual | 351.8 | 1.79 | 1.59 | T.V | - |

Table 2 : Summary of the samples selected for sequencing ; with two blood samples and two tsetse flies each with Trypanosoma congolense (T.c) and Trypanosoma
vivax (T.v) infections respectively.

The samples were sequenced on the Next Seq Illumina sequencing platform resulting in a billion pair-ended reads. 
The genomic data was analyzed in five main phases:

- Quality check and control (Data Preprocessing)
- Mapping
- Variant Calling
- Phylogenetic analysis.
- Omics analysis

## Bioinformatics Analysis 
 
## Data Preprocessing
The raw data was preprocessed prior to downstream process so as to assess the quality of the reads.

## a) Quality Check

`fastqc` was used to check the quality of the reads, and using `multiqc` the html reports of the reads generated by `fastqc` were aggregated to one multiqc_report.html for easier comparison.

The quality report revealed that some of the reads were overrepresented, there was presence of illumina universal adapters in the reads, and warnings in the GC content plot.

The plot below shows %GC content of the reads. Some reads(orange-coloured curves) have warnings on the GC content which could be a problem with the library or presence of contaminants like adapters.


![fastqc_per_sequence_gc_content_plot(1)](https://user-images.githubusercontent.com/77806488/122923800-bb2b3880-d36d-11eb-9abe-e184747194db.png)
Figure 2a : A plot showing the GC content per sequence for the four samples.


The plot below indicates the presence of illumina universal adapters in most of the sequences.


![fastqc_adapter_content_plot](https://user-images.githubusercontent.com/77806488/122922861-bca83100-d36c-11eb-90c3-d5c9d5d1a2fc.png)
Figure 2b: A plot showing the adapter content in the sequences of the four samples before processing the data.

## b) Trimming

Trimming parameters were then inferred from the [multiqc_report.html](). Using `trimmomatic`, low quality reads, adapter sequences and reads shorter than 25bp in length were removed.

    trimmomatic \
    PE -phred33 \
    -threads 8 \
    ${Reads[0]} \
    ${Reads[1]} \
    $fq_1_paired \
    $fq_1_unpaired \
    $fq_2_paired \
    $fq_2_unpaired \
    SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${adapters}:2:40:15

## c) Post-Quality check
The quality of the trimmed reads were further checked using `fastqc` and a similar [multiqc_report]() generated as before.

The plot below indicates low percentage of sequences with adapters after being trimmed out.


![fastqc_adapter_content_plot(1)](https://user-images.githubusercontent.com/77806488/123067317-8d9bc900-d419-11eb-8139-536a1dc2fde6.png)
Figure 2c: A plot showing the adapter content in the sequences of the four samples after trimming the adapter sequences.

## Mapping

The reads were then aligned to the [*Trypanasoma congolense* reference genome](https://tritrypdb.org/common/downloads/Current_Release/TcongolenseIL3000/) using `bwa mem` aligner.

`bwa index` generated index files of the reference genome that are required for alignment.

Upon mapping, there was low mapping coverage as sample Tryp-B87_S2 had 3.26% while Tryp-T70_S4 had 14.4%.

Low mapping coverage is an indication of contamination with other genomes most probably the host genome either *Glossina pallidipes* or *Bos taurus* where the samples were extracted from.


|Sample|Total reads|Reads that mapped|Mapping coverage|
|------|-----------|-----------------|----------------|
|Tryp-B87_S2|84,838,617|11,995,938|3.26%|
|Tryp-T70_S4|63,611,080|2,074,549|14.4%|

Table 3: Mapping coverage of tsetse and blood samples to the T.congolense IL3000 reference genome before decontamination.

In order to confirm for the presence of other genomes in the reads other than the target genome, `fastqscreen` was employed, which mapped all the reads to the reference target genome, and the suspected host genomes.

`fastqscreen` generated plots showing the proportion of reads that mapped to the reference genomes. *Bos taurus* and *Glossina pallidipes* had a high percent of reads mapping to them while a very small percent mapped to the target genome(*T.congolense*).

The plots below show the proportion of reads that mapped to the reference genomes from `fastqscreen`.


![Tryp-B87_S2pre](https://user-images.githubusercontent.com/77806488/123062610-5b886800-d415-11eb-97c8-abc134a680c7.png)
Figure 3a: Plot showing the proportion of Tryp-B87_S2 reads that mapped to various genomes; Glossina pallidipes, T.congolense, T.vivax, Bos taurus 

![Tryp-T70_S4pre](https://user-images.githubusercontent.com/77806488/123061564-6bec1300-d414-11eb-93ff-23163c0955dd.png)
Figure 3b: Plot showing the proportion if Tryp-T70_S4 reads that mapped to various reference genomes; Glossina pallidipes, T.congolense, T.vivax and Bos taurus before
decontamination

![Tryp-T70_S4](https://user-images.githubusercontent.com/77806488/123058747-cafc5880-d411-11eb-8311-89f719cb40bf.png)
Figure 3c: Plot showing the proportion of Tryp-T70_S4 reads that mapped to T.congolense reference genome after decontamination.

Due to low mapping coverage, the reads were further decontaminated so as to remain with only the target genome(*T.congolense*).They were extracted using `samtools` and were then mapped to *T.congolense* upon which the mapping coverage was 100% for both reads.

**Table 1b: Mapping coverage after decontamination.**

| Sample | Total reads | Reads that mapped |Mapping coverage |
|--------|-------------|-------------------|-----------------|
| Tryp-B87_S2 | 11,814,873 | 11,814,873| 100% |
| Tryp-T70_S4 | 2,040,324 | 2,040,324 | 100% |

Table 4: Mapping coverage of tsetse and blood samples to T.congolense IL3000 reference genome after decontamination

*NB/ There was low mapping coverage as the genome size of the target (T.congolense) is relative smaller in comparison to the host genomes (Glossina pallidipes or Bos taurus)*


|Organism|Genome length in bp|
|--------|-------------------|
|*Trypanasoma congolense*|41,372,041|
|*Glossina pallidipes*|357,332,231|
|*Bos taurus*|2,857,605,192 |

Table 5: A comparison of the reference genome sizes amongst T.congolense, Glossina pallidipes, Bos taurus.

The decontaminated reads proceeded to downstream analysis.

## Merging SAM files

Using `samtools merge` the aligned SAM files were merged to one file for further analysis.

`samtools merge Tryp-B87_S2.sam Tryp-T70_S4.sam Merged.sam`

## SAM to BAM conversion
Since SAM files format are large and occupy more disk space, they are converted to binary formats, BAM, that are smaller in size and easier to manipulate.

`samtools view` was used for conversion.

`samtools view -Sb Merged.sam > Merged.bam`

## Sorting
The bam file was then sorted according to genomic coordinates using `samtools sort`.

`samtools sort -O bam -o Merged_sort.bam Merged.bam`

## Indexing
The sorted file was then indexed using `samtools index`.

`samtools index Merged_sort.bam > Merged_sort.bam.bai`

## Deduplication
Duplicates arise during PCR amplification or sequencing for example optical read duplicates. Removing duplicates is essential so as to minimise the occurences of False Positive variants during variant calling arising from duplicate reads.

Duplicates were removed using `gatk Markduplicates` tool which in addition to removing duplicate reads, it also generates a metrics file containing statistics about deduplication.

`gatk Markduplicates --INPUT Merged_sort.bam --OUTPUT Dedups.bam\
--METRICS_FILE metrics.txt --REMOVE_DUPLICATES true`

A [metrics file](https://github.com/mbbu/AAT-project/blob/dev/Results/Deduplication/marked_dups_metrics.txt) indicating the number of duplicate reads removed was generated.


|Library|Unpaired_read duplicates|Read pair duplicates|Optical read duplicates|Percent_duplication|
|-------|-----------|------|------|-----|
|Tryp-T70_S4|40|348014|633|0.059148|
|Tryp-B87_S2|12|35292|57|0.034608|

Table 6: Deduplication metrics summary showing the tsetse samples (T70) had a 0.06% rate of duplication while the blood sample had 0.03% rate of duplication.


## Creation of Sequence Dictionary
`gatk HaplotypeCaller` was the variant caller. It requires a sequence dictionary file, `ref_genome.dict` and `ref_genome.fa.fai` in addition to the reference genome file in fasta format. These auxiliary files contain summary information about the lengths of the chromosomes and where they start.

The sequence dictionary file of the reference genome was created using `gatk CreateSequenceDictionary` tool.

`gatk CreateSequenceDictionary -R T_congolense.fasta`

While the `ref_genome.fa.fai` was created using `samtools faidx`.

`samtools faidx T_congolense.fasta`

## BaseRecalibration
Base Quality Score Recalibration(BQSR) step is recommended by GATK development team as part of the [best practices for Variant calling with GATK](https://www.broadinstitute.org/partnerships/education/broade/best-practices-variant-calling-gatk-1).

BQSR assigns accurate confidence scores to each sequenced base as sequencing machines sometimes assign inaccurate and biased scores which is a major contributor to bad calls during variant analysis.

It involves two steps: first `gatk BaseRecalibrator` tool builds a model based on the `input.bam` file, a set of `known variants` and a `reference genome` hence producing a recalibration file.

`gatk BaseRecalibrator -I marked_duplicates.bam --known-sites knownvariants.vcf -R T.congolense.fasta -O recal.table`

A [recalibration report](https://github.com/mbbu/AAT-project/blob/dev/Results/baserecalibration/recalibration_table) generated after the first step. It contains 5 tables; Argument, Quantization, ReadGroup, Quality Score and Covariates table.

 `knownvariants` file was generated after running the pipeline twice, first time uncalibrated then the second run calibrated using the variants file from the first run. The `knownvariants` file masked bases at sites of expected/real variants so as to avoid counting them as errors.

 The second step in BQSR employs `gatk ApplyBQSR` tool in adjusting the base quality scores of the deduplicated bam file based on the previously generated recalibration table model.

 `gatk ApplyBQSR -R T.congolense.fasta -I marked_duplicates.bam -bqsr recal.table -O recalibrated.bam`


## Variant calling
`gatk HaplotypeCaller` was used in identifying variants in the samples.

The input file is the `recalibrated.bam` file and we specified our output to be in a `.vcf format`. The reference file also needs to be specified for haplotype realignment.

A total of `382950` variants were called using `gatk HaplotypeCaller`.

`gatk HaplotypeCaller -I recalibrated.bam -O variantsGHC.vcf -R T.congolense.fasta`


## Hard-filtering low quality variants

Ideally, `GATK` recommends Variant Quality Score Recalibration (VQSR) as a subsequent step for filtering variants. However, VQSR requires highly validated variant resources such as HapMap datasets to select true positives. These resources are only tailored for human data variant calling practices.

In the absence of those resources, we resorted to hard-filtering so as to remain with high-confident variants.

The [recommended](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) way of getting hard-filtering parameters for one's dataset is to extract annotation values of the variants using `gatk VariantsToTable` from the `.vcf` file then using Rstudio to plot and visualize the annotations.  

`gatk VariantsToTable -V variantsGHC.vcf -F QD -F MQ -F SOR -F FS -F MQRankSum -F ReadPosRankSum -O variants.table`

The above code shows the annotation values which are QualitybyDepth(QD), MappingQuality(MQ), StrandOddsRatio(SOR), FisherStrand(FS), MappingQualityRankSumTest(MQRankSum), and ReadPosRankSumTest(ReadPosRankSum) from the `.vcf file` that were extracted.


| QD | MQ | FS | SOR | MQRankSum | ReadPosRankSum |
|----|----|----|-----|-----------|----------------|  
| 34.62| 60.00 | 0.000 | 0.527 | 0.000 | 1.465 |
|32.92   |60.00   |0.000   |0.330   |0.000|0.887|
|7.79|    53.27|   6.117|   3.003|   -3.028|  -2.112|
|8.31    |52.82   |2.869   |2.712   |-2.928  |-2.680|
|8.51|    52.82|   2.869|   2.712|   -2.928|  -2.802|
|3.80    |55.70   |3.282   |2.276   |-2.275  |-1.878|
|33.49|   49.71|   13.554|  5.283|   -2.546|  -2.448|
|25.36   |48.87   |14.065  |5.421   |-2.591  |-1.922|
|28.73|   47.15|   15.441|  5.784|   -2.926  |-2.189|
|30.97   |47.81   |21.875  |5.782   |-3.208  |-1.345|

Table 7: Called variants annotation values that were used to set plot annotation plots for variant filtering

Visualization of the annotation values was done using Rstudio and the cutoffs threshold were determined by comparing our annotation plots to plots generated by VQSR during filtering.

**Hard-filtering parameters**

1. QualitybyDepth(QD)

QD normalizes variant quality so as to avoid inflation of QUAL scores in cases of deep coverage, thus it is better to QD scores for filtering rather than QUAL(quality scores) and DP(Depth) directly.

For our dataset we used filtered out low quality variants with a QD of less than `6.0` as seen from the QD density plot below.


![QDplot](https://user-images.githubusercontent.com/77806488/123109753-2a735c00-d444-11eb-9d21-74e85abeda8c.png)
Figure 4a: Density plot showing the quality by depth in called variants with most variants with a QD of around 30

2. MappingQuality(MQ)

This annotation value gives the mapping quality of the variants to the reference by taking into account the standard deviation of the mapping qualities. When mapping qualities are good at a site the MQ will be around `60`.

From our dataset, most of our variants as seen in the figure below had a mapping quality of `60` and few variants had a mapping quality of less than `40`. Thus we filtered any low confidence variants with `MQ < 40`.

![MQplot](https://user-images.githubusercontent.com/77806488/123109882-48d95780-d444-11eb-8cda-9ef7c5acaac0.png)
Figure 4b: Density plot showing the mapping quality of called variants with most variants clustering at a mapping quality of 60 and fewer variants at 40.

3. FisherStrand(FS)

FS value indicates strand bias and whether the alternate allele was seen more or less often in the forward or reverse strand than the reference allele. FS deteremines the strand bias based on the Fisher's Exact Test. Ideally, when there is little to no strand bias the FS value is close to `0`.

As from the figure below, we converted the x-axis to log-scale so as to easily note the distribution using FS. As observed most variants had an FS value less than `10` indicating less strand bias, we however set our FS value to filter out any variants which had an FS value greater than `50.0` as they were low confidence variants and most likely false positive.

![logFS](https://user-images.githubusercontent.com/77806488/123110101-74f4d880-d444-11eb-9e36-488caebedead.png)
Figure 4c: Density plot of the fisher strand value of the called variants indicating the strand biasness at a specific site. Most variants had an FS value of <10,
showing less strand biasness

4. StrandOddsRatio (SOR)

SOR also indicates strand bias just like FS values, however it uses a test similar to the symmetric odds ratio test. FS penalizes variants that occur at the ends of exons as reads at ends of exons are mainly covered by reads in one direction. SOR comes in handy, by taking into account the ratios of reads that cover both alleles.

We filtered out low quality variants with an SOR value greater than `3.75` implying that the ratios of the reads covering both alleles at the end of the exons is highly biased. The figure below indicates density plot using SOR values.

![SORplot](https://user-images.githubusercontent.com/77806488/123110224-8c33c600-d444-11eb-9b0e-b7595d154076.png)
Figure 4d: Density plot of Strand Odd Ratio values of called variants which also assesses the strand biasness on each site with most variants clustering at an SOR 
value < 2.5

5. MappingQualityRankSumTest(MQRankSumTest)

MQRankSumTest determines mapping qualities according to the u-based z-approximation from the RankSumTest. It compares the mapping qualities of the reads supporting the alternate and reference alleles. Ideally, to depict little to no difference in the mapping qualities of reads supporting either alleles the MQRankSumTest should be close to `0`.

Even though most of our variants as shown below had a MQRankSumTest close to zero, we filtered out any variants with the MQRankSumTest less than `-5.0` .

![MQRankSum](https://user-images.githubusercontent.com/77806488/123110316-a4a3e080-d444-11eb-8e6f-dd54a9cf7e86.png)
Figure 4e: Density plot of mapping quality ranked sum test values of called variants where most variants have a value close to 0 showing there is little
to no difference in the mapping qualities of reads supporting either of the alleles.

6. ReadPosRankSumTest(ReadPosRankSum)

ReadPosRankSum uses u-based-z-approximation from the Rank Sum Test for site positions within the reads. It compares if the positions of the reference and alternate alleles are different within the reads. Little difference between the positions of the reference and alternate alleles is indicated by values close to `0`.

We filtered out any low quality variants that had a ReadPosRankSum less than `8.0` as it indicates the alternate alleles were found at the ends of the reads more often than the reference alleles.


![ReadPosRankSum](https://user-images.githubusercontent.com/77806488/123110403-bb4a3780-d444-11eb-8c53-7a150c9a1df0.png)
Figure 4f: Density plot of ReadPosRankSum Test values of called variants where most variants have a value close to 0 showing there is little difference between the
positions of the reference and alternate alleles.

**Filtering low confidence variants**

`gatk VariantFiltration` tool was used in filtering out the low confidence variants(FPs) using the hard-filtering parameters determined after visualizing Rplots.

    gatk VariantFiltration -V variantsGHC.vcf \
    -filter "QD < 6.0" --filter-name "QD6" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "FS > 50.0" --filter-name "FS50" \
    -filter "SOR > 3.75" --filter-name "SOR3.75" \
    -filter "MQRankSum < -5.0" --filter-name "MQRankSum5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O variants_tagged.vcf

Low confidence variants that didn't meet the cutoff thresholds were tagged by the filter-name and using `gatk SelectVariants` variants that satisfied all the filtering parameters and had the `PASS` tag in the filter column were selected.

`gatk SelectVariants -V variants_tagged -O variants_filtered -select 'vc.isNotFiltered()'`

`295,087` high confidence variants were left after filtering.

## Normalization and Decomposition.
The filtered VCF file needed to be preprocessed for annotation. This preprocessing involved decomposing and normalizing multiallelic variants that might arise after joint-calling. Decomposing expands multiallelic variants into distinct variant records such that there is one record for each REF/ALT combination. Normalization then aligns ALT and REF alleles and removes any suffixes and prefixes that are shared by the alleles.

301945 variants were identified after decomposition and normalization. This represented the total number of variants from the two samples which were Trypanosoma
congolese collected from host cattle blood and the Tsetse fly vector. ```Bcftools``` was used to split the multisample vcf file to individual sample vcf files
for comparison purposes. Comparison of the two files revealed a large disparity in the resulting variants from the two samples which was not expected considering the
nature of sample collection (the samples were obtained at the same time in regions of close proximity).

This difference was hypothetically labelled as a case of low parasitemia in the blood sample that was collected since it would have meant there was less chance of
obtaining high quantity and quality of *Trypanosoma congolese* from the host.

| Sample | Normalized Variants |
|--------|-------------|
| Tryp-B87_S2 | 6521 |
| Tryp-T70_S4 | 295240 |
| Common variants|315|

Table 8 : A summary of the total variants called in blood sample (Tryp-B87_S2) and Tsetse sample (Tryp-T70_S4) as well as common variants between the two samples.

## Annotation
Annotation was done using the `snpEff` annotater which allows for building off the reference genome into the configuration file for annotation analysis.

From the analysis of the variants;
1. There was a high percentage of high impact variants from the blood sample than the tsetsefly sample.

| Sample | High Impact variants % |
| ------- |---------------------- |
| Tryp-B87_S2|6.314|
|Tryp-T70_S4|0.704|

Table 9: A summary of the high impact variants in blood sample (Tryp-B87_S2) and tsetse sample (Tryp-T70_S4) after using snpEff annotator.

2. Tryp-B87_S2 sample had a high ratio of missense to silent variants as compared to Tryp-T70_S4.

| Sample | Missense / Silent Ratio|
|--------|------------------------|
|Tryp-B87_S2|2.8723|
|Tryp-T70_S4|0.9554|

Table 10: A summary of the missense: silent ratio in blood sample (Tryp-B87_S2) and tsetse sample (Tryp-T70_S4) after using snpEff annotator.

These observations on the high impact variants in the blood as well as a high missense to silent ratio imply that these variants effect by impact and functional class in the blood sample could be enabling the parasite adaptation in evading the reservoir's immune response.
In light of the mixed-type variants present, SNPs were extracted from the `.vcf` files using GATK `SelectVariants`. The variants of interest were those annotated to have
HIGH and MODERATE impact which were extracted from the SNP-containing `.vcf` files on sample basis (separate files for the blood and Tsetse fly samples).Predicted gene products were added with the aid of the Gene Ontology Association File `.gaf` file.
For the blood sample,a summary of the number of genes associated with a particular variant type and their gene product can be found [here](https://docs.google.com/spreadsheets/d/1T2gw39u6QEzc-7KD21ZFPunV-7cLDx12nJyf3qD08e8/edit?usp=sharing) 

## Trypanosoma vivax samples variant calling for comparative analysis
Due to the large disparity in the number of variants called from the blood and Tsetse fly trypanosomal samples, we decided to run a parallel variant calling pipeline
pipeline with the available *Trypanosoma vivax* sample datasets( which were sampled at the same time as the Trypanosoma congolese datasets) as the input reads. This
was done so as to develop a validatory and comparative framework with respect to the large unexpected difference in variant count from the *Trypanosoma congolese*
variant calling pipeline.

The preprocessing procedures that are highlighted at the start of the report were replicated to the *Trypanosoma vivax* samples labelled as `Tryp-B151_S1` and
`Tryp-T8_S3`. This involved screening of potential contaminants and decontamination of the read. The decontaminated reads mapped with 100% coverage against the
Trypanosoma vivax reference genome, `TriTrypDB-52_TvivaxY486_Genome.fasta` The reads then proceeded to the subsequent steps and variant were called against the same
reference genome that was utilized for decontaination. Hard-filtering parameters were also set using graphical models as described previously after which annotation
was done using the snpEff within which the *Trypanosoma vivax* database was built (in the respective configuration file).
The pipeline was run twice for base recalibration purposes to allow for convergence of high confidence variants.

The resulting variant count was as follows:

| Sample | Raw Variants | Filtered variants | Normalized variants |
|--------|-------------|-------------------|-----------------|
| Tryp-B151_S1 | 87,760 | 85,318 | 85,395 |
| Tryp-T8_S3 | 36,504 | 28,050 | 28,456 |

Table 11: A summary of the total variants called in blood sample (Tryp-B151_S1) and tsetse sample (Tryp-T8_S3) infected with T.vivax

As indicated, the results exhibited an inverse variant count when compared to the *Trypanosoma congolese* variants. The presence of significantly more variants
in *Trypanosoma vivax* samples that were extracted from blood than those collected from Glossina pallidipes could be attributed to the fact that the protozoan
parasite's life cycle. The life cycle of Trypansoma parasites within the fly is characterized by three key phases:

1. The ingestion of bloodstream-form parasites.
2. Migration of parasites to the insect's midgut with differentiation into procyclic forms.
3. Migration anteriorly to the proboscis (T. vivax), proventriculus (T. congolense) or salivary gland (T. brucei) and differentiation into epimastigote forms.
3. Differentiation into metacyclic form and innoculation into the vertebrate host upon insect-feeding.

It is however key to note that the *Trypanosoma vivax* life cycle lacks the migratory step which involves movement of parastic forms to the midgut of the Tsetse fly.
This implies that the parasite spends a realtively shorter time in the vector and hence may not acccumulate numerous polymorphic traits while in the insect, possibly
as compared to the genetic variation that would arise when it is present in the infected host. This could then serve to explain the difference in variant count
obtained from the *Trypanosoma vivax* samples.

In a comparative view to the Trypanosoma congolese variant count disaprity, this run confirmed that the decontaminated reads from both T.congo and T.vivax samples
were of good quality annd could hence be used to provide biological inferences that are aligned to the study. This is because the samples were collected at the same
time, extracted using the same kit and sequenced using the same platform (Illumina Hiseq).
THe disparity in the number of T.congo variants called could then possibly be attributed to:

1. Low parasitemia level in the blood sample - varyting parasite load rom the two sample collection sources (host and vector) could have had a direct effect in
   the resulting variant count which hence may cause a wide disarity

2. Multiple genotypes - due to the nature of the Tsetse fly vectors to feed on several hosts, they may ingest different Trypanosoma strains which may not be not have
   been captured when variants were called solely against the Trypanosoma congolese reference genome.



## Taxonomic classification of raw sample datasets
To establish the source of the contaminants that were present in the raw sample datasets, taxonomic classification analysis was done at a metagenomic level.
Metadata and prior decontamination information revealed *Bos taurus* (cattle) and *Glossina pallidipes* (Tsetse fly) as the major contaminants. To assess the presence of
any other possible organism that would have formed the contamination fraction of the samples, classification was done using Kraken which is described by a system for
assigning taxonomic labels to short DNA sequences. Kraken allows for the use of various database sets but since we were interested in a metagenomic based analysis,
the MiniKraken prebuilt database which is constructed from complete bacterial, archeal and viral genomes.

Taxonomic analysis relevaled microbial contamination in the T.congo sample extracted from the Tsetse fly vector, however no microbial contamination was identified
in the T.congo blood sample which was expected since microbe load in blood signifies an infection.

The microbial species identified from the T.congolese raw data samples were as follows:

| Species | Contamination coverage |
|--------|--------------|
| Glossina pallidipes salivary gland hypertrophy virus | 14% |
| Sodalis glossinidius | 5% |
| Alteromonas mediterranea | 2% |

Table 12: A summary of the microbe species identified from T.congolense infected samples after carrying out a metagenomics analysis using the MiniKraken prebuilt
database

These microbial species are majorly endosymbionts which have developed a mutualistic relatioship with the Tsetsefly vector and therefore are expected contaminants.

![snapshot](https://user-images.githubusercontent.com/68329457/124580204-94363180-de58-11eb-93b7-9e2115c30ee8.png)
Figure 5: A pie-chart representation of the microbial contamination distribution in Sample Tryp-T70_S4 generated by Kraken using the MiniKraken prebuilt database.

## Phylogenetic analysis of Trypanosoma congolese samples
To create a distinct picture of how the T.congolese samples collected from Shimba Hills National Park would cluster based the nature of the variants that were called,
we sourced T.congolese datasets that were made available from this [publication](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14271) which served to analyze hybridization between genetic lineages of Trypanosoma congolese.
A summary of the samples that were included for phylogenetic analysis are as follows:

| Sample Identity | Country of origin | Host organism |
|-----------------|------------------|----------------|
| ERR1898732 | Congo | Cattle |
| ERR1898727 | Zambia | Cattle |
| ERR1898728 | Burkina Faso | Cattle |
| ERR1898764 | Zambia | Cattle |
| ERR1898765 | Zambia | Cattle |
| Tryp-B87_S2 | Kenya | Cattle |
| Tryp-T70_S4 | Kenya | Tsetse Fly |

Table 13: A summary of the sample used, location and sample-type from determined dataset for subsequent phylogenetic analysis.

The sample reads served as input to the variant calling pipeline which has been described initially in the report, that is, the GATK Workflow where the variants were
called against the Trypanosoma congolese IL3000 reference sequence. A breakdown of the variants discovered at raw and normalized level are as follows


| Sample | Raw variants | Normalized variants |
|-----------------|------------------|----------------|
| ERR1898732 | 264,491 | 227,993 |
| ERR1898727 | 428,565 | 346,795 |
| ERR1898728 | 764,538 | 578,339 |
| ERR1898764 | 597,943 | 470,172 |
| ERR1898765 | 611,732 | 479,313 |
| Tryp-B87_S2 | 3440 | 1840 |
| Tryp-T70_S4 | 356,839 | 299,094 |

Table 14: A summary of the raw and normalized variants called from samples derived from determined datasets for phylogenetic analysis.

The `.vcf` files containing the discovered variants were then preprocessed in preparation for phylogenetic analysis. The files were converted to the preferred PHYLIP
formart using `vcf2phylip`. The compatible files were then passed into BEAUTi (Bayesian Evolutionary Analysis Utility) which is a program for setting model parameters
for BEAST (Bayesian Evolutionary Analysis Sampling Trees). In BEAUTi, a HKY substitution model was used with base frequencies estimated. An uncorrelated relaxed clock was set
with an exponential relaxed distribution to give sense to the evolution rate. We also used a Markov Chain Monte Carlo(MCMC) chain length of 10000000 for posterior distribution
approximation. BEAUTi returned an XML command line file which serves as the input to BEAST which provides output log files. The phylogenetic tree was visualized by loading the
`.trees` output file to `FigTree`.


![phylotree](https://user-images.githubusercontent.com/68329457/125045706-f3d54c80-e0a5-11eb-916d-e4b196613d4f.jpg)

Figure 6a: Rooted phylogenetic tree visualized using FigTree showing the relationship of variants called in our samples from Shimba Hills Park and 
additional samples from Zambia, Congo and Burkina Faso.

For confirmatory purposes, we inferred evolutionary relationship using another tool based on maximum likelihood estimation. SNPs which were to be utilized for inferring phylogenetic relationship were extrracted from the multisample vcf file using the GATK `SelectVariants` feature.The SNP-containing vcf file was then converted to PHYLIP format using `vcf2phylip`. The `.phy` file then served as input to PHYML which was the tool that was used to estimate maximum likelihood phylogeny from the called SNPs.With respect to parameters set, a HKY85 substitution model was used with an estimated Ts/tv ratio as well as an estimated proportion of invariable sites and gamma distribution. A boostrap value of 1000 was also utilized to ensure a high confidence phylogenetic estimation from the datasets provided.


![phymltree](https://user-images.githubusercontent.com/68329457/125953239-011bd0e4-203d-4256-ada7-fc9d38c034b8.png)
Figure 6b: Rooted phylogenetic tree visualized using PhyML showing the relationship of variants called in our samples from Shimba Hills Park (red) and additional
samples from Zambia, Congo and Burkina Faso(black)


The inferences made from the generated evolutionary tree are:
1. The Tryp-B87_S2 sample obtained from cattle blood in Shimba Hills is hypothesized to have portrayed false divergence due to a seemingly large difference in variant
count. It is also possible that the sample acts as an ancestor to other lineages and has some shared genetic traits as shown in the tree.

2. ERR1898764 and ERR1898765 as well as ERR1898727 and ERR1898728 samples that were obtained from cattle blood clustered together on the same clade as they were both retrieved from the same
geographical region and exhibit a high level of genetic similarity; conlusions that are shared by the authors of the publication as well.

3. Occurrence of Tryp-T70_S4 (Tsetsefly sample from Kenya) and ERR1898732 (blood sample from Congo) on the same clade indicates a level of similarity which would then serve to further imply
that sylvatic transmission events are possible over a given geographical region triggered by events such as migration.  


## Omics Analysis

Gene ontology and enrichment analysis was carried out using the [TriTrypDB](https://tritrypdb.org/tritrypdb/app) database (Release 53 - 21st July 2021). A SnpEff annotater feature
known as `SnpSift` so as to categorize variants based on variant impact identifiers which were HIGH, LOW, MODERATE and MODIFIERS as well as sample type
(blood sample, tsetse fly sample variants and common variants). Custom bash scripts were utilized to extract gene IDs from variant impact files that were generated.  These list of genes organized on variant impact basis files then served as input to the gene ontology and enrichment analysis feature on the TriTrypDB platform. Gene products for the blood and tsetse fly sample as well common genes  were retrieved . Enrichment results based on cellular component, molecular function and biological processes were also retrived based on the categories of interest. 
The retrieved tab separated files were organized in Google sheet spreadsheets with the goal of establishing similarity and differences between the two samples and the implications of the same. Significant ontology results were obtained using p-adjusted values that were provided by the database (Benjamini values) 
A Custom R script was utilized to identify shared genes between the blood and tsetse samples after which separate enrichment analysis was run to establish shared molecular and biological processes.


Annotation and prediction of the variant effects were performed using SnpEff. The predicted effects of the variants were categorized by impact as high (disruptive impact on the protein), moderate (non-synonymous substitution), low (synonymous substitution) and modifiers (an impact on the non-coding region). To get a better understanding of the variants, we categorized them based on impact, functional class and region of variants called in T. congolense positive samples and summarized these results in the tables below.

| Type | Count | Percent |
|------|-------|---------|
| High | 9,470 | 0.833% |
| Low  | 50,715 | 4.461% |
| Moderate | 50,677 | 4.458% |
| Modifier | 1,025,893 | 90.248% |

Table 15a : A summary of T.congolense variants based on the impact of its effects in comparison to the reference genome

| Type | Count | Percent |
|------|-------|---------|
| Missense | 48,684 | 48.625% |
| Nonsense | 738 | 0.737%
| Silent | 50,699 |	50.638% |

Table 15b: A summary of T.congolense variants based on the impact of its effects in comparison to the reference genome

As seen in table 15 (a) 90% of the T. congolense variants are modifiers, while 0.8% are high impact variants. This means that most of the variants have an impact on the non-coding region in the T. congolense genome.  In addition to this, table 15 (b) reveals that 51% of the T. congolense variants are silent mutations while 49% of variants are missense mutations. A missense mutation is a point mutation that results in coding for different protein while a silent mutation does not affect the sequence of the amino acid. Collectively, most of the T. congolense variants could result in changes in amino acids that could change proteins in the non-coding regions in the genome.

![Types of variants in all samples](https://user-images.githubusercontent.com/68329457/128638491-c04677a4-b23a-488f-a366-9e60d56c9d23.png)
Figure 7(a): A summary of the variants called in T. congolense positive blood and tsetse samples indicating a higher percentage of variants in the upstream (37%) and downstream (36%) regions of most genes

The highest number of variants were in the upstream (37.2%) and downstream (36.3%) regions of most genes in both tsetse and blood samples. This was followed by intergenic (16.6%), synonymous (4.5%) and missense (4.3%) variants. The blood samples had a higher number of variants in the frame shift (5.3%) region than in tsetse samples (0.6%) as shown in Figure 7 (b). Additionally, 2.4% of the total blood variants were responsible for a gain in stop codons.

![Type of Variants in B87](https://user-images.githubusercontent.com/68329457/128638515-88a70499-eb56-4974-b55b-e1d715f2f768.png)
![Type of variants in T70](https://user-images.githubusercontent.com/68329457/128638522-406e9801-9191-46b8-afb8-f67bc060eba3.png)
Figure 7(b): A summary of variants called in blood and tsetse samples respectively showing the regions and impact the variants have when compared to the T. congolense IL3000 reference genome

### Enrichment Analysis
 a) Low Impact
 Variants which were analyzed under the low impact category were synonymous variants that are harmless in nature and are unlikely to change protein behavior. The common genes that were least affected by variation were involved in cilia formation as shown in Figure 8(a) below, which implies low modifications in motile ability as well as cytoskeleton synthesis.
 
 ![Functional enrichment analysis of low impact tsetse variants](https://user-images.githubusercontent.com/68329457/128638575-bd6792d1-8c2c-4b9e-9c8d-cb02ba8c461a.png)
Figure 8(a): A summary of the functional enrichment analysis on low impact variants derived from tsetse samples with a p-adjusted value (padj) <0.05

This is further supported by the presence of microtubule associated complex, ciliary pocket collar, mitochondrial proton-transporting ATP synthase complex, ATPase activity and ATP-dependent microtubule motor activity. The presence of ATP associated compounds are necessary in the production of energy used in for the parasite’s motility. Therefore, membrane-bound organelles such as the mitochondria, Golgi apparatus and endoplasmic reticulum are expected to be enriched in this analysis.  ATP-binding cassette (ABC) transporter was the only significantly enriched low-impact variant in the blood variants with a padj of 0.03. ABC transporters are a large superfamily of membrane proteins involved in diverse function in the parasite such as translocating a wide range of substrates across membranes against a concentration gradient in an ATP-dependent process 

![Functional enrichment analysis of synonymous tsetse variants](https://user-images.githubusercontent.com/68329457/128638599-e0e1a183-168f-4a48-960b-18ac5256b6a1.png)
Figure 8(b): A summary of the functional enrichment analysis on synonymous variants derived from tsetse samples with a p-adjusted value (padj) <0.05

Synonymous variants were particularly investigated to paint a vivid picture of the genes that were least affected by variation. Enrichment results at molecular function level revealed that the synonymous variants were localized in genes which give rise to enzymes including nucleases, transferase, synthase and involved in protein and mRNA binding. They also occurred in genes that facilitate motor activity. This further supports that low impact variants are indeed synonymous mutations, as these are essential functions in the parasite to survive in the tsetse fly.

b) High Impact Variants
These include highly disruptive variants such as frameshift mutations as well as stop-lost and stop-gained variants. Variants which are highly disruptive were observed to majorly affect distinct cell organelles. 

![Functional enrichment of high impact tsete variants](https://user-images.githubusercontent.com/68329457/128642451-f1e77c55-a3ff-4faa-8275-dfd353847eb0.png)
Figure 8(c): A summary of the functional enrichment analysis on high impact variants in biological processes (green) and cellular components (blue) derived from tsetse samples with a p-adjusted value (padj) <0.05

Majority of the high-impact variants in tsetse samples are involved in immune responses to the parasite in the tsetse fly with the presence of response to host immune response, evasion of host immune response, response to biotic stimulus in the biological process category of the GO enrichment analysis in Figure 8(c) above. As the parasite migrates from the midgut to the mouth part of the tsetse fly, it encounters a range of immunological and physical barriers that are by-passed by Variant Surface Glycoproteins (VSGs) as well as cell-surface proteins such as Glutamic acid Alanine-rich protein (GARP) (Tetley and Vickerman, 1985, Hehl, 1995). VSGs coat the surface of trypanosomes and are expected to be high-impact variants as they are sequentially replacing their glycoprotein surface coat through antigenic variation to evade the immune responses in the vector.to evade to the immune response (Horn, 2014).

c) Moderate Impact Variants
Moderate impact variants are termed to be non-disruptive but could change the effectiveness of the resulting protein. These may include missense mutations and in frame deletions. The overall analysis on these variants is seen in Figure 8 (d) which indicates an involvement of the genes in the immune response in the tsetse fly. These are similar observations seen in the high impact variants in Figure 8 (c). When the missense variants were analyzed for their functional enrichment, similar observations were seen showing that most of the moderate impact variants are missense mutations.

![Functional enrichment of moderate impact tsetse variants](https://user-images.githubusercontent.com/68329457/128642481-ce0e53bf-837d-4571-8408-d21295386eae.png)
Figure 8(d): A summary of the functional enrichment analysis on moderate impact variants in biological processes (green) and cellular components (blue) derived from tsetse samples with a p-adjusted value (padj) <0.05

d) Modifying variants
Variants in this class are known as non-coding variants or rather variants which are present in non-coding genes. These are genes that are exons, located downstream (occurs in the 3’ region of a gene) and upstream (occurs in the 5’ region of a gene). An overall representation of modifying variants in T. congolense samples have been summarized in Figure 8 (e) below. The figures show a similarity in the cellular component significantly enriched in both blood and tsetse samples which include the plasma membrane, intracellular membrane, cytoplasm, and the nuclear lumen. The blood variants have enriched biological processes in the nucleus as seen with the presence of cellular nitrogen compound metabolic process, nucleobase-containing compound metabolic process and heterocycle metabolic process. This indicates specific activities are targeted in the parasites in the blood samples in comparison to the tsetse samples.

![Functional enrichment of modifying tsetse variants](https://user-images.githubusercontent.com/68329457/128642528-577f4fc9-bdd3-4f9c-8eb5-934de3937613.png)
![Functional enrichment of modifying blood variants](https://user-images.githubusercontent.com/68329457/128642546-cbc80717-c30e-4a6b-b35f-d61d082611cd.png)
Figure 8(e): A summary of the functional enrichment analysis on modifying variants in biological processes (green) and cellular components (blue) derived from tsetse and blood samples with a p-adjusted value (padj) <0.05

On further investigation, the most enriched metabolic pathways in these variants (cut-off padj <0.05) were GDP-mannose biosynthesis, β-1,4-D-mannosyl-N-acetyl-D-glucosamine degradation and β-(1,4)-mannan degradation which are seen to be specifically enriched in the bloodstream forms of the parasite in Trypanosoma brucei (Kuettel, 2012). These sugar nucleotides have been shown to be essential for the survival and for infectivity of trypanosomatid parasites (Roper et al., 2002; 2005; Urbaniak et al., 2006-Kuettel 2012).

![Functional enrichment analysis of tsetse downstream variants](https://user-images.githubusercontent.com/68329457/128642621-9d2639d0-1f4a-4fd7-8702-8a11d21abb56.png)
Figure 8(f): A summary of the functional enrichment analysis on downstream variants in biological processes (green) and cellular components (blue) derived from tsetse samples with a p-adjusted value (padj) <0.05

In the tsetse samples, downstream gene variants were observed to occur in genes whose products were localized in the nucleoplasm, nucleolus, axoneme, phagosome and spliceosomal complex as shown in Figure 8(f) These protein products are distinctively involved in cell motility, endocytosis as well as RNA splicing events.
In the upstream variants derived from tsetse samples are mostly enriched in the regulation of gene expression and metabolic processes in various organelles such as the ribosomes, cilium, cytoskeleton, nucleoplasm, mitochondria, and intracellular membrane as seen in figure 7 (f). MetaCyc and KEGG were used to carry out a metabolic pathway enrichment analysis indicating superpathway of seleno-compound metabolism, NAD/NADH phosphorylation and dephosphorylation, glycine, serine and threonine metabolism, pyruvate metabolism and gluconeogenesis pathways as most enriched in these genes.

![Functional enrichment analysis of tsetse upstream variants](https://user-images.githubusercontent.com/68329457/128642655-2648bceb-3aa9-4ca5-9f58-2260e2cb0dbe.png)
Figure 8(g): A summary of the functional enrichment analysis on upstream variants in biological processes (green) and cellular components (blue) derived from tsetse samples with a p-adjusted value (padj) <0.05

![Functional enrichment analysis of upstream blood variants](https://user-images.githubusercontent.com/68329457/128642674-e6345a5d-ec8e-40d7-9ada-72fc298bc622.png)
Figure 8(h): A summary of the functional enrichment analysis on upstream variants in biological processes (green), molecular functions (red) and cellular components (blue) derived from blood samples with a p-adjusted value (padj) <0.05

Most upstream variants derived from the blood sample occur in genes whose products are predicted to be localized in cytoplasm, nucleus, and mitochondria and are characterized by various enzymes as well as RNA binding activity as shown in Figure 8(h) above. 
Intron variants (variants located within an intron site were localized in fewer genes as compared to upstream gene variants. It is also interesting to note that variations at the intron region were only identified from the tsetse sample. The genes in which the variations occurred enabled the development of cellular organelles, binding of DNA, ATP, and sno(s)RNA (small nucleolar RNA) which are essential for the modification and processing of ribosomal RNA precursors







---
