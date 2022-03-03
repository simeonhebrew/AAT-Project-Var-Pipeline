// Enable DSL 2 syntax

nextflow.enable.dsl = 2

// Process 1: Quality control analysis --tool: fastqc

process QUALITY_CHECK {
    publishDir path: "${params.outdir}/Qc", mode: 'copy'
    tag "Quality Checking:"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_logs"

    script:
    """
    mkdir ${sample_id}_logs
    fastqc -o ${sample_id}_logs -f fastq -t 8 -q ${reads}
    """
}

// Process 2 : MultiQC Summary Report --tool: multiqc

process RAWMULTIQC {
	publishDir path: "${params.outdir}/multiqc_untrimmed", mode: 'copy'
	tag "MultiQC report for raw reads:"
	
	input:
	file(fastqc)

	output:
	file('multiqc_report.html')

	script:
	
	"""
	multiqc . 

	"""
}

// Process 3 : Trimming reads  --tool: trimmomatic

process TRIMMING {
	publishDir path: "${params.outdir}/trimmed", mode: 'copy'
	tag "Trimming reads:"

	input:
	tuple val(sample_id), path(Reads)
	path(trim) 
	path(adapters)

	output:
	tuple path(fq_1_paired), path(fq_2_paired), emit: paired_reads

	script:
	fq_1_paired = sample_id + '_R1.paired.fastq'
	fq_1_unpaired = sample_id +'_R1.unpaired.fastq'
	fq_2_paired = sample_id + '_R2.paired.fastq'
	fq_2_unpaired = sample_id + '_R2.unpaired.fastq'

	"""
	trimmomatic \
	PE -phred33 \
	-threads 4 \
   	${Reads[0]} \
    	${Reads[1]} \
    	$fq_1_paired \
    	$fq_1_unpaired \
    	$fq_2_paired \
    	$fq_2_unpaired \
    	SLIDINGWINDOW:${trim} MINLEN:25 ILLUMINACLIP:${adapters}:2:40:15 

    	"""

}

// Process 4 : Quality checking trimmed reads --tool: FASTQC
process TRIMFASTQC {
    publishDir path: "${params.outdir}/trim_fastqc", mode:'copy'
    tag "Running FASTQC on trimmed reads:"

    input:
    tuple path(Read_R1), path(Read_R2)

    output:
    path "${sample_id}_logs"

    script:
    sample_id = ( Read_R1 =~ /(.+)_R1.paired.fastq/ )[0][1]
    """
    mkdir ${sample_id}_logs
    fastqc -o ${sample_id}_logs -f fastq -t 8 -q ${Read_R1} ${Read_R2}

    """
}

// Process 5 : MultiQC on trimeed reads FASTQC reports --tool: MultiQC
process TRIMMULTIQC{
    publishDir path: "${params.outdir}/multiqc_postfastqc", mode:'copy'
    tag "MultiQC report after trimming:"

    input:
    path(fastqc_out)

    output:
    file('multiqc_report.html')

    script:
    
    """
    multiqc . 

    """
}
// Process 6: Indexing reference genome & alignment --tool: bwa

process ALIGNMENTS {
    publishDir path: "${params.outdir}/aligned", mode: 'copy'
    tag "Alignment:"

    input:
    path(ref_ch)
    tuple path(Read1), path(Read2)
	
    output:
    path "${sample_id}.sam", emit: aligned_sam

    script:

    sample_id = ( Read1 =~ /(.+)_R1.paired.fastq/ )[0][1]
    aligned = "${sample_id}.sam"

    template 'align.sh'

}

// Process 7: Merging sam files --tool: samtools

process MERGESAM {
	publishDir path : "${params.outdir}/mergedsam", mode: 'copy'
	tag "Merging sam files:"
	
	input:
	path mergedsam

	output:
	path "Tryp-T.congolense.sam", emit: ecolisam

	script:
        tryp = "Tryp-T.congolense.sam"

	"""
	samtools merge ${tryp} ${mergedsam}

	"""
}

// Process 8: Converting to bam format --tool: samtools

process CONVERT_TO_BAM {
      publishDir path : "${params.outdir}/conversion", mode: 'copy'
      tag "Conversion"

      input:
      path align_sam

      output:
      path "Tryp-T.congolese_aligned.bam", emit: aligned_bam

      script:
      aligned = "Tryp-T.congolese_aligned.bam"

      """
      samtools view -Sb ${align_sam} > ${aligned}
      """
}

// Process 9: Sorting bam file --tool: samtools 

process SORTING {
      publishDir path: "${params.outdir}/sorting", mode : 'copy'
      tag "Sorting"

      input:
      path align_bam

      output:
      path "Tryp-T.congolese_aligned_sorted.bam", emit: sorted_bam
      path "Tryp-T.congolese_aligned_sorted.bam.bai", emit: sorted_bai
	
      script:
      sort = "Tryp-T.congolese_aligned_sorted.bam"
      sortidx = "Tryp-T.congolese_aligned_sorted.bam.bai"
      """
      samtools sort -O bam -o ${sort} ${align_bam}
      samtools index ${sort} > ${sortidx}

      """
}

// Process 10: Removing Duplicates --tool: GATK & samtools

process REMOVE_DUPLICATES {
      publishDir path: "${params.outdir}/Dedups", mode:'copy'
      tag "Removing Duplicates"

      input:
      path sort_bam

      output:
      path "marked_dups.bam", emit: marked_dups
      path "marked_dups_metrics.txt", emit: met_dedups

      script:
      dedups = "marked_dups.bam"
      metdedups = "marked_dups_metrics.txt"

      """
      gatk MarkDuplicates --INPUT ${sort_bam} --OUTPUT ${dedups} \
      --METRICS_FILE ${metdedups} --REMOVE_DUPLICATES true

      samtools index ${dedups} 
      """
}

// Process 11: Creating sequence dictionary --tool: GATK
process CREATE_SEQ_DICTIONARY {
	publishDir path: "${params.outdir}"
        tag " Creating Sequence Dictionary"

        input:
	path ref_chann

        output:
        path "TriTrypDB-52_TcongolenseIL3000_Genome.dict", emit: tryp_dict
        path "TriTrypDB-52_TcongolenseIL3000_Genome.fasta.fai", emit: tryp_fai

        script:
        fai = "TriTrypDB-52_TcongolenseIL3000_Genome.fasta.fai"

        """
	gatk CreateSequenceDictionary -R ${ref_chann}

        samtools faidx ${ref_chann} > ${fai}
        """
}

// Process 12: Base Recalibration  --tool: GATK

process BASERECALIBRATE {
	publishDir path: "${params.outdir}/baserecalibration", mode:'copy'
	
	tag "Recalibrating bases:"

	input:
	path ref_ref
        path ecoli_dit
	path ecoli_fai
       	path removeddups
	path knownvar

	output:
	path "recalibration_table", emit: recal_tabel
	
	script:
	recal = "recalibration_table"

	"""
	gatk IndexFeatureFile -F ${knownvar}
	gatk BaseRecalibrator -R ${ref_ref} -I ${removeddups} --known-sites ${knownvar} -O ${recal}

	"""
}

// Process 13: ApplyingBQSR --tool: GATK

process BQSR {
	publishDir path: "${params.outdir}/BQSR", mode:'copy'

	tag "Applying BQSR:"

	input:
	path ref_ch
	path ecoli_dict
	path ecoli_fai
	path dupsremoved
	path recaltab

	output:
	path "recalibratedfile.bam", emit: newbam

	script:
	
	bam = "recalibratedfile.bam"

	"""

        gatk ApplyBQSR -R ${ref_ch} -I ${dupsremoved} --bqsr-recal-file ${recaltab} -O ${bam}

        """
}


// Process 14: Calling Variants --tool: GATK

process VARIANT_CALL {
     publishDir path: "${params.outdir}/variants", mode: 'copy'
     tag "Variant call"

     input:
     path recal_bam
     path ref_c
     path ecoli_dict
     path ecoli_fai

     output:
     path "variantsGHC.vcf", emit: variants_ghc 

     script:
     variantsG = "variantsGHC.vcf"

     """
     gatk HaplotypeCaller --native-pair-hmm-threads 8 -I ${recal_bam} -O ${variantsG} -R ${ref_c}	
     
     """
}

// Process 15: Filtering out low quality variants --tool: bcftools 

process VARIANT_FILTER {
      publishDir path: "${params.outdir}/filteredvariants", mode: 'copy'
      tag "filtering variants"

      input:
      path variants

      output:
      path "Tryp-T.congolense_filteredb_variants.vcf", emit: filteredvariants

      script:
      var_filt = "Tryp-T.congolense_filteredb_variants.vcf"
      
      """
      gatk VariantFiltration -V ${variants} -filter "QD < 6.0" --filter-name "QD6" -filter "MQ < 40.0" --filter-name "MQ40" -filter "FS > 50.0" --filter-name "FS50" -filter "SOR > 3.75" --filter-name "SOR3.75" -filter "MQRankSum < -5.0" --filter-name "MQRankSum5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O variants_tagged.vcf

      gatk SelectVariants -V variants_tagged.vcf -O ${var_filt} -select 'vc.isNotFiltered()'
   
      """

}

// Process 16: Variant decomposing and normalizing --tool: vt

process NORMALIZATION {
	publishDir path: "${params.outdir}/normalized", mode: 'copy'
	tag "Decomposing and normalizing variants"

	input:
	path filtvar
	path ref_gen
        path tryp_dict
        path tryp_fai

	output:
	path "Tryp-T.congolense_normalized_var.vcf", emit: norm_var

	script:
	norm = "Tryp-T.congolense_normalized_var.vcf"
	decomp = "Tryp_T.congolense_decomp_var.vcf"

	"""
	vt decompose -s -o ${decomp} ${filtvar}

	vt normalize -r ${ref_gen} -o ${norm} ${decomp}
	"""
}

// Process 17: Variant annotation --tool: snpEFF

process ANNOTATION {
      publishDir path: "${params.outdir}/annotated", mode: 'copy'
      tag "Variant annotation"

      input:
      path normvar

      output:
      path "Tryp-T.congolense_anno_variants.vcf", emit: annotatevar
      path "snpEff_genes.txt", emit: snpEff_genes
      path "snpEff_summary.html", emit: snpEff_html

      script:
      annota_var = "Tryp-T.congolense_anno_variants.vcf"

      """
      # Download respective gff3 file for Trypanosoma congolense IL3000
      # mkdir IL3000 in snpEff/data \
      # && rename the downloaded gff file to genes.gff \
      # save it to snpEff/data/IL3000
      # Edit the snpEff.config file as:
      #Database for Trypanosoma congolense IL3000
        #IL3000.genome : Trypanosoma congolense IL3000
        #IL3000.reference : https://tritrypdb.org/common/downloads/Current_Release/TcongolenseIL3000/fasta/data/TriTrypDB-51_TcongolenseIL3000_Genome.fasta

      # Build reference database

      java -jar /home/hebrewsimeon/Annotater/tools/snpEff/snpEff.jar build -gff3 -v IL3000

      java -jar /home/hebrewsimeon/Annotater/tools/snpEff/snpEff.jar eff IL3000 ${normvar}  > ${annota_var}

      """
}
 
