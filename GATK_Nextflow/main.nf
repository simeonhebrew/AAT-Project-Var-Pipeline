// Enable DSL 2 syntax

nextflow.enable.dsl = 2

// import modules here

include { QUALITY_CHECK; RAWMULTIQC; TRIMMING; TRIMFASTQC; TRIMMULTIQC; ALIGNMENTS; MERGESAM; CONVERT_TO_BAM; SORTING; REMOVE_DUPLICATES; CREATE_SEQ_DICTIONARY; BASERECALIBRATE; BQSR; VARIANT_CALL; VARIANT_FILTER; NORMALIZATION; ANNOTATION } from "./modules/gatkHC.nf"

// set input channels
Channel.fromFilePairs( params.reads, checkIfExists:true)
	.set { read_pairs_ch }
	
Channel.fromPath ( params.reference, checkIfExists:true )
        .set { reference_ch }

Channel.fromPath ( params.knownsites, checkIfExists:true )
	.set { knownsites_ch }

Channel.fromPath ( params.adapters, checkIfExists:true )
	.set { adapters_ch }

Channel.fromPath ( params.slidewindow, checkIfExists:true )
	.set { slidewindow_ch }

Channel.fromPath ( params.qualdepthvalue, checkIfExists:true )
	.set { qual_ch }
 

// Run the workflow
workflow {
// Process 1 Quality Checking
QUALITY_CHECK(read_pairs_ch)

// Process 2 MultiQC for raw reads
RAWMULTIQC(QUALITY_CHECK.out.collect())

// Process 3 Trimming
TRIMMING(read_pairs_ch, slidewindow_ch.collect(), adapters_ch.collect())

// Process 4 Quality Checking trimmed reads
TRIMFASTQC(TRIMMING.out)

// Process 5 MultiQC for trimmed reads
TRIMMULTIQC(TRIMFASTQC.out.collect())

// Process 6 Alignment
ALIGNMENTS(reference_ch.collect(), TRIMMING.out)

// Process 7 Merge Sam files
MERGESAM(ALIGNMENTS.out.collect())

// Process 8 Bam Conversion
CONVERT_TO_BAM(MERGESAM.out)

// Process 9 Sorting
SORTING(CONVERT_TO_BAM.out)

// Process 10 Remove duplicates
REMOVE_DUPLICATES(SORTING.out.sorted_bam)

// Process 11 Creating sequence dictionary
CREATE_SEQ_DICTIONARY(reference_ch)

// Process 12 Baserecalibration
BASERECALIBRATE(reference_ch, CREATE_SEQ_DICTIONARY.out, REMOVE_DUPLICATES.out.marked_dups, knownsites_ch)

// Process 13 BQSR
BQSR(reference_ch, CREATE_SEQ_DICTIONARY.out, REMOVE_DUPLICATES.out.marked_dups, BASERECALIBRATE.out)

// Process 14 Variant calling
VARIANT_CALL(BQSR.out, reference_ch, CREATE_SEQ_DICTIONARY.out)

// Process 15 Variant filter
VARIANT_FILTER(VARIANT_CALL.out)

// Process 16 Normalization
NORMALIZATION(VARIANT_FILTER.out, reference_ch, CREATE_SEQ_DICTIONARY.out)

// Process 17 Annotation
ANNOTATION(NORMALIZATION.out)

}
