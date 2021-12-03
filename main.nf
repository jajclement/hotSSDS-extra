#!/usr/bin/env nextflow
/*
========================================================================================
                        SSDS post-process pipeline version 1.0
                        Pauline Auffret, 2021
                        Contact : pauline.auffret@igh.cnrs.fr
========================================================================================
 SSDS post-process pipeline
 #### Homepage / Documentation
 https://gitlab.igh.cnrs.fr/pauline.auffret/ssdspostprocess
 Contact : pauline.auffret@igh.cnrs.fr
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Compute and Plot general statistics for processed Single-Stranded-DNA-Sequencing (SSDS) data.
The input data MUST come from ssdsnextflowpipeline version 2.0 
(see https://gitlab.igh.cnrs.fr/pauline.auffret/ssdsnextflowpipeline) 
or MUST be organized according to the following classification :

Pipeline overview:


*/

// Construct help message (option --help)
def helpMessage() {
    log.info"""
=============================================================================
  SSDS Pipeline version 2.0 : Align, parse and call hotspots from SSDNA
=============================================================================
    Usage:

    Runs with Nextflow v20.04.1
=============================================================================
Input data parameters:
    --inputcsv                  FILE    PATH TO INPUT CSV FILE (template and default : /work/${USER}/ssdsnextflowpipeline/tests/fastq/input.csv)
    -params_file                FILE    PATH TO PARAMETERS JSON FILE (template and default : /work/${USER}/ssdsnextflowpipeline/conf/mm10.json)

=================================================================================================
""".stripIndent()
}

//***************************************************************************//
//           PRELIMINARY SECTION : GLOBAL VARIABLES ANS SETTINGS             //
//***************************************************************************//
//Show pipeline version
params.version = false
if (params.version){
    println("This is $workflow.manifest.name version $workflow.manifest.version")
    println("$workflow.manifest.description")
    exit 0
}

//Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Define global variables
// Create scratch directory
scrdir = file("${params.scratch}")
result = scrdir.mkdirs()
println result ? "Create scratch directory...OK" : "Cannot create directory: $scrdir"

// Custom name variables
def runName = "${params.name}.$workflow.runName"



//***************************************************************************//
//                                                                           //
//                          BEGINNING PIPELINE                               //
//              
// **************************************************************************//

//***************************************************************************//
//                     SECTION 1 : INPUT SETTINGS                            //
//***************************************************************************//
//Channel
//   .from( params.inputdir.tokenize() )                          // split inputdir parameter in separate file names
//   .map { it -> [ it.split('/')[-1], it ] }                     // get dir name
//   //.map { it -> [ it[0], files(it[1], checkIfExists: true) ] }  // converts the file path string to a file objects
//   .map { it -> it[0,1].flatten() }
//   .set{ ssdsdir_ch }  

//ssdsdir_ch.view { "value: $it" }

process multiBamSummary {
    tag "${params.sample_name}"
    label 'process_basic'
    conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/clustering/matrices",   mode: params.publishdir_mode, pattern: "*.npz*"
    publishDir "${params.outdir}/clustering/log",        mode: params.publishdir_mode, pattern: "*.log*"
//    input:
//        tuple val(ssds_name), val(ssds_dir) from ssdsdir_ch
    output:
        tuple path('*finalpeaks.npz'), path('*peakref.npz') into bamSummary_ch 
        path('*.log')
    when:
        params.with_clustering
    script:
    """
    multiBamSummary BED-file --BED ${params.finalpeaksbed} \
                    --bamfiles ${params.bamfolder}/*.bam \
                    --outRawCounts ${params.sample_name}_rawCounts.txt \
                    --scalingFactors ${params.sample_name}_scalingFactors.txt \
                    --smartLabels \
                    --numberOfProcessors ${task.cpus} \
                    -o ${params.sample_name}_finalpeaks.npz >& ${params.sample_name}_finalpeaks.log 2>&1
    if [ ${params.peakreference} != "None" ]
    then
        multiBamSummary BED-file --BED ${params.peakreference} \
                        --bamfiles ${params.bamfolder}/*.bam ${params.bamreference}/*.bam \
                        --outRawCounts ${params.sample_name}_rawCounts.txt \
                        --scalingFactors ${params.sample_name}_scalingFactors.txt \
                        --smartLabels \
                        --numberOfProcessors ${task.cpus} \
                        -o ${params.sample_name}_peakref.npz >& ${params.sample_name}_peakref.npz.log 2>&1
    else
        touch ${params.sample_name}_peakref.npz
    fi  
    """
}

process plotCorrelation {
    tag "${params.sample_name}"
    label 'process_basic'
    conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/clustering/plots", mode: params.publishdir_mode, pattern: "*.png*"
    publishDir "${params.outdir}/clustering/log",   mode: params.publishdir_mode, pattern: "*.log*"
    input:
        tuple path(finalpeaks_matrix), path(peakref_matrix) from bamSummary_ch
    output:
        path('*.png')
        path('*.log')
    when:
        params.with_clustering
    script:
    """
    plotCorrelation --corData $finalpeaks_matrix \
                    --corMethod ${params.corMethod} \
                    --whatToPlot ${params.corPlot} \
                    -o ${params.sample_name}.finalpeaks.plotCorrelation.${params.corMethod}.${params.corPlot}.png \
                    >& ${params.sample_name}.finalpeaks.plotCorrelation.${params.corMethod}.${params.corPlot}.log 2>&1
    if [ ${params.peakreference} != "None" ]
    then
        plotCorrelation --corData $peakref_matrix \
                        --corMethod ${params.corMethod} \
                        --whatToPlot ${params.corPlot} \
                        -o ${params.sample_name}.peakref.plotCorrelation.${params.corMethod}.${params.corPlot}.png \
                        >& ${params.sample_name}.peakref.plotCorrelation.${params.corMethod}.${params.corPlot}.log 2>&1
    fi
    """
}

process computeMatrix {
    tag "${params.sample_name}"
    label 'process_basic'
    conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/heatmap/matrices", mode: params.publishdir_mode, pattern: "*.matrix"
    publishDir "${params.outdir}/heatmap/log",      mode: params.publishdir_mode, pattern: "*.log"
    output:
        path('*finalpeaks.matrix') into matrix_ch
        path('*.log')
    when:
        params.with_heatmap
    script:
    """
    computeMatrix   reference-point --regionsFileName ${params.finalpeaksbed} \
                    --scoreFileName ${params.bigwigfolder}/${params.bigwigpattern} \
                    --referencePoint=center \
                    --downstream=${params.matrix_downstream} \
                    --upstream=${params.matrix_upstream} \
                    --smartLabels \
                    --numberOfProcessors ${task.cpus} \
                    -o ${params.sample_name}_b${params.matrix_upstream}_a${params.matrix_downstream}.finalpeaks.matrix \
                    >& ${params.sample_name}_finalpeaks.matrix.log 2>&1
    """
}

process plotHeatmap {
    tag "${params.sample_name}"
    label 'process_basic'
    conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/heatmap/plots", mode: params.publishdir_mode, pattern: "*.png"
    publishDir "${params.outdir}/heatmap/log",   mode: params.publishdir_mode, pattern: "*.log"
    input:
        path(matrix) from matrix_ch
    output:
        path('*.png')
        path('*.log')
    when:
        params.with_heatmap
    script:
    """
    plotHeatmap --matrixFile $matrix \
                --refPointLabel=center \
                --legendLocation=best \
                --perGroup \
                --regionsLabel hotspots \
                --heatmapWidth=${params.heatmap_width} \
                --colorMap coolwarm \
                --plotTitle ${params.sample_name} \
                -o ${params.sample_name}_heatmap_files_${params.bigwigpattern}.png \
                >& ${params.sample_name}_heatmap.log 2>&1
    """
}                

Channel
    .fromPath( "${params.bamfolder}/${params.bampattern}" ) 
    .map { tuple( it.name, it ) }
    //.println() 
    .into{ bamF_ch ; bamR_ch }  

process getForwardStrand {
    tag "${bam_id}"
    label 'process_basic'
    conda 'deeptools=3.5.1 bioconda::samtools=1.14'
    publishDir "${params.outdir}/FRbigwig",     mode: params.publishdir_mode, pattern: "*.bigwig"
    publishDir "${params.outdir}/FRbigwig/log", mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(bam_id), path(bam) from bamF_ch
    output:
        tuple val(bam_id) ,path('*bigwig') into fwd_bigwig_ch
        path('*.log')
    when:
        params.with_FR_bigwig
    script:
    """
    # include reads that are 2nd in a pair (128);
    # exclude reads that are mapped to the reverse strand (16)
    samtools view -b -f 128 -F 16 $bam > ${bam.baseName}.fwd1.bam

    # exclude reads that are mapped to the reverse strand (16) and
    # first in a pair (64): 64 + 16 = 80
    samtools view -b -f 80 $bam > ${bam.baseName}.fwd2.bam

    # combine the temporary files
    samtools merge -f ${bam.baseName}.fwd.bam ${bam.baseName}.fwd1.bam ${bam.baseName}.fwd2.bam

    # index the filtered BAM file
    samtools index ${bam.baseName}.fwd.bam

    # run bamCoverage
    bamCoverage -b ${bam.baseName}.fwd.bam -o ${bam.baseName}.fwd.bigwig >& ${params.sample_name}_${bam.baseName}_bamCoverage_fwd.log

    # remove the temporary files
    rm ${bam.baseName}.fwd*.bam*
    """
}

process getReverseStrand {
    tag "${bam_id}"
    label 'process_basic'
    conda 'deeptools=3.5.1 bioconda::samtools=1.14'
    publishDir "${params.outdir}/FRbigwig",     mode: params.publishdir_mode, pattern: "*.bigwig"
    publishDir "${params.outdir}/FRbigwig/log", mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(bam_id), path(bam) from bamR_ch
    output:
        tuple val(bam_id), path('*bigwig') into rev_bigwig_ch
        path('*.log')
    when:
        params.with_FR_bigwig
    script:
    """
    # include reads that map to the reverse strand (128)
    # and are second in a pair (16): 128 + 16 = 144
    samtools view -b -f 144 $bam > ${bam.baseName}.rev1.bam

    # include reads that are first in a pair (64), but
    # exclude those ones that map to the reverse strand (16)
    samtools view -b -f 64 -F 16 $bam > ${bam.baseName}.rev2.bam

    # merge the temporary files
    samtools merge -f ${bam.baseName}.rev.bam ${bam.baseName}.rev1.bam ${bam.baseName}.rev2.bam

    # index the merged, filtered BAM file
    samtools index ${bam.baseName}.rev.bam

    # run bamCoverage
    bamCoverage -b ${bam.baseName}.rev.bam -o ${bam.baseName}.rev.bigwig >& ${params.sample_name}_${bam.baseName}_bamCoverage_rev.log

    # remove temporary files
    rm ${bam.baseName}.rev*.bam*
    """
}


rev_bigwig_ch
    .join(fwd_bigwig_ch)
    .set { FR_bigwig_ch }
    //.println()


process computeMatrixFR {
    tag "${bam_id}"
    label 'process_basic'
    conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/heatmap/matrices", mode: params.publishdir_mode, pattern: "*.matrix.FR"
    publishDir "${params.outdir}/heatmap/log",      mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(bam_id), path(rev_bw), path(fwd_bw) from FR_bigwig_ch
    output:
        path('*finalpeaks.matrix.FR') into matrix_FR_ch
        path('*.log')
    when:
        params.with_FR_heatmap
    script:
    """
    computeMatrix   reference-point --regionsFileName ${params.finalpeaksbed} \
                    --scoreFileName ${fwd_bw} ${rev_bw} \
                    --referencePoint=center \
                    --downstream=${params.matrix_downstream} \
                    --upstream=${params.matrix_upstream} \
                    --smartLabels \
                    --numberOfProcessors ${task.cpus} \
                    -o ${bam_id}_b${params.matrix_upstream}_a${params.matrix_downstream}.finalpeaks.matrix.FR \
                    >& ${bam_id}_finalpeaks.matrix.FR.log 2>&1
    """
}

process plotHeatmapFR {
    tag "${bam_id}"
    label 'process_basic'
    conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/heatmap/plots", mode: params.publishdir_mode, pattern: "*.png"
    publishDir "${params.outdir}/heatmap/log",   mode: params.publishdir_mode, pattern: "*.log"
    input:
        path(matrix) from matrix_FR_ch
    output:
        path('*.png')
        path('*.log')
    when:
        params.with_FR_heatmap
    script:
    """
    plotHeatmap --matrixFile $matrix \
                --refPointLabel=center \
                --legendLocation=best \
                --perGroup \
                --regionsLabel hotspots \
                --heatmapWidth=${params.heatmap_width} \
                --colorMap coolwarm \
                --plotTitle ${bam_id} \
                -o ${bam_id}_heatmap.FR.png \
                >& ${bam_id}_heatmap.FR.log 2>&1
    """
}                








