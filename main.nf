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
Computes and Plots general statistics for processed Single-Stranded-DNA-Sequencing (SSDS) data.
The input data must come from ssdsnextflowpipeline version 2.0 
(see https://gitlab.igh.cnrs.fr/pauline.auffret/ssdsnextflowpipeline) 

Pipeline overview:
// PROCESS 1    : MULTIBAMSUMMARY (DEEPTOOLS)
// PROCESS 2    : PLOTCORRELATION (DEEPTOOLS)
// PROCESS 3    : COMPUTEMATRIX (DEEPTOOLS)
// PROCESS 4    : PLOTHEATMAP (DEEPTOOLS)
// PROCESS 5    : GETFORWARDSTRAND (SAMTOOLS AND DEEPTOOLS)
// PROCESS 6    : GETREVERSESTRAND (SAMTOOLS AND DEEPTOOLS)
// PROCESS 7    : COMPUTEMATRIXFR (DEEPTOOLS)
// PROCESS 8    : PLOTHEATMAPFR (DEEPTOOLS)
// PROCESS 9    : ANNOTATEPEAKS (HOMER)
// PROCESS 10    : PLOT_INTERSECT (INTERVENE)

*/

// Construct help message (option --help)
def helpMessage() {
    log.info"""
=============================================================================
  SSDS post-process pipeline version 1.0 : Computes and Plots general 
  statistics for processed Single-Stranded-DNA-Sequencing (SSDS) data
=============================================================================
    Usage:

            nextflow run main.f -c conf/igh.config -params-file conf/mm10.config -profile conda [options]

   Tested with Nextflow v20.07.1
=============================================================================
Input data parameters:
	-params_file            FILE    PATH TO PARAMETERS JSON FILE (template and default : /work/${USER}/ssdsnextflowpipeline/conf/mm10.json)
	--name			STRING  ANALYSIS NAME (default : "SSDS_postprocess_pipeline")
	--sample_name		STRING  SAMPLE OR GROUP NAME (default : " 
	--publishdir_mode	STRING  MODE FOR EXPORTING PROCESS OUTPUT FILES TO OUTPUT DIRECTORY (default : "copy", must be "symlink", "rellink", "link", "copy", "copyNoFollow","move", see https://www.nextflow.io/docs/latest/process.html)
    	--genomebase	        DIR     PATH TO REFERENCE GENOMES (default : "/poolzfs/genomes")
    	--genome		STRING  REFERENCE GENOME NAME (must correspond to an existing genome in your config file, default : "mm10")
	--finalpeaksbed		FILE	FILTERED PEAKS IN BED FORMAT (for example, coming from ssdsnextflowprocess, in the finalpeaks folder (warning : if finalpeakbed is not None, bamreference cannot be None) ; default : )
	--peakreference		FILE	REFERENCE PEAKS FILE IN BED FORMAT (set to "None" if none provided (warning : if finalpeakbed is not None, bamreference cannot be None) : set of peaks to compare the finalpeaksbed to (for example, hotspots from B6 WT mouse) ; default : )
	--bamfolder		DIR	ABSOLUTE PATH TO BAM FOLDER CONTAINING FILTERED TYPE 1 BAM FILES FROM SSDSNEXTFLOWPIPELINE (for example, coming from ssdsnextflowprocess, in the bwa/filterbam/flag_*/parse_itr/type1/bam folder ; default :
	--bampattern		REGEX	PATTERN FOR MATCHING BAM FILES IN BAMFOLDER (default : "*.bam")
	--bamreference		DIR     ABSOLUTE PATH TO REFERENCE BAM FOLDER CONTAINING FILTERED TYPE 1 BAM FILES FOR REFERENCE (set to "None" if none provided (warning : if bamreference is not None, peakreference cannot be None) 
	--bigwigfolder		DIR     ABSOLUTE PATH TO BIGWIG FOLDER CONTAINING BIGWIG FILES FROM SSDSNEXTFLOWPIPELINE (for example, coming from ssdsnextflowprocess, in the bigwig/kbrick_bigwig/deeptools/binsize* folder ; default : 
	--bigwigpattern		REGEX   PATTERN FOR MATCHING BIGWIG FILES IN BIGWIGFOLDER (default : "*ssDNA_type1.deeptools.RPKM.bigwig")
	--bedfolder		DIR     ABSOLUTE PATH TO BED FOLDER CONTAINING BED TO COMPUTE INTERSECT (default : 
	--bedpattern		REGEX   PATTERN FOR MATCHING BED FILES IN BED FOLDER (default : "*.bed")
	--corMethod		STRING	CORRELATION METHOD FOR DEEPTOOLS PLOTCORRELATION PROCESS (default : spearman; valid options are spearman, pearson)
	--corPlot		STRING	DEEPTOOLS whatToPlot OPTION FOR PLOTCORRELATION PROCESS (default : heatmap; valid options are heatmap, scatterplot)
	--heatmap_width		INT	DEEPTOOLS heatmapWidth OPTION FOR PLOTHEATMAP PROCESS (default : 20)
	--matrix_downstream	INT	DEEPTOOLS downstream OPTION FOR COMPUTEMATRIX PROCESS (default : 2500)
	--matrix_upstream	INT     DEEPTOOLS upstream OPTION FOR COMPUTEMATRIX PROCESS (default : 2500)
	--gtf_id		STRING	GTF FEATURE TYPE FOR HOMER (default : "-gid" for gene_id ; valid options are "-gid" and "" for transcript_id which is HOMER default)
	--figType		STRING	FILE FORMAT FOR THE PLOT (default : pdf ; valid options are pdf,svg,ps,tiff,png)
	--intersect_options	STRING	INTERVENE --bedtools option (default : "-f 1E-9" ; see bedtools intersect documentation https://bedtools.readthedocs.io/en/latest/)
	--intersect_threshold	STRING	INTERVENE --intersect_thresh option (default : 1 ; see bedtools intersect documentation https://bedtools.readthedocs.io/en/latest/)
	--with_clustering	BOOL	If false then clustering step is skipped (default : true)
	--with_heatmap		BOOL	If false then heatmap step is skipped (default : true)
	--with_FR_heatmap	BOOL	If false then heatmap for reverse and forwad step is skipped (default : true)
	--with_FR_bigwig	BOOL	If false then bigwigs for reverse and forwad step is skipped (default : true)
	--with_peak_annot	BOOL	If false then peak annotation step is skipped (default : true)
	--with_plot_intersect	BOOL	If false then bed intersection step is skipped (default : true)

=============================================================================

=============================================================================
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

// Custom name variables
def runName = "${params.name}.$workflow.runName"

// External scripts used in the pipeline
def plot_homer_annotation_script = "${params.src}/plot_homer_annotatepeaks.r"  //Adapted from nf-core chipseq pipeline version 1.2.1

// Check if input files/directories exist
if (params.finalpeaksbed) { println("Checking finalpeaksbed input file...") ; f_finalpeaks = file(params.finalpeaksbed, checkIfExists: true) } ; println("Ok")
if (params.peakreference != "None") { println("Checking peakreference input file...") ; f_peaksref = file(params.peakreference, checkIfExists: true) } ; println("Ok")
if (params.bamfolder) { println("Checking bamfolder directory...") ; d_bam = file(params.bamfolder, checkIfExists: true) } ; println("Ok")
if (params.bigwigfolder) { println("Checking bigwigfolder  directory...") ; d_bw = file(params.bigwigfolder, checkIfExists: true) } ; println("Ok")
if (params.bedfolder) { println("Checking bedfolder  directory...") ; d_bw = file(params.bedfolder, checkIfExists: true) } ; println("Ok")
if (params.bamreference != "None") { println("Checking bamreference directory...") ; d_bamref = file(params.bamreference, checkIfExists: true) } ; println("Ok")
if (params.bamreference != "None" && params.peakreference == "None" ) { println("Checking conformity of reference parameters peakreference and bamreference") ; println("Error : bamreference parameter is set but there is no peakreference parameter ; either set bamreference to None or set up a peakreference file") ; exit 1 } ; else println("Ok")
if (params.bamreference == "None" && params.peakreference != "None" ) { println("Checking conformity of reference parameters peakreference and bamreference") ; println("Error : peakreference parameter is set but there is no bamreference parameter ; either set peakreference to None or set up a bamreference folder") ; exit 1 } ; else println("Ok")

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genome.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Define genome variables
params.genome_fasta = params.genome ? params.genomes[ params.genome ].genome_fasta ?: false : false
params.genomedir = params.genome ? params.genomes[ params.genome ].genomedir ?: false : false
params.genome_name = params.genome ? params.genomes[ params.genome ].genome_name ?: false : false
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
params.genome_gtf = params.genome ? params.genomes[ params.genome ].genome_gtf ?: false : false

// Check input parameters conformity
if(params.publishdir_mode!="copy" && params.publishdir_mode!="symlink" && params.publishdir_mode!="rellink" && params.publishdir_mode!="link" && params.publishdir_mode!="copyNoFollow" && params.publishdir_mode!="move") {
    println("Error : --publishdir_mode must be symlink, rellink, link, copy, copyNoFollow,move, see https://www.nextflow.io/docs/latest/process.html")
    exit 0
}


//***************************************************************************//
//                                                                           //
//                          BEGINNING PIPELINE                               //
//              
// **************************************************************************//

//***************************************************************************//
//                         SECTION 1 : CLUSTERING                            //
//***************************************************************************//
// PROCESS 1    : MULTIBAMSUMMARY (DEEPTOOLS) 
// What it does : computes the read coverages of filtered type1 BAM files at hotspots
// Input        : bed files containing coordinates of hotspots and filtered type1 BAM files from ssdsnextflowpipeline
// and optionaly, the reference hospots bed file (or any bed file containing peaks for comparison purpose)
// Output       : compressed matrix of values as a numpy array (.npz file) 
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html
process multiBamSummary {
    tag "${params.sample_name}"
    label 'process_basic'
    conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/clustering/matrices",   mode: params.publishdir_mode, pattern: "*.npz*"
    publishDir "${params.outdir}/clustering/log",        mode: params.publishdir_mode, pattern: "*.log*"
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
    # If another bed file is given to compare with hotspots 
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
        # If no other bed file is given, creates an empty file just so all the outputs expected by the process exist in the end of the process
        touch ${params.sample_name}_peakref.npz
    fi  
    """
}

// PROCESS 2    : PLOTCORRELATION (DEEPTOOLS) 
// What it does : plots an heatmap and a PCA of correlation of samples based on multibamsummary matrix
// Input        : compressed matrix of values as a numpy array (.npz file) from PROCESS 1 output
// Output       : heatmap and PCA of samples in pdf format 
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html?highlight=plotcorrelation
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
    plotPCA --corData $finalpeaks_matrix \
            --plotTitle ${params.sample_name} \
            --transpose \
            -o ${params.sample_name}.finalpeaks.plotPCA.${params.corMethod}.${params.corPlot}.png \
            >& ${params.sample_name}.finalpeaks.plotPCA.${params.corMethod}.${params.corPlot}.log 2>&1 
    # If another bed file is given to compare with hotspots
    if [ ${params.peakreference} != "None" ]
    then
        plotCorrelation --corData $peakref_matrix \
                        --corMethod ${params.corMethod} \
                        --whatToPlot ${params.corPlot} \
                        -o ${params.sample_name}.peakref.plotCorrelation.${params.corMethod}.${params.corPlot}.png \
                        >& ${params.sample_name}.peakref.plotCorrelation.${params.corMethod}.${params.corPlot}.log 2>&1
        plotPCA --corData $peakref_matrix \
                --plotTitle ${params.sample_name} \
                --transpose \
                -o ${params.sample_name}.peakref.plotPCA.${params.corMethod}.${params.corPlot}.png \
                >& ${params.sample_name}.peakref.plotPCA.${params.corMethod}.${params.corPlot}.log 2>&1
    fi
    """
}

// PROCESS 3    : COMPUTEMATRIX (DEEPTOOLS) 
// What it does : calculates scores per genome regions
// Input        : bed files containing coordinates of hotspots and filtered type1 BAM files and bigwig files from ssdsnextflowpipeline
// Output       : compressed matrix of values (.gzip  file)
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html?highlight=computematrix
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

// PROCESS 4    : PLOTHEATMAP (DEEPTOOLS)
// What it does : creates a heatmap for scores associated with hotspots
// Input        : Matrix file from the computeMatrix process
// Output       : heatmaps in pdf format
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html?highlight=plotheatmap
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
                -o ${params.sample_name}_heatmap.png \
                >& ${params.sample_name}_heatmap.log 2>&1
    """
}                

// Here I create a channel grouping all bam files to an identifier
// The resulting channel is composed of 2 elements : [Id; path-to-bam]
Channel
    .fromPath( "${params.bamfolder}/${params.bampattern}" ) 
    .map { tuple( it.name.split('_')[0..-3].join('_'), it ) }
    //.println() 
    .into{ bamF_ch ; bamR_ch }  

// PROCESS 5    : GETFORWARDSTRAND (SAMTOOLS AND DEEPTOOLS)
// What it does : Gets coverage originated from the forward strand from a bam file
// Input        : Filtered bam file of type1 ssds reads (from ssdsnextflowpipeline)
// Output       : forward coverage in bigwig format
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage
process getForwardStrand {
    tag "${bam_id}"
    label 'process_basic'
    conda 'deeptools=3.5.1 bioconda::samtools=1.14'
    publishDir "${params.outdir}/FRbigwig",     mode: params.publishdir_mode, pattern: "*.bigwig"
    publishDir "${params.outdir}/FRbigwig/log", mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(bam_id), path(bam) from bamF_ch
    output:
        tuple val(bam_id), path('*bigwig') into fwd_bigwig_ch
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

// PROCESS 6    : GETREVERSESTRAND (SAMTOOLS AND DEEPTOOLS)
// What it does : Gets coverage originated from the reverse strand from a bam file
// Input        : Filtered bam file of type1 ssds reads (from ssdsnextflowpipeline)
// Output       : reverse coverage in bigwig format
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage
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

// Here I create a channel linking forward and reverse bigwig files
rev_bigwig_ch
    .join(fwd_bigwig_ch)
    .set { FR_bigwig_ch }
    //.println()

// PROCESS 7    : COMPUTEMATRIXFR (DEEPTOOLS) 
// What it does : calculates scores per genome regions
// Input        : bed files containing coordinates of hotspots and 2 bigwig files per samples
// one containing forward signal (from process 5), the other containing reverse signal (from process 6)
// Output       : compressed matrix of values (.gzip  file) 
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html?highlight=computematrix
process computeMatrixFR {
    tag "${bam_id}"
    label 'process_basic'
    conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/heatmap/matrices", mode: params.publishdir_mode, pattern: "*.matrix.FR"
    publishDir "${params.outdir}/heatmap/log",      mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(bam_id), path(rev_bw), path(fwd_bw) from FR_bigwig_ch
    output:
        tuple val(bam_id), path('*finalpeaks.matrix.FR') into matrix_FR_ch
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

// PROCESS 8    : PLOTHEATMAPFR (DEEPTOOLS)
// What it does : creates a heatmap for scores associated with hotspots
// Input        : Matrix file from the computeMatrixFR process
// Output       : heatmaps in pdf format showing forward and reverse signal around hotspots
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html?highlight=plotheatmap
process plotHeatmapFR {
    tag "${bam_id}"
    label 'process_basic'
    conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/heatmap/plots", mode: params.publishdir_mode, pattern: "*.png"
    publishDir "${params.outdir}/heatmap/log",   mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(bam_id), path(matrix) from matrix_FR_ch
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

// PROCESS 9    : ANNOTATEPEAKS (HOMER AND R)
// What it does : annotate peaks (hotspots) using Homer
// Input        : bed file containing peaks
// Output       : histograms and violin plots
// External tool: plot_homer_annotatepeaks.r script from nf-core chipseq pipeline version 1.2.1 (Ewels et al. 2020) 
// Resource     : http://homer.ucsd.edu/homer/index.html
process annotatePeaks {
    tag "${params.sample_name}"
    label 'process_basic'
    conda 'bioconda::homer=4.11 r::r=3.6.0 conda-forge::r-optparse=1.6.4 conda-forge::r-ggplot2=3.3.0 conda-forge::r-reshape2=1.4.3'
    publishDir "${params.outdir}/annotations/plots", mode: params.publishdir_mode, pattern: "*.pdf"
    publishDir "${params.outdir}/annotations",       mode: params.publishdir_mode, pattern: "*.txt"
    publishDir "${params.outdir}/annotations/log",   mode: params.publishdir_mode, pattern: "*.log"
    output:
        path('*.txt')
        path('*.pdf')
        path('*.log')
    when:
        params.with_peak_annot
    script:
    """
    ## Part 1 : Annotates with Homer
    # If working with mouse data mapped on mm10 genome
    if [ ${params.genome} == "mm10" ]
    then
        #download UCSC mm10 annotation files for Homer
        perl \$CONDA_PREFIX/share/homer/.//configureHomer.pl -install mm10
        annotatePeaks.pl ${params.finalpeaksbed} \
                     mm10 \
                     1> ${params.sample_name}_annotatePeaks_homer_${params.genome}.txt \
                     2> ${params.sample_name}_annotatePeaks_homer_${params.genome}.log
    else
        annotatePeaks.pl ${params.finalpeaksbed} \
                     ${params.genome_fasta} \
                     -${params.gtf_id} \
                     -gtf ${params.genome_gtf} \
                     1> ${params.sample_name}_annotatePeaks_homer_${params.genome}.txt \
                     2> ${params.sample_name}_annotatePeaks_homer_${params.genome}.log
    fi
    ## Part 2 : Plots with R
    ${plot_homer_annotation_script} -i ${params.sample_name}_annotatePeaks_homer_${params.genome}.txt \
                                    -s ${params.sample_name} \
                                    -o ./ \
                                    -p ${params.sample_name}.plothomer \
                                    >& ${params.sample_name}.plothomer.log 2>&1
    """
}

// PROCESS 10    : PLOT_INTERSECT (INTERVENE)
// What it does : Compute overlap of bed sets using intervene
// Input        : list of bed files
// Output       : overlapping bed files and upset plot of overlaps
// Resource     : https://intervene.readthedocs.io/en/latest/modules.html
process plotIntersect {
    tag "${params.sample_name}"
    label 'process_basic'
    conda 'bioconda::intervene=0.6.4'
    publishDir "${params.outdir}/intersect/plots",   mode: params.publishdir_mode, pattern: "*.${params.figType}"
    publishDir "${params.outdir}/intersect/overlap", mode: params.publishdir_mode, pattern: "*.bed"
    publishDir "${params.outdir}/intersect/log",     mode: params.publishdir_mode, pattern: "*.log"
    output:
        path('*.bed')
        path('*.${params.figType}')
        path('*.log')
    when:
        params.with_plot_intersect
    script:
    """
    #Run intervene to compute the intesecting sets
    intervene upset -i ${params.bedfolder}/${params.bedpattern} \
                    --scriptonly --figtype ${params.figType} \
                    --bedtools-options ${params.intersect_options} \
                    --overlap-thresh ${params.intersect_threshold} \
                    --project ${params.sample_name} \
                    --save-overlaps >& ${params.sample_name}_intervene.log 2>&1

    #There is a bug (need to investigate but no time now) : need to edit the R script generated by intervene before execution otherwise the execution is halted
    sed -i "s/, mainbar.y.label =\"No. of Intersections\", sets.x.label =\"Set size\")/)/g" ${params.sample_name}_upset.R >& ${params.sample_name}_intervene_sed.log 2>&1

    #Finally, esecute Rscript to plot the upset plot
    Rscript ${params.sample_name}_upset.R >& ${params.sample_name}_intervene_rscript.log 2>&1
    """
}

