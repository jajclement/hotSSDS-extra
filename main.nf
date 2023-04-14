#!/usr/bin/env nextflow
/*
========================================================================================
                        SSDS post-process pipeline version 1.0
                        Pauline Auffret, 2021
                        Contact : Pauline.Auffret@ifremer.fr
			Secondary contact : julie.clement@univ-perp.fr
========================================================================================
 SSDS post-process pipeline
 #### Homepage / Documentation
 https://gitlab.igh.cnrs.fr/pauline.auffret/ssdspostprocess
 Contact : Pauline.Auffret@ifremer.fr
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
// PROCESS 10   : PLOT_INTERSECT (INTERVENE)
// PROCESS 11   : GENERAL REPORT (MULTIQC)
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
        --sample_name           STRING  SAMPLE OR GROUP NAME (default : "DMC1-ChIP")       
	--finalpeaksbed		FILE	FILTERED PEAKS IN BED FORMAT (for example, coming from ssdsnextflowprocess, in the finalpeaks folder (warning : if finalpeakbed is not None, bamreference cannot be None) ; default : )
	--peakreference		FILE	REFERENCE PEAKS FILE IN BED FORMAT (set to "None" if none provided (warning : if finalpeakbed is not None, bamreference cannot be None) : set of peaks to compare the finalpeaksbed to (for example, hotspots from B6 WT mouse) ; default : )
	--bamfolder		DIR	ABSOLUTE PATH TO BAM FOLDER CONTAINING FILTERED TYPE 1 BAM FILES FROM SSDSNEXTFLOWPIPELINE (for example, coming from ssdsnextflowprocess, in the bwa/filterbam/flag_*/parse_itr/type1/bam folder ; default :
	--bampattern		REGEX	PATTERN FOR MATCHING BAM FILES IN BAMFOLDER (default : "*.bam")
	--bamreference		DIR     ABSOLUTE PATH TO REFERENCE BAM FOLDER CONTAINING FILTERED TYPE 1 BAM FILES FOR REFERENCE (set to "None" if none provided (warning : if bamreference is not None, peakreference cannot be None) 
	--bigwigfolder		DIR     ABSOLUTE PATH TO BIGWIG FOLDER CONTAINING BIGWIG FILES FROM SSDSNEXTFLOWPIPELINE (for example, coming from ssdsnextflowprocess, in the bigwig/kbrick_bigwig/deeptools/binsize* folder ; default : 
	--bigwigpattern		REGEX   PATTERN FOR MATCHING BIGWIG FILES IN BIGWIGFOLDER (default : "*ssDNA_type1.deeptools.RPKM.bigwig")
	--bedfolder		DIR     ABSOLUTE PATH TO BED FOLDER CONTAINING BED TO COMPUTE INTERSECT (default : ) 
	--bedpattern		REGEX   PATTERN FOR MATCHING BED FILES IN BED FOLDER (default : "*.bed")
	--fragfolder		DIR	ABSOLUTE PATH TO THE FOLDER CONTAINING BED FILE WITH FRAGMENTS AFTER ITR (default: )
	--fragpattern		REGEX	PATTERN FOR MATCHING FRAGMENT BED FILES IN FRAGFOLDER (default: "*.bed")

Genome parameters:	
        --genomebase            DIR     PATH TO REFERENCE GENOMES (default : "/poolzfs/genomes")
        --genome                STRING  REFERENCE GENOME NAME (must correspond to an existing genome in your config file, default : "mm10")
        --genomedir             DIR     PATH TO GENOME DIRECTORY (required if your reference genome is not present in your config file)
        --genome_fasta          FILE    PATH TO FILE GENOME FASTA FILE WITH PREEXISTING INDEX FILES FOR BWA (required if your reference genome is not present in your config file)
        --genome_gtf
	--genome_size		FILE	PATH TO GENOME SIZE a two-column textfile containing chromosome name in the first colum and chromosome size in the second column

Tools specific parameters:
        --corMethod		STRING	CORRELATION METHOD FOR DEEPTOOLS PLOTCORRELATION PROCESS (default : spearman; valid options are spearman, pearson)
	--corPlot		STRING	DEEPTOOLS whatToPlot OPTION FOR PLOTCORRELATION PROCESS (default : heatmap; valid options are heatmap, scatterplot)
	--heatmap_width		INT	DEEPTOOLS heatmapWidth OPTION FOR PLOTHEATMAP PROCESS (default : 20)
	--matrix_downstream	INT	DEEPTOOLS downstream OPTION FOR COMPUTEMATRIX PROCESS (default : 2500)
	--matrix_upstream	INT     DEEPTOOLS upstream OPTION FOR COMPUTEMATRIX PROCESS (default : 2500)
	--gtf_id		STRING	GTF FEATURE TYPE FOR HOMER (default : "-gid" for gene_id ; valid options are "-gid" and "" for transcript_id which is HOMER default)
	--figType		STRING	FILE FORMAT FOR THE PLOT (default : pdf ; valid options are pdf,svg,ps,tiff,png)
	--intersect_options	STRING	INTERVENE --bedtools option (default : "-f 1E-9" ; see bedtools intersect documentation https://bedtools.readthedocs.io/en/latest/)
	--intersect_threshold	STRING	INTERVENE --intersect_thresh option (default : 1 ; see bedtools intersect documentation https://bedtools.readthedocs.io/en/latest/)

Pipeline parameters:
	--with_clustering	BOOL	If false then clustering step is skipped (default : true)
	--with_heatmap		BOOL	If false then heatmap step is skipped (default : true)
	--with_FR_heatmap	BOOL	If false then heatmap for reverse and forwad step is skipped (default : true ; cannot be true if --with_FR_bigwig is false))
	--with_FR_bigwig	BOOL	If false then bigwigs for reverse and forwad step is skipped (default : true)
	--with_peak_annot	BOOL	If false then peak annotation step is skipped (default : true)
	--with_plot_intersect	BOOL	If false then bed intersection step is skipped (default : true)
        --with_report           BOOL    If false then general report with multiqc  step is skipped (default : true)
        --publishdir_mode       STRING  MODE FOR EXPORTING PROCESS OUTPUT FILES TO OUTPUT DIRECTORY (default : "copy", must be "symlink", "rellink", "link", "copy", "copyNoFollow","move", see https://www.nextflow.io/docs/latest/process.html)

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


// Pipeline parameters information
def paramsSection() {
log.info """
==========================================================================
   SSDS post-process pipeline version 1.0 : Computes and Plots general
   statistics for processed Single-Stranded-DNA-Sequencing (SSDS) data
==========================================================================
** Main parameters ** 
Run name                            : ${params.name}
Sample name                         : ${params.sample_name}   
Filtered final peaks                : ${params.finalpeaksbed}
Reference hotspots                  : ${params.peakreference}
Filtered Type 1 reads               : ${params.bamfolder}/${params.bampattern}
Reference Type 1 reads              : ${params.bamreference}/${params.bampattern}
Bigwigs                             : ${params.bigwigfolder}/${params.bigwigpattern}
Bed files to overlap                : ${params.bedfolder}/${params.bedpattern}
Genome name                         : ${params.genome}
PublishDir mode                     : ${params.publishdir_mode}

** Tools parameters **
Correlation method for clustering   : ${params.corMethod}
Correlation plot type               : ${params.corPlot}
Heatmap width                       : ${params.heatmap_width}
Matrix downstream                   : ${params.matrix_downstream}
Matrix Upstream                     : ${params.matrix_upstream}
GTF feature base                    : ${params.gtf_id}
Figure type                         : ${params.figType}
Bedtools intersect options          : ${params.intersect_options}
Bedtools intersect threshold        : ${params.intersect_threshold}

** Pipeline steps **
Clustering                          : ${params.with_clustering}
Heatmap                             : ${params.with_heatmap}
Forward and Reverse  bigwigs        : ${params.with_FR_bigwig}
Forward and Reverse heatmaps        : ${params.with_FR_heatmap}
Peak annotation                     : ${params.with_peak_annot}
Plot intersection                   : ${params.with_plot_intersect}
HTML report			    : ${params.with_report}    
""".stripIndent()
}
// Print parameters in log report
paramsSection()

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
    //conda 'deeptools=3.5.1'
    conda "${baseDir}/conda_yml/environment_deeptools.yml"
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
    conda "${baseDir}/conda_yml/environment_deeptools.yml"
    //conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/clustering/plots", mode: params.publishdir_mode, pattern: "*.png*"
    publishDir "${params.outdir}/clustering/log",   mode: params.publishdir_mode, pattern: "*.log*"
    input:
        tuple path(finalpeaks_matrix), path(peakref_matrix) from bamSummary_ch
    output:
        path('*.png')
        path('*.log')
        val("ok") into clustering_ok
    when:
        params.with_clustering
    script:
    """
    plotCorrelation --corData $finalpeaks_matrix \
                    --corMethod ${params.corMethod} \
                    --whatToPlot ${params.corPlot} \
                    -o ${params.sample_name}.finalpeaks.plotCorrelation.${params.corMethod}.${params.corPlot}_mqc.png \
                    >& ${params.sample_name}.finalpeaks.plotCorrelation.${params.corMethod}.${params.corPlot}.log 2>&1
    plotPCA --corData $finalpeaks_matrix \
            --plotTitle ${params.sample_name} \
            --transpose \
            -o ${params.sample_name}.finalpeaks.plotPCA.${params.corMethod}.${params.corPlot}_mqc.png \
            >& ${params.sample_name}.finalpeaks.plotPCA.${params.corMethod}.${params.corPlot}.log 2>&1 
    # If another bed file is given to compare with hotspots
    if [ ${params.peakreference} != "None" ]
    then
        plotCorrelation --corData $peakref_matrix \
                        --corMethod ${params.corMethod} \
                        --whatToPlot ${params.corPlot} \
                        -o ${params.sample_name}.peakref.plotCorrelation.${params.corMethod}.${params.corPlot}_mqc.png \
                        >& ${params.sample_name}.peakref.plotCorrelation.${params.corMethod}.${params.corPlot}.log 2>&1
        plotPCA --corData $peakref_matrix \
                --plotTitle ${params.sample_name} \
                --transpose \
                -o ${params.sample_name}.peakref.plotPCA.${params.corMethod}.${params.corPlot}_mqc.png \
                >& ${params.sample_name}.peakref.plotPCA.${params.corMethod}.${params.corPlot}.log 2>&1
    fi
    """
}

//***************************************************************************//
//                         SECTION 2 : HEATMAPS                              //
//***************************************************************************//
// PROCESS 3    : COMPUTEMATRIX (DEEPTOOLS) 
// What it does : calculates scores per genome regions
// Input        : bed files containing coordinates of hotspots and filtered type1 BAM files and bigwig files from ssdsnextflowpipeline
// Output       : compressed matrix of values (.gzip  file)
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html?highlight=computematrix
process computeMatrix {
    tag "${params.sample_name}"
    label 'process_basic'
    conda "${baseDir}/conda_yml/environment_deeptools.yml"
    //conda 'deeptools=3.5.1'
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
    conda "${baseDir}/conda_yml/environment_deeptools.yml"
    //conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/heatmap/plots", mode: params.publishdir_mode, pattern: "*.png"
    publishDir "${params.outdir}/heatmap/log",   mode: params.publishdir_mode, pattern: "*.log"
    input:
        path(matrix) from matrix_ch
    output:
        path('*.png')
        path('*.log')
        val("ok") into heatmap_ok
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
                -o ${params.sample_name}_heatmap_mqc.png \
                >& ${params.sample_name}_heatmap.log 2>&1
    """
}                

// Here I create a channel grouping all bam files to an identifier
// The resulting channel is composed of 2 elements : [Id; path-to-bam]
Channel
    .fromPath( "${params.fragfolder}/${params.fragpattern}" ) 
    .map { tuple( it.name.split('_')[0..-3].join('_'), it ) }
    //.println() 
    .into{ fragF_ch ; fragR_ch }  

// PROCESS 5    : GETFORWARDSTRAND (SAMTOOLS AND DEEPTOOLS)
// What it does : Gets coverage originated from the forward strand from a bam file
// Input        : Filtered bam file of type1 ssds reads (from ssdsnextflowpipeline)
// Output       : forward coverage in bigwig format
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage
process getForwardStrand {
    tag "${frag_id}"
    label 'process_basic'
    conda "${baseDir}/conda_yml/environment_deeptools.yml"
    //conda 'deeptools=3.5.1 bioconda::samtools=1.14'
    publishDir "${params.outdir}/FRbigwig",     mode: params.publishdir_mode, pattern: "*.bigwig"
    publishDir "${params.fragfolder}/FRbed" , 	mode: params.publishdir_mode, pattern: "*.watson.bed"
    input:
        tuple val(frag_id), path(frag) from fragF_ch
    output:
        tuple val(frag_id), path('*bigwig') into fwd_bigwig_ch
    when:
        params.with_FR_bigwig
    script:
    """
    # filtering out the watson reads from fragment bedfile (+ strand)
    grep "+" $frag > ${frag.baseName}.watson.bed

    # preparing a bedgraph from watson-selected fragments, sorted by coordinates
    genomeCoverageBed -i ${frag.baseName}.watson.bed -g ${params.genome_size} -bga | sort -k1,1 -k2,2n > ${frag.baseName}.watson.bedgraph

    # converting bedgraph to bigwig
    bedGraphToBigWig ${frag.baseName}.watson.bedgraph ${params.genome_size} ${frag.baseName}.watson.bigwig


    """
}

// PROCESS 6    : GETREVERSESTRAND (SAMTOOLS AND DEEPTOOLS)
// What it does : Gets coverage originated from the reverse strand from a bam file
// Input        : Filtered bam file of type1 ssds reads (from ssdsnextflowpipeline)
// Output       : reverse coverage in bigwig format
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage
process getReverseStrand {
    tag "${frag_id}"
    label 'process_basic'
    conda "${baseDir}/conda_yml/environment_deeptools.yml"
    //conda 'deeptools=3.5.1 bioconda::samtools=1.14'
    publishDir "${params.outdir}/FRbigwig",     mode: params.publishdir_mode, pattern: "*.bigwig"
    publishDir "${params.fragfolder}/FRbed" ,   mode: params.publishdir_mode, pattern: "*.crick.bed"

    input:
        tuple val(frag_id), path(frag) from fragR_ch
    output:
        tuple val(frag_id), path('*bigwig') into rev_bigwig_ch
        val("ok") into FRbigwig_ok
    when:
        params.with_FR_bigwig
    script:
    """
    # filtering out the watson reads from fragment bedfile (+ strand)
    grep "-" $frag > ${frag.baseName}.crick.bed

    # preparing a bedgraph from watson-selected fragments, sorted by coordinates
    genomeCoverageBed -i ${frag.baseName}.crick.bed -g ${params.genome_size} -bga | sort -k1,1 -k2,2n > ${frag.baseName}.crick.bedgraph

    # converting bedgraph to bigwig
    bedGraphToBigWig ${frag.baseName}.crick.bedgraph ${params.genome_size} ${frag.baseName}.crick.bigwig


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
    tag "${frag_id}"
    label 'process_basic'
    conda "${baseDir}/conda_yml/environment_deeptools.yml"
    //conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/heatmap/matrices", mode: params.publishdir_mode, pattern: "*.matrix.FR"
    publishDir "${params.outdir}/heatmap/log",      mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(frag_id), path(rev_bw), path(fwd_bw) from FR_bigwig_ch
    output:
        tuple val(frag_id), path('*finalpeaks.matrix.FR') into matrix_FR_ch
        path('*.log')
    when:
        params.with_FR_heatmap && params.with_FR_bigwig
    script:
    """
    computeMatrix   reference-point --regionsFileName ${params.finalpeaksbed} \
                    --scoreFileName ${fwd_bw} ${rev_bw} \
                    --referencePoint=center \
                    --downstream=${params.matrix_downstream} \
                    --upstream=${params.matrix_upstream} \
                    --smartLabels \
                    --numberOfProcessors ${task.cpus} \
                    -o ${frag_id}_b${params.matrix_upstream}_a${params.matrix_downstream}.finalpeaks.matrix.FR \
                    >& ${frag_id}_finalpeaks.matrix.FR.log 2>&1
    """
}

// PROCESS 8    : PLOTHEATMAPFR (DEEPTOOLS)
// What it does : creates a heatmap for scores associated with hotspots
// Input        : Matrix file from the computeMatrixFR process
// Output       : heatmaps in pdf format showing forward and reverse signal around hotspots
// Resources    : https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html?highlight=plotheatmap
process plotHeatmapFR {
    tag "${frag_id}"
    label 'process_basic'
    conda "${baseDir}/conda_yml/environment_deeptools.yml"
    //conda 'deeptools=3.5.1'
    publishDir "${params.outdir}/heatmap/plots", mode: params.publishdir_mode, pattern: "*.png"
    publishDir "${params.outdir}/heatmap/log",   mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(frag_id), path(matrix) from matrix_FR_ch
    output:
        path('*.png')
        path('*.log')
        val("ok") into FRheatmap_ok
    when:
        params.with_FR_heatmap && params.with_FR_bigwig
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
                -o ${frag_id}_heatmap.FR_mqc.png \
                >& ${frag_id}_heatmap.FR.log 2>&1
    """
}                
//***************************************************************************//
//                         SECTION 3 : PEAK ANNOTATION                       //
//***************************************************************************//
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
        val("ok") into peakannot_ok
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
                     ${params.gtf_id} \
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
//***************************************************************************//
//                         SECTION 3 : PEAK COMPARISON                       //
//***************************************************************************//
// PROCESS 10    : PLOT_INTERSECT (INTERVENE)
// What it does : Compute overlap of bed sets using intervene
// Input        : list of bed files
// Output       : overlapping bed files and upset plot of overlaps
// Resource     : https://intervene.readthedocs.io/en/latest/modules.html
process plotIntersect {
    tag "${params.sample_name}"
    label 'process_basic'
    //conda 'bioconda::intervene=0.6.4'
    conda "${params.conda_intervene}"
    publishDir "${params.outdir}/intersect/plots",   mode: params.publishdir_mode, pattern: "*"
    //publishDir "${params.outdir}/intersect/overlap", mode: params.publishdir_mode, pattern: "*.bed"
    publishDir "${params.outdir}/intersect/log",     mode: params.publishdir_mode, pattern: "*.log"
    output:
        path('*')
        val("ok") into intersect_ok
    when:
        params.with_plot_intersect
    script:
    """
    #Run intervene to compute the intesecting sets
    intervene upset -i ${params.bedfolder}/${params.bedpattern} \
                    --scriptonly --figtype ${params.figType} \
                    --overlap-thresh ${params.intersect_threshold} \
                    --project ${params.sample_name}_mqc \
                    -o ./${params.sample_name}_mqc \
                    --save-overlaps >& ${params.sample_name}_intervene.log 2>&1
                    #--bedtools-options ${params.intersect_options} \

    #There is a bug (need to investigate but no time now) : need to edit the R script generated by intervene before execution otherwise the execution is halted
    sed -i "s/, mainbar.y.label =\\"No. of Intersections\\", sets.x.label =\\"Set size\\")/)/g" ./${params.sample_name}_mqc/${params.sample_name}_mqc_upset.R >& ${params.sample_name}_intervene_sed.log 2>&1

    #Finally, esecute Rscript to plot the upset plot
    Rscript ./${params.sample_name}_mqc/${params.sample_name}_mqc_upset.R >& ${params.sample_name}_intervene_rscript.log 2>&1
    """
}

//***************************************************************************//
//                         SECTION 4 : GENERAL REPORT                        //
//***************************************************************************//
if (!params.with_clustering) { clustering_ok = Channel.value( 'ok' ) }
if (!params.with_heatmap) { heatmap_ok = Channel.value( 'ok' ) }
if (!params.with_FR_heatmap) { FRheatmap_ok = Channel.value( 'ok' ) }
if (!params.with_FR_bigwig) { FRbigwig_ok = Channel.value( 'ok' ) }
if (!params.with_peak_annot) { peakannot_ok = Channel.value( 'ok' ) }
if (!params.with_plot_intersect) { intersect_ok = Channel.value( 'ok' ) }

process generalReport {
    tag "${params.sample_name}"
    label 'process_basic'
    conda "${baseDir}/conda_yml/environment_multiqc.yml"
    publishDir "${params.outdir}/multiqc",   mode: params.publishdir_mode 
    input:
        val('clustering_ok') from clustering_ok.collect().ifEmpty([]) 
        val('heatmap_ok') from heatmap_ok.collect().ifEmpty([])   
        val('FRheatmap_ok') from FRheatmap_ok.collect().ifEmpty([])
        val('FRbigwig_ok') from FRbigwig_ok.collect().ifEmpty([])
        val('peakannot_ok') from peakannot_ok.collect().ifEmpty([])
        val('intersect_ok') from intersect_ok.collect().ifEmpty([])
    output:
        path('*')
    when:
        params.with_report
    script:
    """
    multiqc --export -c ${params.multiqc_configfile} -n ${params.sample_name}.multiqc.report \
        ${params.outdir}/*
    """
} 

//***************************************************************************//
//                                                                           //
//                          END OF PIPELINE !!                               //
//                                                                           //
//***************************************************************************//

// PRINT LOG MESSAGE ON COMPLETION        
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Command line: $workflow.commandLine"
    println "Script ID: $workflow.scriptId"
    println "Run name: $workflow.runName"

}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    println "Error report: ${workflow.errorReport}"
}




