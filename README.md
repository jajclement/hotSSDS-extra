# :chart_with_upwards_trend: **SSDS post-process pipeline version 1.0** :bar_chart:
**Compute and Plot general statistics for processed Single-Stranded-DNA-Sequencing (SSDS) data**

## **Context**  
**SSDS method** was originally published by [Khil et al., 2012](https://genome.cshlp.org/content/22/5/957.long).  
The objective is to map **double-strand breaks** (DSBs) along the genome.   
In this method, chromatin is extracted from adult testes and then immunoprecipitated with an antibody against **DMC1 protein**, which is a meiosis-specific recombinase. DMC1 covers the **single-stranded DNA** resulting from the **resection of double-strand breaks (DSBs)**. SSDS uses the ability of single-stranded DNA to form hairpins. 

The raw reads from DMC1 ChIP-Seq sequencing need to be processed using [SSDS nextflow pipeline](https://gitlab.igh.cnrs.fr/pauline.auffret/ssdsnextflowpipeline).
In this pipeline, the final peaks, the type1 mapped reads and bigwig produced by **SSDS nextflow pipeline** will be used to compute and plot general statistics data about the signal.

The pipeline uses [Nextflow]( https://www.nextflow.io/) >= 21.10.0  

**In details, the pipeline is composed of 10 processes :**            
- **PROCESS 1**    : MULTIBAMSUMMARY (DEEPTOOLS)
- **PROCESS 2**    : PLOTCORRELATION (DEEPTOOLS)
- **PROCESS 3**    : COMPUTEMATRIX (DEEPTOOLS)
- **PROCESS 4**    : PLOTHEATMAP (DEEPTOOLS)
- **PROCESS 5**    : GETFORWARDSTRAND (SAMTOOLS AND DEEPTOOLS)
- **PROCESS 6**    : GETREVERSESTRAND (SAMTOOLS AND DEEPTOOLS)
- **PROCESS 7**    : COMPUTEMATRIXFR (DEEPTOOLS)
- **PROCESS 8**    : PLOTHEATMAPFR (DEEPTOOLS)
- **PROCESS 9**    : ANNOTATEPEAKS (HOMER AND R)
- **PROCESS 10**   : PLOT_INTERSECT (INTERVENE)
- **PROCESS 11**   : GENERAL_REPORT (MULTIQC)


## **How to run the pipeline on IGH cluster**
### 0. Requirements
Get [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) if not installed on your system. First, download the installation script (for example in ``/home/${USER}/work/bin`` directory) :
````sh
mkdir -p /home/${USER}/work/bin
cd /home/${USER}/work/bin
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
````
Then execute the installation script and follow the prompts on the installer screens :
````sh
bash Miniconda3-latest-Linux-x86_64.sh
````

### 1. Get the pipeline and set up conda environment
First you need to clone the pipeline in your working directory (in the following instructions, ``/home/${USER}/work/`` will refer to your working directory. Please substitute with the according path if different) :
````sh
cd /home/${USER}/work
git clone https://gitlab.igh.cnrs.fr/pauline.auffret/ssdspostprocess.git
cd ssdspostprocess
````
Then install pipeline conda environment through sbatch script : 
```` 
sbatch -p computepart -J "install_conda_env" --export=ALL --mem 5G -t 5-0:0 --wrap "bash src/install_pipeline.sh"
```` 
Please use ``bash src/install_pipeline.sh -h`` to see details.   
This will create 1 conda environment : **nextflow21** containing **nextflow version 21.10.0**.  
You can use your own conda environment as long as it contains that version of Nextflow.   

### 2. Pipeline configuration 
There are currently **3 configuration files** :
- ````./conf/igh.config```` contains cluster resources requirements & reference genomes info. You don't need to edit this file unless you want to custom the requirements for CPU/memory usage and compute queue (see [IGH cluster documentation](https://kojiki.igh.cnrs.fr/doku.php?id=cluster,)). If a new genome is available on the cluster and does not appear in this file, please contact me or Aubin Thomas, manager of bioinformatics resources at the IGH. If you are running the pipeline on an other computing cluster, you need to specify the relevant configuration file.
- ````./nextflow.config```` contains default pipeline parameters (it's better to **not** edit this file, default parameters will be overwritten by your custom json parameter file, see next point).
- ````./conf/mm10.json```` custom parameters file, you can use your own json file following the same template.    
    
The complete list of parameters is accessible through the command :
````
cd /home/${USER}/work/ssdspostprocess
conda activate nextflow21
nextflow run main.nf --help
````
````
=============================================================================
Input data parameters:
        -params_file            FILE    PATH TO PARAMETERS JSON FILE (template and default : /work/demassyie/ssdsnextflowpipeline/conf/mm10.json)
        --name                  STRING  ANALYSIS NAME (default : "SSDS_postprocess_pipeline")
        --sample_name           STRING  SAMPLE OR GROUP NAME (default : "DMC1-ChIP")
        --finalpeaksbed         FILE    FILTERED PEAKS IN BED FORMAT (for example, coming from ssdsnextflowprocess, in the finalpeaks folder (warning : if finalpeakbed is not None, bamreference cannot be None) ; default : )
        --peakreference         FILE    REFERENCE PEAKS FILE IN BED FORMAT (set to "None" if none provided (warning : if finalpeakbed is not None, bamreference cannot be None) : set of peaks to compare the finalpeaksbed to (for example, hotspots from B6 WT mouse) ; default : )
        --bamfolder             DIR     ABSOLUTE PATH TO BAM FOLDER CONTAINING FILTERED TYPE 1 BAM FILES FROM SSDSNEXTFLOWPIPELINE (for example, coming from ssdsnextflowprocess, in the bwa/filterbam/flag_*/parse_itr/type1/bam folder ; default :
        --bampattern            REGEX   PATTERN FOR MATCHING BAM FILES IN BAMFOLDER (default : "*.bam")
        --bamreference          DIR     ABSOLUTE PATH TO REFERENCE BAM FOLDER CONTAINING FILTERED TYPE 1 BAM FILES FOR REFERENCE (set to "None" if none provided (warning : if bamreference is not None, peakreference cannot be None)
        --bigwigfolder          DIR     ABSOLUTE PATH TO BIGWIG FOLDER CONTAINING BIGWIG FILES FROM SSDSNEXTFLOWPIPELINE (for example, coming from ssdsnextflowprocess, in the bigwig/kbrick_bigwig/deeptools/binsize* folder ; default :
        --bigwigpattern         REGEX   PATTERN FOR MATCHING BIGWIG FILES IN BIGWIGFOLDER (default : "*ssDNA_type1.deeptools.RPKM.bigwig")
        --bedfolder             DIR     ABSOLUTE PATH TO BED FOLDER CONTAINING BED TO COMPUTE INTERSECT (default :
        --bedpattern            REGEX   PATTERN FOR MATCHING BED FILES IN BED FOLDER (default : "*.bed")

Genome parameters:
        --genomebase            DIR     PATH TO REFERENCE GENOMES (default : "/poolzfs/genomes")
        --genome                STRING  REFERENCE GENOME NAME (must correspond to an existing genome in your config file, default : "mm10")
        --genomedir             DIR     PATH TO GENOME DIRECTORY (required if your reference genome is not present in your config file)
        --genome_fasta          FILE    PATH TO FILE GENOME FASTA FILE WITH PREEXISTING INDEX FILES FOR BWA (required if your reference genome is not present in your config file)
        --genome_gtf

Tools specific parameters:
        --corMethod             STRING  CORRELATION METHOD FOR DEEPTOOLS PLOTCORRELATION PROCESS (default : spearman; valid options are spearman, pearson)
        --corPlot               STRING  DEEPTOOLS whatToPlot OPTION FOR PLOTCORRELATION PROCESS (default : heatmap; valid options are heatmap, scatterplot)
        --heatmap_width         INT     DEEPTOOLS heatmapWidth OPTION FOR PLOTHEATMAP PROCESS (default : 20)
        --matrix_downstream     INT     DEEPTOOLS downstream OPTION FOR COMPUTEMATRIX PROCESS (default : 2500)
        --matrix_upstream       INT     DEEPTOOLS upstream OPTION FOR COMPUTEMATRIX PROCESS (default : 2500)
        --gtf_id                STRING  GTF FEATURE TYPE FOR HOMER (default : "-gid" for gene_id ; valid options are "-gid" and "" for transcript_id which is HOMER default)
        --figType               STRING  FILE FORMAT FOR THE PLOT (default : pdf ; valid options are pdf,svg,ps,tiff,png)
        --intersect_options     STRING  INTERVENE --bedtools option (default : "-f 1E-9" ; see bedtools intersect documentation https://bedtools.readthedocs.io/en/latest/)
        --intersect_threshold   STRING  INTERVENE --intersect_thresh option (default : 1 ; see bedtools intersect documentation https://bedtools.readthedocs.io/en/latest/)

Pipeline parameters:
        --with_clustering       BOOL    If false then clustering step is skipped (default : true)
        --with_heatmap          BOOL    If false then heatmap step is skipped (default : true)
        --with_FR_heatmap       BOOL    If false then heatmap for reverse and forwad step is skipped (default : true ; cannot be true if --with_FR_bigwig is false))
        --with_FR_bigwig        BOOL    If false then bigwigs for reverse and forwad step is skipped (default : true)
        --with_peak_annot       BOOL    If false then peak annotation step is skipped (default : true)
        --with_plot_intersect   BOOL    If false then bed intersection step is skipped (default : true)
        --with_report           BOOL    If false then general report with multiqc  step is skipped (default : true)
        --publishdir_mode       STRING  MODE FOR EXPORTING PROCESS OUTPUT FILES TO OUTPUT DIRECTORY (default : "copy", must be "symlink", "rellink", "link", "copy", "copyNoFollow","move", see https://www.nextflow.io/docs/latest/process.html)

=============================================================================

````

**WARNING** the option ``--with-report`` is not stable due to conda issue, set it to false if not working.   

One **important thing to note**, in Nextflow command lines, the **native options** are preceded with one **single hyphen** (e.g. ``-profile``), while **parameters specific to SSDS pipeline** are preceded with **2 hyphens** (e.g. ``--genome 'mm10'``).   

### 4. Reference genome
The reference genome should be in the ``/poolzvs/genomes`` directory on IGH cluster. Currently, the available genomes are mm10, hg19, hg38, sacCer2, sacCer3, dm3, dm6.

You can use ````--genome 'mm10'```` and you won't need to worry about the other genome parameters like fasta path etc.

But if you want to use **another reference**, you will need to set the following parameters : 
- absolute path to genome ````--genomedir /path/to/genome````
- absolute path to fasta file. **Indexes for BWA SHOULD EXIST in the same directory** ````--genome_fasta /path/to/genome.fa````
- the name of the genome ````--genome_name mm11````
- absolute path to the gtf file ````--genome_gtf /path/to/genome.gtf````

You can create your own json parameter file.

### 4. Run the pipeline !
   
**Run through bash script ``bash run_pipeline.sh [options]`` (RECOMMENDED)**    

Usage: bash run_pipeline.sh [options]
Options :
    * ``-h`` display help message     
    * ``-g`` Absolute path to the genome config file (default : /home/${USER}/work/ssdspostprocess/conf/mm10.json)    
    * ``-p`` Absolute path to ssds postprocess pipeline base directory (default : /home/${USER}/work/ssdspostprocess)    
    * ``-b`` Absolute path to base directory where the output directory will be created (default : /home/${USER}/work/results)    
    * ``-n`` Analysis name (default : SSDS_postprocess) INFO : by default, this parameter will match the --name option in nextflow command line    
    * ``-c`` Absolute path to IGH cluster configuration file (default : /home/${USER}/work/ssdspostprocess/conf/igh.config)    
    * ``-a`` Absolute path to conda environment for nextflow (default : /home/${USER}/work/bin/miniconda3/envs/nextflow_dev)    
    * ``-o`` Optional arguments for the pipeline (for example ``"--corMethod spearman --with_heatmap false"`` ;  default : "")    
    * ``-w`` Valid Nextflow Tower token (default : None ; if not None, then the option ``-with-tower`` has to be added in -o parameter))    
    * ``-t`` set to 1 if running pipeline on test data located in /home/${USER}/ssdspostprocess/tests/fastq (default : 0)   
    * ``-f`` set to 1 to force pipeline to run without checking resume/output directory (default : 0)   
INFO : the output directory will be located in the base directory and will be named after the analysis name parameter with the .outdir suffix (default /home/${USER}/work/results/SSDS_postprocess.outdir)    

Please provide all files and directories with **absolute paths**.   
The command line should look like :
 ````
bash run_pipeline.sh -g custom_params.json -n your-analysis-name -b my-base-directory -o "-resume"
````
The results will be located in a directory named after the analyis name (-n argument) suffixed with .outdir, in the base directory.

### 5. Monitor the pipeline with Nextflow Tower
You can use ``-with-tower`` option to monitor your jobs through [nextflow tower web interface](https://tower.nf/).   
You first need to sign in to get your key, then add it to your parameters with the ``-w "your_key"`` option (see 4.1)     

### 6. Output
Tree overview of the output folder composition :
.
├── nxfReports   
├── slurm   
├── annotations  
│   ├── plots   
│   └── log   
├── FRbigwig   
│   └── log   
├── multiqc   
│   ├── *.multiqc.report_plots   
│   └── *.multiqc.report_data   
├── clustering    
│   ├── plots   
│   ├── log   
│   └── matrices   
└── heatmap   
    ├── plots   
    ├── log   
    └── matrices
   

:christmas_tree:
