# :tropical_fish: **SSDS post-process pipeline version 1.0** :bar_chart:
**Compute and Plot general statistics for processed Single-Stranded-DNA-Sequencing (SSDS) data**

## **Context**  
**SSDS method** was originally published by [Khil et al., 2012](https://genome.cshlp.org/content/22/5/957.long).  
The objective is to map **double-strand breaks** (DSBs) along the genome.   
In this method, chromatin is extracted from adult testes and then immunoprecipitated with an antibody against **DMC1 protein**, which is a meiosis-specific recombinase. DMC1 covers the **single-stranded DNA** resulting from the **resection of double-strand breaks (DSBs)**. SSDS uses the ability of single-stranded DNA to form hairpins. 

The raw reads from DMC1 ChIP-Seq sequencing need to be processed using [SSDS nextflow pipeline](https://gitlab.igh.cnrs.fr/pauline.auffret/ssdsnextflowpipeline).
In this pipeline, the final peaks, the type1 mapped reads and bigwig produced by **SSDS nextflow pipeline** will be used to compute and plot general statistics data about the signal.

The pipeline uses [Nextflow]( https://www.nextflow.io/) > 20.04.1  

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

### 1. Get the pipeline and set up conda environment
First you need to clone the pipeline in your working directory (in the following instructions, ``/home/${USER}/work/`` will refer to your working directory. Please substitute with the according path if different) :
````sh
cd /home/${USER}/work
git clone https://gitlab.igh.cnrs.fr/pauline.auffret/ssdsnextflowpipeline.git
cd ssdsnextflowpipeline
````
Then install pipeline conda environment through sbatch script : 
```` 
sbatch -p computepart -J "install_conda_env" --export=ALL --mem 5G -t 5-0:0 --wrap "bash src/install_pipeline.sh"
```` 
Please use ``bash src/install_pipeline.sh -h`` to see details.   
This will create 1 conda environment : **nextflow_dev** containing **nextflow version 20.04.01**.  
You can use your own conda environment as long as it contains that version of Nextflow.   

