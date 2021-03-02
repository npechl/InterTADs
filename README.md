# InterTADs

InterTADs is an open-source tool written in [R](https://www.r-project.org/), for integrating multi-omics data (e.g. DNA methylation, expression, mutation) from the same physical source (e.g. patient) taking into account the chromatin configuration of the genome, i.e. the topologically associating domains (TADs).

## Installation

You can simply clone the repository by using [git](https://git-scm.com/):

```
git clone https://github.com/BiodataAnalysisGroup/InterTADs.git
```

Before running any scripts, make sure the following packages are installed in your machine:

```R
install.packages(c("data.table", 
                   "tidyverse", 
                   "gplots", 
                   "png", 
                   "gghalves"))

devtools::install_github("stephenturner/annotables")
```

...and from [Bioconductor](https://www.bioconductor.org/):

```R
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene", 
                       "TxDb.Hsapiens.UCSC.hg38.knownGene", 
                       "GenomicRanges", 
                       "org.Hs.eg.db", 
                       "systemPipeR", 
                       "karyoploteR"))
```

## Usage

There are three main scripts for integrating your multi-omics data:

* ```Data_Integration.R```
* ```TADiff.R```
* ```Visualization.R```

### Data Integration

For the Data Integration part, all datasets are separated into two folders, ```freq``` and ```counts```, based on the information they are carrying (frequency or score count values). 

The two folders are placed into a directory, along with a meta-data file which provides information about the mapping between the columns for each dataset. For more details regarding the structure of this file please see [here](Datasets/meta-data.csv).

The script allows the user to define different folder (or file) names. Moreover, the user can choose a folder name for the output table and a option about the Human Genome that is being used (accepted values are [`hg19`](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/) or [`hg38`](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39)).

Once every input is provided, the script can be run by:

```
source("Data_Integration.R")
```

### TADiff

For the TADiff part, the paths to the input and output folders must be provided. Also a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file is needed containing information about the TADs. In order to run the script:

```
source("TADiff.R")
```

### Visualization

For the visualization of the results, the paths to input and output data need to be provided:

```
source("Visualization.R")
```

## Data 

The proposed method was evaluated on data from Chronic lymphocytic leukemia (DNA methylation and expression values). The datasets have been deposited in the [ArrayExpress](www.ebi.ac.uk/arrayexpress) database at EMBL‐EBI under the accession numbers E‐MTAB‐6955 and E‐MTAB‐6962, respectively.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the [MIT](https://opensource.org/licenses/MIT) License - see the [LICENSE](LICENSE) file for details