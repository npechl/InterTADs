# nikopechThesis/Genomic Data Integration

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nikopech/data-integration-model.git/master?urlpath=rstudio)

## Getting Started

Before running the algorithm, the following steps need to be followed

### Prerequisites

In order to run the integration algorithm, the following packages are required:

- IlluminaHumanMethylation450kanno.ilmn12.hg19 (https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylation450kanno.ilmn12.hg19.html)
- GenomicFeatures (https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
- systemPipeR (https://bioconductor.org/packages/release/bioc/html/systemPipeR.html)
- biomaRt (https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- vcfR (https://cran.r-project.org/web/packages/vcfR/index.html)
- stringr / stringi
- data.table
- gtools

For plotting, the following packages are also required:

- karyoploteR (http://bioconductor.org/packages/release/bioc/html/karyoploteR.html)
- magick (https://cran.r-project.org/web/packages/magick/index.html)
- png (https://cran.r-project.org/web/packages/png/index.html)

### Installing

Download the files in a directory, set the directory as your working directory, and run the following command into Rstudio's console

```
source("bioCombine.R")
```

## Documentation

```
bioCombine(biodata, colCmb = NULL, scale = 100, chromosomes = NULL, txdb)
```
| Property    | Type            | Default | Description |
|:------------|:----------------|:--------|:------------|
| ```biodata``` | list of strings | required | List of all file paths that are going to be given as input to the algorithm (for the format of this variable, please see bellow). |
| ```colCmb``` | data table or data frame          | NULL     | Table for column connection. If NULL, column integration is not required (for the format of this variable, please see bellow). |                                                                                                                             
| ```scale``` | Numeric         | 100      | Number used for scaling of inputs. |                                                                                              
| ```chromosomes``` | string          | NULL     | Chromosome that is going to be used (given as a number, not 'chr1').|  
| ```txdb``` | TxDB object | required     | Database to be used for getting genomic features (gene id, gene locus).| 

* ```biodata```: Α list of strings that shows all paths to file-inputs. Each element of the list must be a unique table.
* ```colCmb```: Α table used for column connection between input tables:
  + First column refers to the index of the table to which operations are going to be applied.
  + Second column contains all column names that are going be processed.
  + Third column contains the operation that is going to be applied.
  + Fourth column contains weights that are going to be used.
  + Fifth column contains new names for the columns.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
