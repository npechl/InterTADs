# Genomic Data Integration

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fpsom/nikopechThesis.git/master?urlpath=rstudio)

## Getting Started

Before running the algorithm, please follow the instructions below

### Prerequisites

In order to run the integration algorithm, the following packages are required:

- [```IlluminaHumanMethylation450kanno.ilmn12.hg19```](https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylation450kanno.ilmn12.hg19.html)
- [```biomaRt```](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- [```vcfR```](https://cran.r-project.org/web/packages/vcfR/index.html)
- [```stringr```](https://cran.r-project.org/web/packages/stringr/index.html)/[```stringi```](https://cran.r-project.org/web/packages/stringi/index.html)
- [```data.table```](https://cran.r-project.org/web/packages/data.table/index.html)
- [```gtools```](https://cran.r-project.org/web/packages/gtools/index.html)

For plotting, the following packages are also required:

- [```karyoploteR```](http://bioconductor.org/packages/release/bioc/html/karyoploteR.html)
- [```magick```](https://cran.r-project.org/web/packages/magick/index.html)
- [```png```](https://cran.r-project.org/web/packages/png/index.html)

### Installing

1. Download the files in a directory

2. Set the directory as your working directory

3. Run the following command at R console

```
source("bioCombine.R")
```

## Documentation

```
bioCombine(biodata, colCmb = NULL, scale = NULL, chromosomes = NULL, txdb)
```
| Property    | Type            | Default | Description |
|:------------|:----------------|:--------|:------------|
| ```biodata``` | list of strings | required | List of all file paths that are going to be given as input to the algorithm (for the format of this variable, please see bellow). |
| ```colCmb``` | data table or data frame          | NULL     | Table for column connection. <br /> If NULL, column integration is not required (for the format of this variable, please see bellow). |
| ```scale``` | Numeric         | NULL      | Number used for scaling of inputs. <br /> If NULL, no scaling will be applied. |                                                                         
| ```chromosomes``` | string          | NULL     | Chromosome that is going to be used (given as a number, not 'chr1'). <br /> If NULL, all chromosomes found in data will be used.|  
| ```txdb``` | TxDB object | required     | Database to be used for getting genomic features (gene id, gene locus).| 

* ```biodata```: It is a list of strings that shows all paths to file-inputs. Each element of the list must be a table.
* ```colCmb```: It is a table showing information about the column connection that is going to be applied between input tables.
  + First column refers to the index of the table to which operations are going to be applied.
  + Second column contains all column names that are going be processed.
  + Third column contains the operation that is going to be applied.
  + Fourth column contains weights that are going to be used.
  + Fifth column contains new names for the columns.

## Usage

1. Create a TxDB object. Find instructions [here](https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf):

    ```mydb = makeTxDBfromGFF("path/to/GFF/file")```

2. Construct paths to files: 

    ```
    paths = list(c("path/to/CSV/file"),
                 c("path/to/TXT/file1", "path/to/TXT/file2"),
                 c("path/to/VCF/indel/file1", "path/to/VCF/snp/file1", 
                   "path/to/VCF/indel/file2", "path/to/VCF/snp/file2"))
    ```
                    
3. Call integration function:

    ```post_table = bioCombine(biodata = paths, chromosomes = 22, txdb = mydb)```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
