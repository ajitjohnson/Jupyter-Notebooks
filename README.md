# List of custom script plugins for analysis
#### Ajit Johnson Nirmal

## R
#### Use these scripts by:
library(devtools) </br>
library(roxygen2) </br>
source_url("") </br>

## Merge multiple transcripts into one gene
genesummary (x)
source_url("https://raw.githubusercontent.com/ajitjohnson/Jupyter-Notebooks/master/r_scripts/collapse_to_one_gene.R")

## After normalization, this function can be used to plot the before & after comparision (barplot)
plotba (x,y)
source_url("https://raw.githubusercontent.com/ajitjohnson/Jupyter-Notebooks/master/r_scripts/plot_before_after_normalization.R")

## Go enrichment analysis
## Requires output from a DESeq analysis with a column name "padj"
goenrichment (x)
source_url("https://raw.githubusercontent.com/ajitjohnson/Jupyter-Notebooks/master/r_scripts/goenrichment.R")


# Python
### Use these scripts by:
### Run file_url and then run the following command
import wget
exec(open(wget.download(file_url)).read())

## Preprocess CycIF data: Normalization
## Requires a pandas dataframe with cells on rows and markers/genes on columns
mi_pp_normalize (x)
file_url = 'https://raw.githubusercontent.com/ajitjohnson/Jupyter-Notebooks/master/py_scripts/mi_pp_normalize.py'
