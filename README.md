## List of custom script plugins for analysis
Ajit Johnson Nirmal

## R
#### Use these scripts by:
~~~~
library(devtools)
library(roxygen2)
source_url ("")
~~~~

### Merge multiple transcripts into one gene
~~~~
source_url("https://raw.githubusercontent.com/ajitjohnson/Jupyter-Notebooks/master/r_scripts/collapse_to_one_gene.R")
genesummary (x)
~~~~

### After normalization, this function can be used to plot the before & after comparision (barplot)
~~~~
source_url("https://raw.githubusercontent.com/ajitjohnson/Jupyter-Notebooks/master/r_scripts/plot_before_after_normalization.R")
plotba (x,y)
~~~~

### Go enrichment analysis
Requires output from a DESeq analysis with a column name "padj"
~~~~
source_url("https://raw.githubusercontent.com/ajitjohnson/Jupyter-Notebooks/master/r_scripts/goenrichment.R")
goenrichment (x)
~~~~
---
---
---
## Python
#### Use these scripts by:
Run file_url and then run the following command
~~~~
import wget
exec(open(wget.download(file_url)).read())
~~~~

### Preprocess CycIF data: Normalization
Requires a pandas dataframe with cells on rows and markers/genes on columns.
~~~~
file_url = 'https://raw.githubusercontent.com/ajitjohnson/Jupyter-Notebooks/master/py_scripts/mi_pp_normalize.py'
exec(open(wget.download(file_url)).read())
mi_pp_normalize (x)
~~~~

### Preprocess CycIF data: Convert dataframe into AnnData (Annotated Data)
Requires a pandas dataframe with cells on rows and markers/genes on columns.
~~~~
file_url = 'https://raw.githubusercontent.com/ajitjohnson/Jupyter-Notebooks/master/py_scripts/mi_pp_anndata.py'
exec(open(wget.download(file_url)).read())
mi_pp_anndata (x)
~~~~
