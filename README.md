## flu-evidence

This repository contains the code used for flu data analysis first published in

Baguelin, Marc, Stefan Flasche, Anton Camacho, Nikolaos Demiris, Elizabeth Miller, and W. John Edmunds. ‘Assessing Optimal Target Populations for Influenza Vaccination Programmes: An Evidence Synthesis and Modelling Study’. PLoS Med 10, no. 10 (8 October 2013): e1001527. doi:10.1371/journal.pmed.1001527.

The code is distributed under a GPLv3 license. Please remember to cite the above manuscript when using this code.

## Install

This code is now mainly developed as an R package and we recommend it to install it as such.

If you haven’t done so yet, first clone this repository:

```
git clone https://github.com/MJomaba/flu-evidence-synthesis.git
cd flu-evidence-synthesis
```

If you have already cloned it then cd into the directory and update it:

```
git pull
```
Install the dependencies for the package, first start R and in R run:

```
install.packages(c("Rcpp", "BH", "RcppEigen"))
```

Now build the R package (from the commandline):

```
R CMD build fluEvidenceSynthesis
```

Then open R again and install the created package and its dependencies. In R run:

```
install.packages("fluEvidenceSynthesis_1.0.tar.gz",repos=NULL)
```

## Usage examples

Coming soon!
