#!/bin/bash
R -e 'library(devtools);devtools::document()';
R -e 'Rcpp::compileAttributes(".",verbose=TRUE)';
R CMD build ../fluEvidenceSynthesis
