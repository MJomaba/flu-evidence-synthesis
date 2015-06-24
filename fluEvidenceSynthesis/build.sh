#!/bin/bash
R -e 'Rcpp::compileAttributes(".",verbose=TRUE)';
R CMD build ../fluEvidenceSynthesis
