#!/bin/bash

# Always run at least once
#R -e 'Rcpp::compileAttributes(".",verbose=TRUE);library(devtools);devtools::check()'
R -e 'Rcpp::compileAttributes(".",verbose=TRUE);library(devtools);devtools::test()'

# This script depends on inotify-hookable
inotify-hookable -t 1000 -i '*.o' -i '\..*' -w data -w src -w tests -f DESCRIPTION -w R -c "R -e 'Rcpp::compileAttributes(\".\",verbose=TRUE);library(devtools);devtools::test()' && sleep 10"
#inotify-hookable -t 1000 -i '*.o' -i '\..*' -w data -w src -w tests -f DESCRIPTION -w R -c "R -e 'Rcpp::compileAttributes(\".\",verbose=TRUE);library(devtools);devtools::check()' && sleep 10"
