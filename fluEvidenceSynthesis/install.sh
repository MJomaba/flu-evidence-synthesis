#!/bin/bash

R_EXEC=R

for i in "$@"
do
case $i in
    -p=*|--path=*)
    R_PATH="${i#*=}"
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done

if [ -n "$R_PATH" ]; then
    #echo "R_PATH  = ${R_PATH}";
    R_EXEC="${R_PATH}/${R_EXEC}"
fi
#echo "R_EXEC = ${R_EXEC}";

if command -v ${R_EXEC} 2>/dev/null; then
    ${R_EXEC} -e 'library(devtools);devtools::document()';
    ${R_EXEC} -e 'Rcpp::compileAttributes(".",verbose=TRUE)';
    # Should actually move needed header in inst/include/fluEvidenceSynthesis.h insted of using echo
    #echo -e "#include \"rcppwrap.hh\"\n$(cat src/RcppExports.cpp)" > src/RcppExports.cpp;
    ${R_EXEC} CMD INSTALL ../fluEvidenceSynthesis
else
    echo "Cannot locate the executable ${R_EXEC}. Make sure it is installed. If installed, but not in path use --path=/path/to to specify the path, e.g.:";
    echo ""
    echo "install.sh --path=/usr/local/bin"
    echo ""
fi
