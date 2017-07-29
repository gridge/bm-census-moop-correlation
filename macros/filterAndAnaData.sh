#!/bin/bash
# Helper to quickly filter data using R and analyze with ROOT program

if [ -z "$4" ]; then
    echo "Usage: $0 input_census_data input_moop_map category category_values [category_exception [output_label]]"
    echo "Launch from the main repository folder (assumes analysis executable in 'build/', outputs are stored in the 'results/' folder."
    echo " input_census_data: CVS input file with census data"
    echo " input_moop_map: input Moop map XML file"
    echo " category: category label to select"
    echo " category_values: string with format as R vector with category values to consider for category == 0 (e.g. 'c(\"yes\")')"
    echo " category_exception: string with format as R vector with category values to *not* consider (exclude)"
    echo " output_label: output label (default: same as category)"
    echo ""
    exit
fi

#Start filtering data
inputData=$1
inputMoop=$2
category=$3
catValues=$4
catExclude='c()' #default: empty R vector
outputLabel="${category}"

#echo $inputData
#echo $inputMoop
#echo $category
#echo $catValues
#echo $catExclude

if ! [ -z "$5" ]; then
    catExclude=$5
fi
if ! [ -z "$6" ]; then
    outputLabel=$6
fi

#Start filtering input data
echo "====> Filtering input Census data"
Rscript -e 'source("macros/parseCensusCSV.R")' -e "parseCensusCSV(inputFile=\"data/Census2016Geomapwtlabels.csv\", outputFile=\"results/CensusData-${outputLabel}.tsv\", yCatLabel=\"${category}\", yCat0Values=${catValues}, yCatExclude=${catExclude})"

echo 
echo "====> Start data analysis program"
./build/analyzeCensusData --output "results/CensusAnalysis-${outputLabel}.root" --category "${category}" "${inputMoop}" "results/CensusData-${outputLabel}.tsv"

