#!/bin/bash

# Remove the non-org files from the 0-directory
rm -rf 0/alpha1.gz 0/gamma.gz 0/U.gz 0/pd.gz 0/p_rgh.gz 0/alpha1 0/pd 0/p_rgh 0/U

# Remove the waveProperties file
rm -f constant/waveProperties

# Different syntax between Mac and Linux in determining the OF-version number
if [ `uname` = "Darwin" ]
then
    version=`echo $WM_PROJECT_VERSION | sed -e 's/\./\'$'\n/g' -e 's/-/\'$'\n/g' | grep "[0-9]" | head -2 | tr -d '\n'`
else
    version=`echo $WM_PROJECT_VERSION | sed -e 's/\./\n/g' -e 's/-/\n/g' | grep "[0-9]" | head -2 | tr -d '\n'`
fi

# alpha1 and U have the same name in all versions
cp 0/alpha1.org 0/alpha1
cp 0/U.org 0/U

# Copying fvSchemes, fvSolution and pressure fields
if [ "$version" = "16" ]
then
    cp system/fvSolution.16 system/fvSolution
    cp system/fvSchemes.16 system/fvSchemes
    cp 0/pd.org 0/pd
elif [ "$version" = "17" ]
then
    cp system/fvSolution.17 system/fvSolution
    cp system/fvSchemes.17 system/fvSchemes
    cp 0/p_rgh.org 0/p_rgh
else
    cp system/fvSolution.21 system/fvSolution
    cp system/fvSchemes.21 system/fvSchemes
    cp 0/p_rgh.org 0/p_rgh
fi
