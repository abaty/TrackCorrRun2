#!/bin/bash
if [ $# -ne 1 ]; then
    echo $0: usage: script isPP
    exit 1
fi

isPP=""
if [ "$1" = "pp" ]; then
  isPP=$1
fi
inDir="/mnt/hadoop/cms/store/user/abaty/temporaryStorage/"
inFile="*Hists_job*.root"

echo copying histograms over
echo make sure to copy the TrkCorrInputFile.txt and weights file over as well!
mv $inDir$isPP$inFile . 
