if [ $# -ne 1 ]
then
  echo "Usage: ./run.sh <condor_iteration>"
  exit 1
fi

sleep $1

echo | awk -v i=$1 -v '{print "./calcCorr.exe "i" "0}' 
echo | awk -v i=$1 -v '{print "./calcCorr.exe "i" "0}' | bash

echo | awk -v user=$USER '{print "mv *Hists*.root /mnt/hadoop/cms/store/user/"user"/temporaryStorage/"}'
echo | awk -v user=$USER '{print "mv *Hists*.root /mnt/hadoop/cms/store/user/"user"/temporaryStorage/"}' | bash
rm *.root
echo "job done successfully"
