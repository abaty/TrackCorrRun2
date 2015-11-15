if [ $# -ne 0 ]
then
  echo "Usage: ./psort.sh <trackqual> <file-list> <tag> <nmin> <nmax> <pttrigmin> <pttrigmax> <ptassmin> <ptassmax>"
  exit 1
fi

now="Corrections_$(date +"%Y_%m_%d__%H_%M_%S")"
njobs=1

mkdir $now
cp ../TrkCorrInputFile.txt $now
cp ../src/calcCorr.C $now
cp ../src/Settings.h $now
cp ../src/prepStep.C $now
cp ../src/makeSkim.C $now
cp ../src/getWeights.C $now
cp ../src/iterate.C $now

cp run.sh $now

#insert blacklist code
blacklist=""
for i in $(cat /net/hisrv0001/home/abaty/condor_blacklist/condor_blacklist.txt); do
  
  blacklist=$blacklist$i" \&\& "
done
blacklist=$blacklist"endoflinetag"
blacklist=$(echo $blacklist | sed "s@ \\\&\\\& endoflinetag@@g")
echo "blacklist: "$blacklist 
cat run.condor | sed "s@blacklist_here@$blacklist@g" > $now/run.condor

cat $now/TrkCorrInputFile.txt | sed "s@SEDTARGETDATE@$(date +"%Y_%m_%d__%H_%M_%S")@g" > $now/TrkCorrInputFile.txt
cat $now/run.condor | sed "s@SEDTARGETDATE@$(date +"%Y_%m_%d__%H_%M_%S")@g" > $now/run.condor

cat $now/run.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@user_flag@$USER@g" |  sed "s@arglist@@g" | sed "s@njobs@$njobs@g" > $now/run.condor
sleep 1
cd $now
g++ prepStep.C $(root-config --cflags --libs) -Wall -O2 -o "prepStep.exe"
g++ calcCorr.C $(root-config --cflags --libs) -Wall -O2 -o "calcCorr.exe"
#g++ prepStep.C $(root-config --cflags --libs) -Werror -Wall -O2 -o "prepStep.exe"
#g++ calcCorr.C $(root-config --cflags --libs) -Werror -Wall -O2 -o "calcCorr.exe"
echo finished compilation
echo
sleep 1
cat run.condor
echo 
echo running prepStep in 10s
sleep 10
./prepStep.exe
#echo condor_submit $now/run.condor
#echo
# condor_submit $now/run.condor

