#!/bin/bash -f

#indir="/store/user/miheejo/ExpressPhysicsPA/ExpressPhysicsPA_test1/161121_183233"
indir="/store/user/miheejo/MinBias/Pyquen_WToMuNu_test1/161122_141948"
final="ExpressPhysicsPA"

for idx in 0 
do
  /afs/cern.ch/project/eos/installation/cms/bin/eos ls $indir/000$idx/ | grep root > list$idx
  list2=$(awk -v p=$indir/000$idx/ '{print "root://eoscms//eos/cms"p$1}' list$idx)
  echo $list2
  hadd /tmp/miheejo/$final\_$idx\.root $list2 >& haddlog$idx 
done

for idx in 0 
do
  rm -f list$idx
done
