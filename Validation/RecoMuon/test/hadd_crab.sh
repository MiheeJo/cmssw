#!/bin/bash -f

indir="/store/user/miheejo/muonTrackValidator_pp/Pythia8_Jpsi/JpsiMM_5p02TeV_TuneCUETP8M1/muonTrackValidator_pp/161028_112946/0000/"
final="JpsiMM_pp_validation.root"

/afs/cern.ch/project/eos/installation/cms/bin/eos ls $indir | grep root > list
list2=$(awk -v p=$indir '{print "root://eoscms//eos/cms"p$1}' list)
echo $list2

hadd /tmp/miheejo/$final $list2

rm -f list
