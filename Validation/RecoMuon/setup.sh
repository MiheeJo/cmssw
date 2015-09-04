#!/bin/bash
$ git cms-addpkg CondFormats/HIObjects
$ git cms-addpkg DataFormats/HeavyIonEvent
$ git cms-addpkg HeavyIonsAnalysis/Configuration
$ git cms-addpkg HiAnalysis/HiTagAndProbe
$ git cms-addpkg RecoHI/HiCentralityAlgos
$ git cms-addpkg Validation/RecoMuon

$ scram b

# Running script:
$ cmsRun Validation/RecoMuon/test/muonValidation.py 
