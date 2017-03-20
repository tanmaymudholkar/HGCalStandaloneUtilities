#!/bin/bash

make analyzeHepMC &&
    ./analyzeHepMC.out ~/public/research/hgcal_BHStudies/mcConversion_900pre6/CMSSW_9_0_0_pre6/src/Minbias_14TeV_9950Events.dat Minbias_14TeV_9950Events_Distributions.root &&
    ./analyzeHepMC.out /afs/cern.ch/work/t/tmudholk/public/pythia/Pythia140305_000001.dat from_PD_5000Events_Distributions_1.root &&
    ./analyzeHepMC.out /afs/cern.ch/work/t/tmudholk/public/pythia/Pythia140305_000099.dat from_PD_5000Events_Distributions_2.root &&
    ./analyzeHepMC.out /afs/cern.ch/user/t/tmudholk/public/research/pythia/pythia8223/examples/pythiaTest_softQCDNonDiffractive.dat pythiaStandaloneTest_softQCDNonDiffractive.root

# make analyzeHepMC && ./analyzeHepMC.out ~/public/research/hgcal_BHStudies/mcConversion_900pre6/CMSSW_9_0_0_pre6/src/Minbias_14TeV_9950Events.dat Minbias_14TeV_9950Events_Distributions_test.root
