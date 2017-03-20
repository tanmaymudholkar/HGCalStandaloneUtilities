#!/bin/bash

make compareHistograms &&
    ./compareHistograms.out Minbias_14TeV_9950Events_Distributions.root from_PD_5000Events_Distributions_1.root from_PD_5000Events_Distributions_2.root pythiaStandaloneTest_softQCDNonDiffractive.root
    # ./compareHistograms.out from_PD_5000Events_Distributions_1.root from_PD_5000Events_Distributions_2.root
    # ./compareHistograms.out Minbias_14TeV_9950Events_Distributions.root from_PD_5000Events_Distributions_1.root from_PD_5000Events_Distributions_2.root
    # ./compareHistograms.out pythiaStandaloneTest_softQCDNonDiffractive.root
