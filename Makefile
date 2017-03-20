geometryReader: geometryReader.C
	g++ -Wall -Wextra `root-config --cflags --glibs` -o geometryReader.out geometryReader.C

analyzeHepMC: analyzeHepMC.C
	g++ -Wall -Wextra `root-config --cflags --glibs` -L/afs/cern.ch/user/t/tmudholk/public/research/pythia/hepmc/x86_64-slc6-gcc46-opt/lib -lHepMC -o analyzeHepMC.out analyzeHepMC.C

compareHistograms: compareHistograms.C
	g++ -std=c++0x -Wall -Wextra `root-config --cflags --glibs` -o compareHistograms.out compareHistograms.C

clean:
	rm *.out
