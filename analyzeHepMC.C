#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"

using namespace HepMC;

bool isCharged(int pdgId) {
  int unsignedpdgId = (pdgId > 0? pdgId : -pdgId);
  if (unsignedpdgId == 21 || unsignedpdgId == 22) {
    return false;
  }
  return true;
}

int main(int argc, char ** argv) {

  if (argc != 3) {
    std::cout << "Please enter names of input and output files" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string inputFile = argv[1];
  std::string outputFile = argv[2];

  std::cout << " -- Starting program..." << std::endl;

  // const unsigned nZ = 100;
  // const unsigned z0 = 100;
  // const unsigned step = 2;

  //TFile *fout = TFile::Open("PLOTS/output_hepmc_vtxorig_MB.root","RECREATE");
  //TFile *fout = TFile::Open("PLOTS/output_hepmc_vtxmodif.root","RECREATE");
  TFile *fout = TFile::Open(outputFile.c_str(),"RECREATE");
  fout->cd();
  // TH1F *hvtx_x = new TH1F("hvtx_x",";x (mm);vertices",200,-0.06,0.06);
  // TH1F *hvtx_y = new TH1F("hvtx_y",";y (mm);vertices",200,-0.03,0.03);
  TH1F *hvtx_x = new TH1F("hvtx_x",";x (mm);vertices",200,0.,0.);
  TH1F *hvtx_y = new TH1F("hvtx_y",";y (mm);vertices",200,0.,0.);
  TH1F *hvtx_z = new TH1F("hvtx_z",";z (mm);vertices",200,-200,200);
  TH1F *hvtx_t = new TH1F("hvtx_t",";t (ns);vertices",5000,-250,250);
  // TH1F *hTotNParticles = new TH1F("hTotNParticles", ";NParticles;events", 1000, -0.5, 1000.5);
  TH1F *hTotNParticles = new TH1F("hTotNParticles", ";NParticles;frac_events", 250, -0.5, 250.5);
  TProfile *hEvtMultVsAvgpT = new TProfile("hEvtMultVsAvgpT", ";pT_avg(GeV);evtMultiplicity", 54, 0.1, 0.65);
  TH1F *hEta = new TH1F("hEta", ";eta;frac_finalStateParticles", 250, -5., 5.);
  TH1F *hEtaAbs = new TH1F("hEtaAbs", ";etaAbs;frac_finalStateParticles", 250, 0., 5.);
  TH1F *hEnergyWeightedEta = new TH1F("hEnergyWeightedEta", ";eta;frac_finalStateParticles", 250, -5., 5.);
  TH1F *hEnergyWeightedEtaAbs = new TH1F("hEnergyWeightedEtaAbs", ";etaAbs;frac_finalStateParticles", 250, 0., 5.);
  TH1F *hp = new TH1F("hp", ";p(GeV);frac_finalStateParticles", 200, 0., 400.);
  TH1F *hp_avg = new TH1F("hp_avg", ";p_avg(GeV);frac_events", 200, 0., 400.);
  TH1F *hpT = new TH1F("hpT", ";pT(GeV);frac_finalStateParticles", 100, 0., 0.);
  TH1F *hpT_avg = new TH1F("hpT_avg", ";pT_avg(GeV);frac_events", 100, 0., 0.);
  TH1F *hdNdpTVspT = new TH1F("hdNdpTVspT", ";pT(GeV);#frac{dN}{dp_{T}}", 199, 0., 2.);
  // TH1F *hpT_gt30_avg = new TH1F("hpT_gt30_avg", ";pT_gt30_avg(GeV);frac_events", 100, 0., 0.);
  TH1F *hChargeMultiplicity = new TH1F("hChargeMultiplicity", ";chargeMultiplicity;frac_finalStateParticles", 500, -0.5, 500.5);
  // TH2F *hvtx_tvsz = new TH2F("hvtx_tvsz",";z (mm); t (ps);vertices",1000,-200,200,1200,-600,600);
  // TH1F *hvtx_t_z[nZ];
  // for (unsigned i(0);i<nZ;++i){
  //   std::ostringstream label;
  //   int zmin = -1*z0+i*step;
  //   int zmax = -1*z0+(i+1)*step;
  //   label << "hvtx_t_z_" << zmin << "_" << zmax;
  //   hvtx_t_z[i] = new TH1F(label.str().c_str(),";t (ps);vertices",120,-600,600);
  // }
  TH1F *hProton = new TH1F("hProton",";p (GeV);particles",1000,0,100);
  TH1F *hNeutron = new TH1F("hNeutron",";p (GeV);particles",1000,0,100);
  TH1F *hPipm = new TH1F("hPipm",";p (GeV);particles",1000,0,100);
  TH1F *hProtonLog = new TH1F("hProtonLog",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *hNeutronLog = new TH1F("hNeutronLog",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *hPipmLog = new TH1F("hPipmLog",";log(p) (log(GeV));particles",100,-2,2);

  // const unsigned nFiles = 1;//6;

  // for (unsigned i(0);i<nFiles;++i){

  // specify an input file
  // std::ostringstream lname;
  //lname << "/afs/cern.ch/work/p/pdauncey/public/Pythia140305_";
  //if (i==0) lname << "000000";
  //else if (i<10) lname << "00000" << i;
  //else if (i<100) lname << "0000" << i;
  //else if (i<1000) lname << "000" << i;
  //else if (i<10000) lname << "00" << i;
  //else if (i<100000) lname << "0" << i;
  //else lname << i;
  //lname << ".dat";
  //lname << "/afs/cern.ch/work/a/amagnan/public/HepMCFiles/ggHgg_origVtx.dat";
  //lname << "/afs/cern.ch/work/a/amagnan/public/HepMCFiles/ggHgg_modifyVtx.dat";
  //lname << "/afs/cern.ch/work/a/amagnan/public/HepMCFiles/vertexHLLHC.dat";
  //lname << "/afs/cern.ch/work/p/pdauncey/public/Pythia140305_000000.dat";
  // lname << "/afs/cern.ch/work/a/amagnan/public/HepMCFiles/ggHgg_1428658356.dat";
  // lname << "/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/mcConversion_900pre6/CMSSW_9_0_0_pre6/src/Minbias_14TeV_9950Events.dat";
  // HepMC::IO_GenEvent ascii_in(lname.str().c_str(),std::ios::in);
  HepMC::IO_GenEvent ascii_in(inputFile.c_str(),std::ios::in);

  int totNEvents;
    
  // get the first event
  HepMC::GenEvent* evt = ascii_in.read_next_event();
  // loop until we run out of events
  while ( evt ) {
    unsigned ievt =  evt->event_number();
    if (ievt%100==0) std::cout << "Processing Event Number "
				 << ievt
				 << std::endl;


    // GenVertex * parent = 0;
    HepMC::GenEvent::vertex_const_iterator q = evt->vertices_begin();

    //for (; q != evt->vertices_end(); ++q ){

    //if ((*q)->position().x()!=0 || (*q)->position().y()!=0) continue;

    double z = (*q)->position().z();
    double t = (*q)->position().t();
    hvtx_x->Fill((*q)->position().x());
    hvtx_y->Fill((*q)->position().y());
    hvtx_z->Fill(z);
    hvtx_t->Fill(t);
    //hvtx_tvsz->Fill(z,t*1000);
    // hvtx_tvsz->Fill(z,t);

    // if (fabs(z)<z0){
    //   unsigned idx = static_cast<unsigned>((z+z0)*1./step);
    //   if (idx>(nZ-1)) continue;
    //   hvtx_t_z[idx]->Fill(t*1000);
    // }

    /*std::cout << " -- vtx pos: " << (*q)->position().x() << " " << (*q)->position().y() << " " << (*q)->position().z() << " nParticles: in=" << (*q)->particles_in_size() << " " << (*q)->particles_out_size()
      << std::endl;
      for ( HepMC::GenVertex::particles_in_const_iterator p
      = (*q)->particles_in_const_begin(); p != (*q)->particles_in_const_end(); ++p ){

      std::cout << " ---- in particle " << (*p)->pdg_id() << " status " << (*p)->status()
      << std::endl;
	  
      }
      for ( HepMC::GenVertex::particles_out_const_iterator p
      = (*q)->particles_out_const_begin(); p != (*q)->particles_out_const_end(); ++p ){

      std::cout << " ---- out particle " << (*p)->pdg_id() << " status " << (*p)->status()
      << std::endl;
	  
      }*/
    //}

    // hTotNParticles->Fill(evt->particles_size());
    
    // analyze the event
    HepMC::GenEvent::particle_const_iterator lPart = evt->particles_begin();
    //unsigned counter = 0;
      
    //std::cout << " -- Number of particles: " << evt->particles_size()
    //<< std::endl;

    
    double pT_avg = 0.;
    // double pT_gt30_avg = 0.;
    double p_avg = 0.;
    int nParticlesFinalState = 0;
    // int nParticlesPTgt30 = 0;
    int chargeMultiplicity = 0;
    for (; lPart!=evt->particles_end();++lPart){
      //std::cout << counter << " " 
      //<< (*lPart)->pdg_id() <<  " " 
      //<< (*lPart)->status() 
      //<< std::endl;
      //counter++;
      if ((*lPart)->status()!=1) continue;
      if ((*lPart)->momentum().eta() < 0) continue;
      int pdgId = (*lPart)->pdg_id();
      double pZ = (*lPart)->momentum().pz();
      double pX = (*lPart)->momentum().px();
      double pY = (*lPart)->momentum().py();
      // double p = sqrt(pow((*lPart)->momentum().px(),2)+pow((*lPart)->momentum().py(),2)+pow((*lPart)->momentum().pz(),2));
      double p = sqrt(pX*pX + pY*pY + pZ*pZ);
      //if (fabs((*lPart)->momentum().eta())>2.8 && fabs((*lPart)->momentum().eta())<3.0){
      // double pT = sqrt(pow((*lPart)->momentum().px(),2)+pow((*lPart)->momentum().py(),2));
      double pT = sqrt(pX*pX + pY*pY);
      // double absTheta = atan(pT/p);
      // double absEta = -log(tan(0.5*absTheta));
      // bool etaPositive = (pZ >= 0);
      // if (!(etaPositive)) std::cout << "Aha! eta is not positive for at least some pZ: " << pZ << std::endl;
      // else std::cout << "Bummer! eta is positive for this pZ: " << pZ << std::endl;
      // double absEta = atanh(pZ/p);
      double eta = atanh(pZ/p);
      // double eta = ( etaPositive ? absEta : -absEta);
      // std::cout << "Here, eta = " << eta << std::endl;
      // if (p>10.) std::cout << "pT = " << pT << "; p = " << p << std::endl;
      nParticlesFinalState += 1;
      pT_avg += pT;
      p_avg += p;
      // if (pT > 30.) {
      //   nParticlesPTgt30 += 1;
      // }
      hEta->Fill((*lPart)->momentum().eta());
      hEtaAbs->Fill(fabs(eta));
      hEnergyWeightedEta->Fill((*lPart)->momentum().eta(), p);
      hEnergyWeightedEtaAbs->Fill(fabs(eta), p);
      hpT->Fill(pT);
      hp->Fill(p);
      hdNdpTVspT->Fill(pT);
      if (fabs((*lPart)->momentum().eta())<2.5){
        if (fabs(pdgId)==2212) {hProton->Fill(p);hProtonLog->Fill(log10(p));}
        else if (fabs(pdgId)==2112) {hNeutron->Fill(p);hNeutronLog->Fill(log10(p));}
        else if (fabs(pdgId)==211) {hPipm->Fill(p);hPipmLog->Fill(log10(p));}
      }
      if (isCharged(pdgId)) {
        chargeMultiplicity += 1;
      }
      else {
        if (pdgId != 22) std::cout << "Found non-photonic neutral particle with pdg id " << pdgId << " in the final state." << std::endl;
      }
    }
    if (nParticlesFinalState > 0) {
      hTotNParticles->Fill(nParticlesFinalState);
      pT_avg = pT_avg/nParticlesFinalState;
      p_avg = p_avg/nParticlesFinalState;
      hChargeMultiplicity->Fill(chargeMultiplicity);
      hpT_avg->Fill(pT_avg);
      hp_avg->Fill(p_avg);
      hEvtMultVsAvgpT->Fill(pT_avg, 1.0*nParticlesFinalState);
    }
    // if (nParticlesPTgt30 > 0) {
    //   pT_gt30_avg = pT_gt30_avg/nParticlesPTgt30;
    //   hpT_gt30_avg->Fill(pT_gt30_avg);
    // }

    // delete the created event from memory
    delete evt;
    // read the next event
    ascii_in >> evt;
    ievt++;
    totNEvents++;
  }
  // }

  hdNdpTVspT->Scale(1./totNEvents);

  fout->Write();
  
  return 0;

}
