#include <iostream>
#include <cstdlib>
#include <map>

#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "Rtypes.h"
#include "TColor.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLegendEntry.h"

std::map<Int_t, Int_t> colorsForComparisonPlots = {{0, kBlack}, {1, kBlue}, {2, kRed}, {3, kGreen}};
// std::map<Int_t, Int_t> colorsForComparisonPlots = {{0, kBlue}, {1, kRed}};
// std::map<Int_t, Int_t> colorsForComparisonPlots = {{0, kGreen}};
std::map<Int_t, std::string> namesForComparisonPlots = {{0, "Minbias_14TeV"}, {1, "from_PD1"}, {2, "from_PD2"}, {3, "pythiaStandalone"}};
// std::map<Int_t, std::string> namesForComparisonPlots = {{0, "from_PD1"}, {1, "from_PD2"}};
// std::map<Int_t, std::string> namesForComparisonPlots = {{0, "pythiaStandalone"}};
const std::string outputFileName = "comparison_HepMC.root";
const std::string outputPlotsFolder = "plots_all4";
// const std::string outputPlotsFolder = "plots_pythiaStandalone";

void normalizeHistograms(std::vector<TH1F*> inputHistograms) {
  for (unsigned histogramCounter = 0; histogramCounter < inputHistograms.size(); ++histogramCounter) {
    double histIntegral = inputHistograms[histogramCounter]->Integral("width");
    inputHistograms[histogramCounter]->Scale(1./histIntegral);
  }
}

std::pair<double, double> getCombinedYRangeProfiles(std::vector<TProfile*> inputHistograms) {
  double ymin = inputHistograms[0]->GetMinimum();
  double ymax = inputHistograms[0]->GetMaximum();

  for (unsigned histogramsToCompareCounter = 1; histogramsToCompareCounter < inputHistograms.size(); ++histogramsToCompareCounter) {
    double histymin = inputHistograms[histogramsToCompareCounter]->GetMinimum();
    double histymax = inputHistograms[histogramsToCompareCounter]->GetMaximum();
    if (histymin < ymin) ymin = histymin;
    if (histymax > ymax) ymax = histymax;
  }

  std::pair<double, double> combinedYRange;
  if (ymin > 0) {
    combinedYRange.first = 0.;
  }
  else {
    combinedYRange.first = ymin - 0.1*(ymax - ymin);
  }
  combinedYRange.second = ymax + 0.1*(ymax - ymin);
  return combinedYRange;
}

std::pair<double, double> getCombinedYRange(std::vector<TH1F*> inputHistograms) {
  double ymin = inputHistograms[0]->GetMinimum();
  double ymax = inputHistograms[0]->GetMaximum();

  for (unsigned histogramsToCompareCounter = 1; histogramsToCompareCounter < inputHistograms.size(); ++histogramsToCompareCounter) {
    double histymin = inputHistograms[histogramsToCompareCounter]->GetMinimum();
    double histymax = inputHistograms[histogramsToCompareCounter]->GetMaximum();
    if (histymin < ymin) ymin = histymin;
    if (histymax > ymax) ymax = histymax;
  }

  std::pair<double, double> combinedYRange;
  if (ymin > 0) {
    combinedYRange.first = 0.;
  }
  else {
    combinedYRange.first = ymin - 0.1*(ymax - ymin);
  }
  combinedYRange.second = ymax + 0.1*(ymax - ymin);
  return combinedYRange;
}

void createProfilePlots(std::vector<TFile*> inputFiles, TFile *outputFile, std::string nameToFetch) {
  std::vector<TProfile*> profilesToPlot;
  for (unsigned fileCounter = 0; fileCounter < inputFiles.size(); ++fileCounter) {
    TProfile *profileToPlot = nullptr;
    TFile *inputFile = inputFiles[fileCounter];
    inputFile->GetObject(nameToFetch.c_str(), profileToPlot);
    if (profileToPlot) {
      profilesToPlot.push_back(profileToPlot);
      std::cout << "File counter: " << fileCounter << ": opened profile " << nameToFetch << std::endl;
    }
    else {
      std::cout << "File counter: " << fileCounter << ": unable to open profile " << nameToFetch << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  std::pair<double, double> combinedYRange = getCombinedYRangeProfiles(profilesToPlot);
  
  TCanvas *outputCanvas = new TCanvas(("c_" + nameToFetch).c_str(), ("c_" + nameToFetch).c_str(), 1280, 1024);
  TLegend *legend = new TLegend(0.1, 0.85, 0.25, 1.0, "Colors:");
  TProfile *baseProfile = profilesToPlot[0];
  TLegendEntry *baseLegendEntry = legend->AddEntry(baseProfile, namesForComparisonPlots[0].c_str(), "");
  baseLegendEntry->SetTextColor(colorsForComparisonPlots[0]);
  baseProfile->SetLineColor(colorsForComparisonPlots[0]);
  baseProfile->GetYaxis()->SetRangeUser(combinedYRange.first, combinedYRange.second);
  baseProfile->Draw();

  for (unsigned profilesToCompareCounter = 1; profilesToCompareCounter < profilesToPlot.size(); ++profilesToCompareCounter) {
    TProfile *profileToCompare = profilesToPlot[profilesToCompareCounter];
    TLegendEntry *legendEntry = legend->AddEntry(profileToCompare, namesForComparisonPlots[profilesToCompareCounter].c_str(), "");
    legendEntry->SetTextColor(colorsForComparisonPlots[profilesToCompareCounter]);
    profileToCompare->SetLineColor(colorsForComparisonPlots[profilesToCompareCounter]);
    profileToCompare->Draw("SAME");
  }
  legend->SetTextSize(0.02);
  legend->Draw();
  outputFile->WriteTObject(outputCanvas);
  outputCanvas->SaveAs((outputPlotsFolder + "/" + nameToFetch + ".png").c_str());
  delete outputCanvas;
}

void createPlots(std::vector<TFile*> inputFiles, TFile *outputFile, std::string nameToFetch, bool normalizationRequired) {
  std::vector<TH1F*> histogramsToPlot;
  for (unsigned fileCounter = 0; fileCounter < inputFiles.size(); ++fileCounter) {
    TH1F *histogramToPlot = nullptr;
    TFile *inputFile = inputFiles[fileCounter];
    inputFile->GetObject(nameToFetch.c_str(), histogramToPlot);
    if (histogramToPlot) {
      histogramsToPlot.push_back(histogramToPlot);
      std::cout << "File counter: " << fileCounter << ": opened histogram " << nameToFetch << std::endl;
    }
    else {
      std::cout << "File counter: " << fileCounter << ": unable to open histogram " << nameToFetch << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (normalizationRequired) normalizeHistograms(histogramsToPlot);
  std::pair<double, double> combinedYRange = getCombinedYRange(histogramsToPlot);
  
  TCanvas *outputCanvas = new TCanvas(("c_" + nameToFetch).c_str(), ("c_" + nameToFetch).c_str(), 1280, 1024);
  TLegend *legend = new TLegend(0.1, 0.85, 0.25, 1.0, "Colors:");
  TH1F *baseHistogram = histogramsToPlot[0];
  TLegendEntry *baseLegendEntry = legend->AddEntry(baseHistogram, namesForComparisonPlots[0].c_str(), "");
  baseLegendEntry->SetTextColor(colorsForComparisonPlots[0]);
  baseHistogram->SetLineColor(colorsForComparisonPlots[0]);
  baseHistogram->GetYaxis()->SetRangeUser(combinedYRange.first, combinedYRange.second);
  baseHistogram->Draw();

  for (unsigned histogramsToCompareCounter = 1; histogramsToCompareCounter < histogramsToPlot.size(); ++histogramsToCompareCounter) {
    TH1F *histogramToCompare = histogramsToPlot[histogramsToCompareCounter];
    TLegendEntry *legendEntry = legend->AddEntry(histogramToCompare, namesForComparisonPlots[histogramsToCompareCounter].c_str(), "");
    legendEntry->SetTextColor(colorsForComparisonPlots[histogramsToCompareCounter]);
    histogramToCompare->SetLineColor(colorsForComparisonPlots[histogramsToCompareCounter]);
    histogramToCompare->Draw("SAME");
  }
  legend->SetTextSize(0.02);
  legend->Draw();
  outputFile->WriteTObject(outputCanvas);
  outputCanvas->SaveAs((outputPlotsFolder + "/" + nameToFetch + ".png").c_str());
  delete outputCanvas;
}

int main(int argc, char ** argv) {
  if (argc < 2) {
    std::cout << "Please enter name of base input root file" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  gStyle->SetOptStat(0);

  std::cout << "Starting program..." << std::endl;

  std::vector<std::string> inputFileNames;
  unsigned nFilesToCompare = static_cast<unsigned>(argc - 1);
  
  for (unsigned inputFileNameCounter = 1; inputFileNameCounter <= nFilesToCompare; ++inputFileNameCounter) {
    inputFileNames.push_back(std::string(argv[inputFileNameCounter]));
  }

  std::vector<TFile*> inputFiles;
  for (unsigned inputFileNameCounter = 0; inputFileNameCounter < nFilesToCompare; ++inputFileNameCounter) {
    TFile *inputFile = TFile::Open(inputFileNames[inputFileNameCounter].c_str());
    inputFiles.push_back(inputFile);
  }
  
  TFile *outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");
  TFile *baseFile = inputFiles[0];
  TList *listOfKeys = baseFile->GetListOfKeys();
  TListIter keyIterator(listOfKeys);
  TObject *componentObject = nullptr;
  while ( (componentObject = (TObject*)keyIterator.Next()) ) {
    std::string nameToFetch = std::string(componentObject->GetName());
    std::cout << "Comparing object with name: " << nameToFetch << std::endl;
    bool plotProfileInsteadOfHistogram = (nameToFetch == "hEvtMultVsAvgpT");
    if (plotProfileInsteadOfHistogram) {
      createProfilePlots(inputFiles, outputFile, nameToFetch);
    }
    else {
      if (nameToFetch == "hdNdpTVspT") {
        createPlots(inputFiles, outputFile, nameToFetch, false);
      }
      else {
        createPlots(inputFiles, outputFile, nameToFetch, true);
      }
    }
  }
  outputFile->Close();
}
