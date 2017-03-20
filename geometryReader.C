#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TH2F.h"
#include "TH2Poly.h"
// #include "TAxis.h"
#include "TStyle.h"
#include "TMath.h"
// #include "TSystem.h"
// #include "Rtypes.h"
#include "TList.h"
#include "TCollection.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TCollection.h"
// #include "TMath.h"
#include "Math/Point3D.h"
#define MAX_FRACTIONAL_TOLERANCE 0.001
#define STARTLAYER 36
#define CONST_PI 3.14159265

Double_t minRVal = -1.;
Double_t maxRVal = 2500.;

const std::string inputLayerZPositionsFileName = "/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/PFCal/PFCalEE/analysis/macros/MinBias/layerZPositions.dat";
Double_t *layerZPositions;
unsigned nLayers = 52;

template <typename T>
void printVector(std::vector<T> inputVector) {
  for (std::vector<Double_t>::iterator inputVectorIterator = inputVector.begin(); inputVectorIterator != inputVector.end(); ++inputVectorIterator) {
    std::cout << "Vector element: " << *inputVectorIterator << std::endl;
  }
}

void printAxisBinning(TAxis *inputAxis) {
  unsigned nBins = inputAxis->GetNbins();
  for (unsigned binCounter = 0; binCounter <= 1+nBins; ++binCounter) {
    std::cout << "Bin number " << binCounter << " ranges from " << inputAxis->GetBinLowEdge(binCounter) << " to " << inputAxis->GetBinUpEdge(binCounter) << " and its center is " << inputAxis->GetBinCenter(binCounter) << std::endl;
  }
}

void fillLayerZPositions() {
  std::vector<Double_t> zPositions;
  std::string line;
  std::ifstream inputFile(inputLayerZPositionsFileName.c_str());
  if (inputFile.is_open()) {
    while(std::getline(inputFile, line)) {
      zPositions.push_back(atof(line.c_str()));
    }
    inputFile.close();
  }
  else {
    std::cout << "ERROR: Unable to open input file to read in layer z positions" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  layerZPositions = new Double_t[zPositions.size()];
  for (unsigned layerCounter = 0; layerCounter < zPositions.size(); ++layerCounter) {
    layerZPositions[layerCounter] = zPositions[layerCounter];
    // std::cout << "z position [" << layerCounter << "] = " << layerZPositions[layerCounter] << std::endl;
  }
  if (zPositions.size() != nLayers) {
    std::cout << "Something has gone terribly wrong... number of layers present in input file for z positions of layers is not the same as the number of layers obtained from initializing detector with this design" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  zPositions.clear();
}

// void PolarTH2()
// {
//   gStyle->SetOptStat(0);
//   Double_t rmax(1.);

//   TH2D* pol_his = new TH2D("polarHist", "polarHist", 20, 0., 2.*TMath::Pi(), 20, 0., rmax); // the polar data. X maps to theta, Y maps to R

//   Double_t theta, r;
//   //fill the histogram with something nice.
//   for(Int_t i=1; i<=pol_his->GetNbinsX(); i++)
//     {
//       theta = pol_his->GetXaxis()->GetBinCenter(i);
//       for(Int_t j=1; j<=pol_his->GetNbinsY(); j++)
// 	{
//           r = pol_his->GetYaxis()->GetBinCenter(j);
//           pol_his->SetBinContent(i, j, r*cos(2.*theta) );
// 	}
//     }

//   //make a dummy histogram to set the axis range from -rmax to rmax in both directions
//   TH2D* dummy_his = new TH2D("dummy", "histo title", 100, -rmax, rmax, 100, -rmax, rmax);

//   TCanvas* c1 = new TCanvas("theCanvas", "theCanvas", 600, 600);
//   dummy_his->Draw("COL"); // draw the dummy histogram first
//   pol_his->Draw("COL POL SAME"); // now draw the data histogram. If it has "SAME" it will use the first histogram ranges
// }

// void readInGeometryFromFile(std::string fileName, std::vector<Double_t> &geometryVector, unsigned layerCounter) {
//   std::cout << "Starting to read in file with name " << fileName << std::endl;
//   geometryVector.clear();
//   std::ifstream inputFile;
//   inputFile.open(fileName);
//   if (inputFile.is_open()) {
//     // std::cout << "Opened file! Now reading in geometry..." << std::endl;
//     std::string layerIdentifier;
//     Double_t centerR, area, outerR, innerR;
//     Double_t outerR_old = 0;
//     bool isFirstLine = true;
//     while (inputFile >> layerIdentifier >> centerR >> area >> outerR >> innerR) {
//       // std::cout << "layerIdentifier: " << layerIdentifier << "; centerR = " << centerR << "; area = " << area << "; outerR = " << outerR << "; innerR = " << innerR << std::endl;
//       geometryVector.push_back(innerR);
//       if (isFirstLine) {
//         isFirstLine = false;
//         Double_t innermostEta = (ROOT::Math::XYZPoint(0, innerR, layerZPositions[layerCounter]).eta());
//         std::cout << "Innermost eta = " << innermostEta << std::endl;
//       }
//       else if((fabs(outerR_old - innerR)/innerR) > MAX_FRACTIONAL_TOLERANCE) {
//         std::cout << "ERROR: old outer R, " << outerR << ", does not match new inner R, " << innerR << std::endl;
//         // std::exit(EXIT_FAILURE);
//       }
//       outerR_old = outerR;
//     }
//     geometryVector.push_back(outerR);
//     Double_t outermostEta = (ROOT::Math::XYZPoint(0, outerR, layerZPositions[layerCounter]).eta());
//     std::cout << "Outermost eta = " << outermostEta << std::endl;
//   }
//   else {
//     std::cout << "Could not open file!" << std::endl;
//   }
//   inputFile.close();
// }
std::vector<Double_t> readInGeometryFromFile(std::string inputFileName) {
  std::cout << "Starting to read in file with name " << inputFileName << std::endl;
  std::vector<Double_t> geometryVector;
  std::ifstream inputFile;
  inputFile.open(inputFileName);
  if (inputFile.is_open()) {
    std::cout << "Opened file! Now reading in geometry..." << std::endl;
    // std::string layerIdentifier;
    Int_t layerNumber;
    Double_t centerR, area, mipsig, sipm_noise, sig_over_noise, power, fluence, dose, outerR, innerR;
    Double_t outerR_old = 0;
    bool isFirstLine = true;
    // while (inputFile >> layerIdentifier >> centerR >> area >> outerR >> innerR) {
    while (inputFile >> layerNumber >> centerR >> area >> mipsig >> sipm_noise >> sig_over_noise >> power >> fluence >> dose >> outerR >> innerR) {
      // std::cout << "layerIdentifier: " << layerIdentifier << "; centerR = " << centerR << "; area = " << area << "; outerR = " << outerR << "; innerR = " << innerR << std::endl;
      std::cout << "layerNumber: " << layerNumber
                << "; centerR = " << centerR
                << "; area = " << area
                << "; mipsig = " << mipsig
                << "; sipm_noise = " << sipm_noise
                << "; sig_over_noise = " << sig_over_noise
                << "; power = " << power
                << "; fluence = " << fluence
                << "; dose = " << dose
                << "; outerR = " << outerR
                << "; innerR = " << innerR << std::endl;
      
      geometryVector.push_back(innerR);
      if (isFirstLine) {
        isFirstLine = false;
      }
      else if((fabs(outerR_old - innerR)/innerR) > MAX_FRACTIONAL_TOLERANCE) {
        std::cout << "ERROR: old outer R, " << outerR << ", does not match new inner R, " << innerR << std::endl;
        std::exit(EXIT_FAILURE);
      }
      outerR_old = outerR;
    }
    geometryVector.push_back(outerR);
  }
  else {
    std::cout << "Could not open file " << inputFileName << " for reading!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  inputFile.close();
  std::cout << "Successfully read geometry!" << std::endl;
  return geometryVector;
}

void setLayerBinning(TH2F *h_layerGeometry, std::vector<Double_t> geometryVector, Int_t nPhiBins) {
  std::vector<Double_t> phiBinEdges;
  for (Int_t phiBinCounter = 0; phiBinCounter <= nPhiBins; ++phiBinCounter) {
    phiBinEdges.push_back(2*TMath::Pi()*(phiBinCounter*1./nPhiBins));
  }
  h_layerGeometry->SetBins(nPhiBins, &phiBinEdges[0], -1+geometryVector.size(), &geometryVector[0]);
}

// void setLayerBinningPoly(TH2Poly *h_testPoly, std::vector<Double_t> geometryVector, Int_t nPhiBins) {
//   std::cout << "Setting poly layer binning..." << std::endl;
//   for (Int_t phiBinCounter = 0; phiBinCounter < nPhiBins; ++phiBinCounter) {
//     // phiBinEdges.push_back(2*TMath::Pi()*(phiBinCounter*1./nPhiBins));
//     Double_t phiBoundaryLow = 2*TMath::Pi()*(phiBinCounter*1./nPhiBins);
//     Double_t phiBoundaryHigh = 2*TMath::Pi()*((phiBinCounter+1)*1./nPhiBins);
//     for (unsigned geometryVectorCounter = 0; geometryVectorCounter < -1+geometryVector.size(); ++geometryVectorCounter) {
//       Double_t radialBoundaryLow = geometryVector[geometryVectorCounter];
//       Double_t radialBoundaryHigh = geometryVector[geometryVectorCounter+1];
//       h_testPoly->AddBin(phiBoundaryLow, radialBoundaryLow, phiBoundaryHigh, radialBoundaryHigh);
//       if (minRVal == -1.) minRVal = radialBoundaryLow;
//       else if (radialBoundaryLow < minRVal) minRVal = radialBoundaryLow;
//       if (maxRVal == -1.) maxRVal = radialBoundaryHigh;
//       else if (radialBoundaryHigh > maxRVal) maxRVal = radialBoundaryHigh;
//     }
//   }
//   // setLayerBinning((TH2F*)h_testPoly, geometryVector, nPhiBins);
//   //Code dustbin:
//   // std::cout << "Initialized phi bin edges" << std::endl;
//   // h_testPoly->SetBins(nPhiBins, &phiBinEdges[0], -1+geometryVector.size(), &geometryVector[0]);
// }

void drawXValuesPoly(TH2Poly *inputTH2Poly) {
  std::cout << "Drawing x values..." << std::endl;
  TIter next(inputTH2Poly->GetBins());
  std::cout << "Variable next initialized" << std::endl;
  TObject *obj=0;
  // Int_t nBinsTot = inputTH2Poly->GetNumberOfBins();
  TH2PolyBin *polyBin = 0;
    
  while ((obj=next())){
    // for (Int_t binCounter = 1; binCounter <= nBinsTot; ++binCounter) {
    polyBin=(TH2PolyBin*)obj;
    
    int id = polyBin->GetBinNumber();
    Double_t phi = (polyBin->GetXMax()+polyBin->GetXMin())/2.;
    Double_t r = (polyBin->GetYMax()+polyBin->GetYMin())/2.;
    inputTH2Poly->SetBinContent(id, r*cos(phi));
    std::cout << "id = " << id << "; r = " << r << "; phi = " << phi << std::endl;
  }

  // std::cout << "Filled th2poly" << std::endl;
  TCanvas *c_output = new TCanvas("c_output", "c_output", 1024, 768);
  // TH2F *dummyHist = new TH2F("dummyHist", "dummyHist", 100, -1.1*maxRVal, 1.1*maxRVal, 100, -1.1*maxRVal, 1.1*maxRVal);
  // dummyHist->Draw("COL");
  // inputTH2Poly->Draw("COLZ POL SAME");
  // inputTH2Poly->GetXaxis()->SetRangeUser(0, 2*TMath::Pi());
  // inputTH2Poly->GetYaxis()->SetRangeUser(minRVal, maxRVal);
  inputTH2Poly->Draw("COLZ");
  c_output->SaveAs("test_xValues_poly.png");
  delete c_output;
}

void drawYValuesPoly(TH2Poly *inputTH2Poly) {
  TIter next(inputTH2Poly->GetBins());
  TObject *obj=0;
  TH2PolyBin *polyBin = 0;
    
  while ((obj=next())){
    polyBin=(TH2PolyBin*)obj;
    int id = polyBin->GetBinNumber();
    Double_t phi = (polyBin->GetXMax()+polyBin->GetXMin())/2.;
    Double_t r = (polyBin->GetYMax()+polyBin->GetYMin())/2.;
    inputTH2Poly->SetBinContent(id, r*sin(phi));
  }
  TCanvas *c_output = new TCanvas("c_output", "c_output", 1024, 768);
  // c_output->SetGrid();
  // inputTH2Poly->GetXaxis()->SetRangeUser(0, 2*TMath::Pi());
  // inputTH2Poly->GetYaxis()->SetRangeUser(minRVal, maxRVal);
  // TH2F *dummyHist = new TH2F("dummyHist", "dummyHist", 100, -1.1*maxRVal, 1.1*maxRVal, 100, -1.1*maxRVal, 1.1*maxRVal);
  // dummyHist->Draw("COL");
  // inputTH2Poly->Draw("COLZ POL SAME");
  inputTH2Poly->Draw("COLZ");
  c_output->SaveAs("test_yValues_poly.png");
  delete c_output;
}

void drawIDValuesPoly(TH2Poly *inputTH2Poly) {
  TIter next(inputTH2Poly->GetBins());
  TObject *obj=0;
  TH2PolyBin *polyBin = 0;
    
  while ((obj=next())){
    polyBin=(TH2PolyBin*)obj;
    int id = polyBin->GetBinNumber();
    // Double_t phi = (polyBin->GetXMax()+polyBin->GetXMin())/2.;
    // Double_t r = (polyBin->GetYMax()+polyBin->GetYMin())/2.;
    inputTH2Poly->SetBinContent(id, id*1.0);
  }
  TCanvas *c_output = new TCanvas("c_output", "c_output", 1024, 768);
  // inputTH2Poly->GetXaxis()->SetRangeUser(0, 2*TMath::Pi());
  // inputTH2Poly->GetYaxis()->SetRangeUser(minRVal, maxRVal);
  // c_output->SetGrid();
  // TH2F *dummyHist = new TH2F("dummyHist", "dummyHist", 100, -1.1*maxRVal, 1.1*maxRVal, 100, -1.1*maxRVal, 1.1*maxRVal);
  // dummyHist->Draw("COL");
  // inputTH2Poly->Draw("COLZ POL SAME");
  inputTH2Poly->Draw("COLZ");
  c_output->SaveAs("test_idValues_poly.png");
  delete c_output;
}

void drawXValues(TH2F *inputHistogram, std::string layerSuffix) {
  Int_t nBinsX = inputHistogram->GetXaxis()->GetNbins();
  Int_t nBinsY = inputHistogram->GetYaxis()->GetNbins();
  for (Int_t binXCounter = 1; binXCounter <= nBinsX; ++binXCounter) {
    Double_t phi = inputHistogram->GetXaxis()->GetBinCenter(binXCounter);
    for (Int_t binYCounter = 1; binYCounter <= nBinsY; ++binYCounter) {
      Double_t r = inputHistogram->GetYaxis()->GetBinCenter(binYCounter);
      inputHistogram->SetBinContent(binXCounter, binYCounter, r*cos(phi));
    }
  }
  TCanvas *c_output = new TCanvas("c_output", "c_output", 1024, 768);
  TH2F *dummyHist = new TH2F("X values", "X values", 100, -1.1*maxRVal, 1.1*maxRVal, 100, -1.1*maxRVal, 1.1*maxRVal);
  dummyHist->Draw("COL");
  inputHistogram->Draw("COLZ POL SAME");
  c_output->SaveAs(("test_xValues_" + layerSuffix + ".png").c_str());
  delete c_output;
}

void drawYValues(TH2F *inputHistogram) {
  Int_t nBinsX = inputHistogram->GetXaxis()->GetNbins();
  Int_t nBinsY = inputHistogram->GetYaxis()->GetNbins();
  for (Int_t binXCounter = 1; binXCounter <= nBinsX; ++binXCounter) {
    Double_t phi = inputHistogram->GetXaxis()->GetBinCenter(binXCounter);
    // std::cout << "Switching to mean phi = " << phi << std::endl;
    for (Int_t binYCounter = 1; binYCounter <= nBinsY; ++binYCounter) {
      Double_t r = inputHistogram->GetYaxis()->GetBinCenter(binYCounter);
      // std::cout << "Switching to mean r = " << r << std::endl;
      inputHistogram->SetBinContent(binXCounter, binYCounter, r*sin(phi));
    }
  }
  TCanvas *c_output = new TCanvas("c_output", "c_output", 1024, 768);
  // c_output->SetGrid();
  TH2F *dummyHist = new TH2F("Y values", "Y values", 100, -1.1*maxRVal, 1.1*maxRVal, 100, -1.1*maxRVal, 1.1*maxRVal);
  dummyHist->Draw("COL");
  inputHistogram->Draw("COLZ POL SAME");
  c_output->SaveAs("test_yValues.png");
  delete c_output;  
}

double rad(double x, double y) {
  return sqrt(x*x + y*y);
}

double phi(double x, double y) {
  // if (x > 0 && y > 0) return ((180.*atan(y/x))/TMath::Pi());
  // if (x < 0 && y > 0) return ((180.*atan(y/x))/TMath::Pi());
  double result = (ROOT::Math::XYZPoint(x, y, 0.).phi());
  if (result < 0) result += 2*TMath::Pi();
  return result;
}

void geometryReader(std::string layerSuffix, unsigned layerNumber) {
  gStyle->SetOptStat(0);

  // std::cout << "Reading in layer Z positions..." << std::endl;
  // fillLayerZPositions();
  
  std::cout << "Starting to read in data for layer number..." << layerNumber << std::endl;
  std::vector<Double_t> testGeometryVector;
  testGeometryVector = readInGeometryFromFile("FH_BH_Geometry/geometry_" + layerSuffix + ".txt");
  // testGeometryVector = readInGeometryFromFile("FH_BH_Geometry/geometry_" + layerSuffix + ".txt", layerNumber);
  TH2F *h_layerGeometry = new TH2F("h_testGeometry", "h_testGeometry", 1, 0., 0., 1, 0., 0.);
  if (layerNumber >= 40) setLayerBinning(h_layerGeometry, testGeometryVector, 288);
  else setLayerBinning(h_layerGeometry, testGeometryVector, 360);
  std::cout << "Checking... X axis binning:" << std::endl;
  printAxisBinning(h_layerGeometry->GetXaxis());
  std::cout << "Checking... Y axis binning:" << std::endl;
  printAxisBinning(h_layerGeometry->GetYaxis());
  printVector(testGeometryVector);
  // TH2Poly *h_testPoly = new TH2Poly();
  // h_testPoly->SetNameTitle("h_testGeometry", "h_testGeometry");
  // std::cout << "Initialized poly" << std::endl;
  // setLayerBinningPoly(h_testPoly, testGeometryVector, 288);

  drawXValues(h_layerGeometry, layerSuffix);
  // drawYValues(h_layerGeometry);
  // // TH2Poly *h_testPoly = (TH2Poly*)h_layerGeometry;

  // std::cout << "Now checking poly..." << std::endl;
  // // TH2Poly *h_testPoly = new TH2Poly("h_testGeometry", "h_testGeometry", 1, 0., 0., 1, 0., 0.);
  
  // // std::cout << "Layer binning successfully set" << std::endl;

  // h_testPoly->GetXaxis()->SetRangeUser(0, 2*TMath::Pi());
  // h_testPoly->GetYaxis()->SetRangeUser(minRVal, maxRVal);
  // drawXValuesPoly(h_testPoly);
  // drawYValuesPoly(h_testPoly);
  // drawIDValuesPoly(h_testPoly);
  // // TH2F *testClosure = (TH2F*)testPoly;
  // // drawXValues(testClosure);
  // // drawYValues(testClosure);

  // std::cout << "Testing...10-1 is:" << std::endl;
  // unsigned ten = 10;
  // unsigned one = 1;
  // std::cout << std::to_string(ten-one) << std::endl;
  // std::cout << "Checking... sqrt(5) = " << sqrt(5) << std::endl;
  // std::cout << "Checking... for (x,y) = (1,0), r = " << rad(1,0) << ", phi = " << phi(1,0) << std::endl;
  // std::cout << "Checking... for (x,y) = (1,1), r = " << rad(1,1) << ", phi = " << phi(1,1) << std::endl;
  // std::cout << "Checking... for (x,y) = (0,1), r = " << rad(0,1) << ", phi = " << phi(0,1) << std::endl;
  // std::cout << "Checking... for (x,y) = (-1,1), r = " << rad(-1,1) << ", phi = " << phi(-1,1) << std::endl;
  // std::cout << "Checking... for (x,y) = (-1,0), r = " << rad(-1,0) << ", phi = " << phi(-1,0) << std::endl;
  // std::cout << "Checking... for (x,y) = (-1,-1), r = " << rad(-1,-1) << ", phi = " << phi(-1,-1) << std::endl;
  // std::cout << "Checking... for (x,y) = (0,-1), r = " << rad(0,-1) << ", phi = " << phi(0,-1) << std::endl;
  // std::cout << "Checking... for (x,y) = (1,-1), r = " << rad(1,-1) << ", phi = " << phi(1,-1) << std::endl;
}

#ifndef __CINT__
int main(int argc, char* argv[]) {
  if (argc != 3) std::cout << "Please enter layer suffix and layer number!" << std::endl;
  geometryReader(std::string(argv[1]), static_cast<unsigned>(std::atoi(argv[2])));
  return 0;
}
#endif
