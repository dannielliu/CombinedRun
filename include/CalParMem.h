//#include <boost/filesystem/path.hpp>
//#include <boost/filesystem/operations.hpp>
#ifndef CalParMem_h
#define CalParMem_h

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include <TString.h>
#include "Constant.h"
#include <stdio.h>
#include <iostream>

using namespace std;
using namespace cosmictest; 

class CalPar {
public:
//for all
  Long64_t timestart;
  Long64_t timestop;
  Char_t FileName[100];
  TBranch *b_timestart;
  TBranch *b_timestop;
  TBranch *b_FileName;
   
  //for Pedestal
  Double_t PedMean[nGID];
  Double_t PedSigma[nGID];
  TBranch *b_PedMean;
  TBranch *b_PedSigma;
  //for MIPs
  Double_t MipMPV[nGbar];
  Double_t MipGsigma[nGbar];
  TBranch *b_MipMPV;
  TBranch *b_MipGsigma;
  //get 
  Long64_t runtime;
};
#endif 
