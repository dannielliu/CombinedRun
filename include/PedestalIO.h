//#include <boost/filesystem/path.hpp>
//#include <boost/filesystem/operations.hpp>
#ifndef PedestalIO_h
#define PedestalIO_h

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

void SetPedPar(Double_t PedPar[][2016],char *filename,Long64_t time[2]){
  char RootFile[80]="./Cal/Calibration.root";
  // Double_t time[2]={0,100};
  TString   RF=(TString)RootFile;
  Double_t PedMean[2016];
  Double_t PedSigma[2016];
  Double_t timestart=0;
  Double_t timestop=0;
  char FileName[80];
  memset(FileName,0,sizeof(FileName));
  for(int i=0;i<(int)(sizeof(PedPar[0])/sizeof(Double_t));i++){
   PedMean[i]=PedPar[0][i];
   PedSigma[i]=PedPar[1][i];
  }
  timestart=time[0];
  timestop=time[1];
  strcpy(FileName,filename);
  

  TFile *outFile;
  TTree *PedTree;
  //boost::filesystem::path data_dir(RootFile);
//  if(boost::filesystem::exists("nima.root"))
  if(gSystem->FindFile("./",RF,kFileExists))
  {
  outFile=new TFile(RootFile,"UPDATE");
  PedTree=(TTree*)outFile->Get("PedTree");
  TBranch *b_PedMean;
  TBranch *b_PedSigma;
  TBranch *b_timestart;
  TBranch *b_timestop;
  TBranch *b_FileName;
  PedTree->SetBranchAddress("PedMean",&PedMean,&b_PedMean);
  PedTree->SetBranchAddress("PedSigma",&PedSigma,&b_PedSigma);
  PedTree->SetBranchAddress("timestart",&timestart,&b_timestart);
  PedTree->SetBranchAddress("timestop",&timestop,&b_timestop);
  PedTree->SetBranchAddress("FileName",&FileName,&b_FileName);
  }
  else{ 
  outFile=new TFile(RootFile,"NEW");
  PedTree=new TTree("PedTree","Peddata");
  PedTree->Branch("PedMean",&PedMean,"PedMean[2016]/D");
  PedTree->Branch("PedSigma",&PedSigma,"PedSigma[2016]/D");
  PedTree->Branch("timestart",&timestart,"timestart/D");
  PedTree->Branch("timestop",&timestop,"timestop/D");
  PedTree->Branch("FileName",&FileName,"FileName/C");
  }		  //event++;
  PedTree->Fill();
 
  outFile->Write("",TObject::kOverwrite);
  outFile->Close();
  delete outFile;
}
#endif
