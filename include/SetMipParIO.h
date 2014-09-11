#ifndef SetMipParIO_h
#define SetMipParIO_h

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
void SetMipPar(Double_t MIPsPar[][672],char *M_filename,Long64_t M_time[2]){
  char RootFile[80]="./Cal/MIPs.root";
  TString RF=(TString)RootFile;
  Double_t fWidth[672];
  Double_t fMp[672];
  Double_t fArea[672];
  Double_t fGSigma[672];
  Long64_t timestart=0;
  Long64_t timestop=0;
  char FileName[80];
  memset(FileName,0,sizeof(FileName));
  for(int i=0;i<672;i++){
    fWidth[i]=MIPsPar[0][i];
    fMp[i]=MIPsPar[1][i];
    fArea[i]=MIPsPar[2][i];
    fGSigma[i]=MIPsPar[3][i];
  }
  timestart=M_time[0];
  timestop=M_time[1];
  strcpy(FileName,M_filename);


  TFile *outFile;
  TTree *MipTree;
  if(gSystem->FindFile("./",RF,kFileExists))
    {
      outFile=new TFile(RootFile,"UPDATE");
      MipTree=(TTree*)outFile->Get("MipTree");
      TBranch *b_Width;
      TBranch *b_Mp;
      TBranch *b_Area;
      TBranch *b_GSigma;
      TBranch *b_timestart;
      TBranch *b_timestop;
      TBranch *b_FileName;
      MipTree->SetBranchAddress("fWidth",fWidth,&b_Width);
      MipTree->SetBranchAddress("fMp",fMp,&b_Mp);
      MipTree->SetBranchAddress("fArea",fArea,&b_Area);
      MipTree->SetBranchAddress("fGSigma",fGSigma,&b_GSigma);
      MipTree->SetBranchAddress("timestart",&timestart,&b_timestart);
      MipTree->SetBranchAddress("timestop",&timestop,&b_timestop);
      MipTree->SetBranchAddress("FileName",FileName,&b_FileName);
}
  else{
    outFile=new TFile(RootFile,"NEW");
    MipTree=new TTree("MipTree","MipTree");
    MipTree->Branch("fWidth",&fWidth,"fWidth[672]/D");
    MipTree->Branch("fMp",&fMp,"fMp[672]/D");
    MipTree->Branch("fArea",&fArea,"fArea[672]/D");
    MipTree->Branch("fGSigma",&fGSigma,"fGSigma[672]/D");
    MipTree->Branch("timestart",&timestart,"timestart/L");
    MipTree->Branch("timestop",&timestop,"timestop/L");
    MipTree->Branch("FileName",&FileName,"FileName/C");
  }
  MipTree->Fill();
  outFile->Write("",TObject::kOverwrite);
  outFile->Close();
  delete outFile;
}
#endif
