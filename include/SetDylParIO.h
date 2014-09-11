#ifndef SetDylParIO_h
#define SetDylParIO_h

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
void SetDylPar(Double_t DyRelPar[][672],char *D_filename,Long64_t D_time[2]){
  char RootFile[80]="./Cal/DyRel.root";
  TString RF=(TString)RootFile;
  Double_t Dy5_Dy2[672];
  Double_t Dy8_Dy5[672];
  Long64_t timestart=0;
  Long64_t timestop=0;
  char FileName[80];
  memset(FileName,0,sizeof(FileName));
  for(int i=0;i<672;i++){
    Dy5_Dy2[i]=DyRelPar[0][i];
    Dy8_Dy5[i]=DyRelPar[1][i];
  }
  timestart=D_time[0];
  timestop=D_time[1];
  strcpy(FileName,D_filename);

  TFile *outFile;
  TTree *DylTree;
  if(gSystem->FindFile("./",RF,kFileExists))
    {
      outFile=new TFile(RootFile,"UPDATE");
      DylTree=(TTree*)outFile->Get("DylTree");
      TBranch *b_Dy5_Dy2;
      TBranch *b_Dy8_Dy5;
      TBranch *b_timestart;
      TBranch *b_timestop;
      TBranch *b_FileName;
      DylTree->SetBranchAddress("Dy5_Dy2",Dy5_Dy2,&b_Dy5_Dy2);
      DylTree->SetBranchAddress("Dy8_Dy5",Dy8_Dy5,&b_Dy8_Dy5);
      DylTree->SetBranchAddress("timestart",&timestart,&b_timestart);
      DylTree->SetBranchAddress("timestop",&timestop,&b_timestop);
      DylTree->SetBranchAddress("FileName",FileName,&b_FileName);

}
  else{
    outFile=new TFile(RootFile,"NEW");
    DylTree=new TTree("DylTree","DylTree");
    DylTree->Branch("Dy5_Dy2",&Dy5_Dy2,"Dy5_Dy2[672]/D");
    DylTree->Branch("Dy8_Dy5",&Dy8_Dy5,"Dy8_Dy5[672]/D");
    DylTree->Branch("timestart",&timestart,"timestart/L");
    DylTree->Branch("timestop",&timestop,"timestop/L");
    DylTree->Branch("FileName",&FileName,"FileName/C");

  }
  DylTree->Fill();
  outFile->Write("",TObject::kOverwrite);
  outFile->Close();
  delete outFile;

}
#endif
