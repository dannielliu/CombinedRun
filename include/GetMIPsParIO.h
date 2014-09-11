#ifndef GetMIPsParIO_h
#define GetMIPsParIO_h

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

void GetMIPsPar(Double_t MIPsPar[][672],char *M_filename,Long64_t M_time[2],Long64_t Miptime){
  TFile *f=new TFile("./Cal/MIPs.root");
  TTree *ktree=(TTree*)f->Get("MipTree");
  Long64_t stimestart=0;
  Long64_t stimestop=0;
  Double_t sWidth[672];
  Double_t sMp[672];
  Double_t sArea[672];
  Double_t sGSigma[672];
  char sFileName[80];
  memset(sFileName,0,sizeof(sFileName));
  Long64_t timegap;
  Long64_t Mintime=60000000000;
  Int_t ib=0;

  TBranch *kWidth=ktree->GetBranch("fWidth");
 TBranch *kMp=ktree->GetBranch("fMp")
 TBranch *kArea=ktree->GetBranch("fArea");
 TBranch *kGSigma=ktree->GetBranch("fGSigma");
 TBranch *ktimestart=ktree->GetBranch("timestart");
 TBranch *ktimestop=ktree->GetBranch("timestop");
 TBranch *kFileName=ktree->GetBranch("FileName");

 ktimestart->SetAddress(&stimestart);
 ktimestop->SetAddress(&stimestop);
 kWidth->SetAddress(&sWidth);
 kMp->SetAddress(&sMp);
 kArea->SetAddress(&sArea);
 kGSigma->SetAddress(&sGSigma);
 kFileName->SetAddress(&sFileName);

 Int_t nentries=(Int_t)ktree->GetEntries();
 for(int i=0;i<nentries;i++){
   ktree->GetEntry(i);
   timegap=Miptime-stimestart;
   if(timegap>=0){
     if(timegap<=Mintime){
       Mintime=timegap;
       ib=i;
     }
   }
 }
 ktree->GetEntry(ib);
 strcpy(M_filename,sFileName);
 M_time[0]=stimestart;
 M_time[1]=stimestop;
 for(Int_t j=0;j<672;j++){
   MIPsPar[0][j]=sWidth[j];
   MIPsPar[1][j]=sMp[j];
   MIPsPar[2][j]=sArea[j];
   MIPsPar[3][j]=sGSigma[j];
 }
 f->Close();
}
#endlif
