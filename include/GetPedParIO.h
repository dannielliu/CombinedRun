//void GetPedPar(Long64_t Miptime){
#ifndef GetPedParIO_h
#define GetPedParIO_h

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"


void GetPedPar(Double_t zPedPar[][2016],char *zFileName,Long64_t ztimess[2],Long64_t Miptime){  
  //TFile *kfile=new TFile("kfile","RECREATE");
 TFile *f=new TFile("./Cal/Pedestal.root");
 TTree *ktree=(TTree*)f->Get("PedTree");     //????????????????

  Long64_t ftimestart=50;
  Long64_t ftimestop=50;
  Double_t fPedMean[2016];
  Double_t fPedSigma[2016];
  // char *fFileName=NULL;

  char fFileName[80];
  memset(fFileName,0,sizeof(fFileName));
  Long64_t timegap;
  //Long64_t Miptime;
  Long64_t Mintime=60000000000;
  Int_t ib=0;
  
  // Long64_t Maxtime=-6000000;
  // Long64_t Midval=0;
  TBranch *kPedMean=ktree->GetBranch("PedMean");
  TBranch *kPedSigma=ktree->GetBranch("PedSigma");
  TBranch *ktimestart=ktree->GetBranch("timestart");
  TBranch *ktimestop=ktree->GetBranch("timestop");
  TBranch *kFileName=ktree->GetBranch("FileName");
  //PedTree->SetBranchAddress("ktimestart",&timestart);
  //PedTree->SetBranchAddress("ktimestop",&ktimestop);
  ktimestart->SetAddress(&ftimestart);
  ktimestop->SetAddress(&ftimestop);
  kPedMean->SetAddress(fPedMean);
  kPedSigma->SetAddress(fPedSigma);
  kFileName->SetAddress(fFileName);

  Int_t nentries=(Int_t)ktree->GetEntries();
  // ktree->GetEntry(0);
  for(int i=0;i<nentries;i++){
    ktree->GetEntry(i);
     timegap=Miptime-ftimestart;
    // if(mistime[nentries]>0){
    if (timegap>=0){    
      if (timegap<=Mintime){
	Mintime=timegap;
	ib=i;
      }    
    }
    //}
    //else{
    // if(mistime[nentries]>Maxtime)
    //	mistime[nentries]=mistime[netries];
    // } 
  }
  ktree->GetEntry(ib);
  //zFileName=(char *)malloc(strlen(fFileName)+1);
  strcpy(zFileName,fFileName);
  //cout<<*zFileName<<endl;
  ztimess[0]=ftimestart;
  ztimess[1]=ftimestop;
  for (Int_t j=0;j<2016;j++){
    zPedPar[0][j]=fPedMean[j];
    zPedPar[1][j]=fPedSigma[j];
  }
 

    //kPedMean[2016]
    //kPedSigma[2016]->GetEntry(ib);
  //kfile->cd();
  //kfile->Write();
  // kfile->Close();
  f->Close();
}
#endif
