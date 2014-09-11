//#include <boost/filesystem/path.hpp>
//#include <boost/filesystem/operations.hpp>
#ifndef CalParIO_h
#define CalParIO_h

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
#include <CalParMem.h>
#include <string.h>
using namespace std;
using namespace cosmictest; 

class CalParIO {
public:
  TFile *dataIO;//for root IO
  TTree *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t fCurrent; //!current Tree number in a TChain
  char *IOpath;
  
  //function
  CalParIO();
  virtual ~CalParIO();
  virtual void SetPath(const char* path);//style:pedestal,mips,dyrel,att,temp
  
  virtual void InitPedPar(CalPar *calped);
  virtual void SetPedPar();
  virtual void GetPedPar(CalPar *calped);
  virtual void FinPedPar();//Rearrange the event order

  virtual void InitMipPar(CalPar *calmip);
  virtual void SetMipPar();
  virtual void GetMipPar(CalPar *calmip);
  virtual void FinMipPar();//Rearrange the event order
}; 
CalParIO::CalParIO()
{}
CalParIO::~CalParIO()
{}
void CalParIO::SetPath(const char* path){
//  memset(IOpath,0,sizeof(IOpath));
//  cout<<path<<"~~~~~~~~~~~~~"<<endl;
  IOpath=(char*)malloc(100);//
  strcpy(IOpath,path);
}
void CalParIO::InitPedPar(CalPar *calped){
//char PedRoot[80]="Pedestal.root";
  TString RootFile(IOpath);
  if(gSystem->FindFile("./",RootFile,kFileExists))
  {
  cout<<"Find the root file:"<<IOpath<<">> update Or read"<<endl;
//  dataIO=new TFile(PedRoot,"UPDATE");
  dataIO=new TFile(IOpath,"update");
  fChain=(TTree*)dataIO->Get("fChain");
  fChain->SetBranchAddress("PedMean",calped->PedMean,&calped->b_PedMean);
  fChain->SetBranchAddress("PedSigma",calped->PedSigma,&calped->b_PedSigma);
  fChain->SetBranchAddress("timestart",&calped->timestart,&calped->b_timestart);
  fChain->SetBranchAddress("timestop",&calped->timestop,&calped->b_timestop);
  fChain->SetBranchAddress("FileName",calped->FileName,&calped->b_FileName);
   }
  else{ 
  cout<<"Can not find the file:"<<IOpath<<">> new"<<endl;
  dataIO=new TFile(IOpath,"NEW");
  fChain=new TTree("fChain","Peddata");
  fChain->Branch("PedMean",calped->PedMean,"calped->PedMean[2016]/D");
  fChain->Branch("PedSigma",calped->PedSigma,"calped->PedSigma[2016]/D");
  fChain->Branch("timestart",&calped->timestart,"calped->timestart/L");
  fChain->Branch("timestop",&calped->timestop,"calped->timestop/L");
  fChain->Branch("FileName",calped->FileName,"calped->FileName[80]/B");
  }		  //event++;
} 
void CalParIO::SetPedPar(){
  fChain->Fill();
//  dataIO->Close();
} 
void CalParIO::GetPedPar(CalPar *calped){
  Int_t nentries=(Int_t)fChain->GetEntries();
//cout<<"Entries:"<<nentries<<endl;
  for(int i=0;i<nentries;i++){
    fChain->GetEntry(i);  
 // cout<<i<<endl;
    if(calped->timestop>=calped->runtime)
    break;
    }
  }
void CalParIO::FinPedPar(){
  dataIO->Write("",TObject::kOverwrite);
  dataIO->Close();
}

//for MIPs
void CalParIO::InitMipPar(CalPar *calmip){
  TString RootFile(IOpath);
  if(gSystem->FindFile("./",RootFile,kFileExists))
  {
  cout<<"Find the root file:"<<IOpath<<">> update Or read"<<endl;
  dataIO=new TFile(IOpath,"update");
  fChain=(TTree*)dataIO->Get("fChain");
  fChain->SetBranchAddress("MipMPV",calmip->MipMPV,&calmip->b_MipMPV);
  fChain->SetBranchAddress("MipGsigma",calmip->MipGsigma,&calmip->b_MipGsigma);
  fChain->SetBranchAddress("timestart",&calmip->timestart,&calmip->b_timestart);
  fChain->SetBranchAddress("timestop",&calmip->timestop,&calmip->b_timestop);
  fChain->SetBranchAddress("FileName",calmip->FileName,&calmip->b_FileName);
   }
  else{ 
  cout<<"Can not find the file:"<<IOpath<<">> new"<<endl;
  dataIO=new TFile(IOpath,"NEW");
  fChain=new TTree("fChain","Mipdata");
  fChain->Branch("MipMPV",calmip->MipMPV,"calmip->MipMPV[616]/D");
  fChain->Branch("MipGsigma",calmip->MipGsigma,"calmip->MipGsigma[616]/D");
  fChain->Branch("timestart",&calmip->timestart,"calmip->timestart/L");
  fChain->Branch("timestop",&calmip->timestop,"calmip->timestop/L");
  fChain->Branch("FileName",calmip->FileName,"calmip->FileName[80]/B");
  }		  //event++;
} 
void CalParIO::SetMipPar(){
  fChain->Fill();
//  dataIO->Close();
} 
void CalParIO::GetMipPar(CalPar *calmip){
  Int_t nentries=(Int_t)fChain->GetEntries();
//cout<<"Entries:"<<nentries<<endl;
  for(int i=nentries-1;i>0;i--){
    fChain->GetEntry(i);  
 // cout<<i<<endl;
    if(calmip->timestart<=calmip->runtime)
    break;
    }
// dataIO->Close();
  }
void CalParIO::FinMipPar(){
  dataIO->Write("",TObject::kOverwrite);
  dataIO->Close();
}
#endif
