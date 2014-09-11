#define DAC_cxx
#include <TStyle.h>
#include <sstream>
#include <string.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <utility>
#include <stdio.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <RootConvert.h>
#include "TFile.h"
#include <TF1.h>
#include "TInterpreter.h"
#include <TStyle.h>
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"



using namespace std;
using namespace cosmictest;
void RootConvert::Loop(char *RootFileName){
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //output .root .eps files
  //string namesuffix="_pedestal.eps";
  string namesuffix="_calibration";
  string outputfile=(string)RootFileName;
  outputfile.replace(strlen(RootFileName)-5,5,namesuffix);

  int Lname=sizeof(outputfile.c_str());
  //Lname-sizeof("./Raw2ROOT/")=Lname-11
  string out;
  out.assign(outputfile,11,Lname-11);
  outputfile=out;
  //pedestal par cout 
  string TXTFile0="./DAC/"+out+".txt";
  outputfile=out+".eps";
  string  EPSFile="./DAC/"+outputfile;

  string TXTFile1="./DAC/DAC.txt";
  cout<<"Output file 1 "<<EPSFile<<endl;
  cout<<"Output file 2 "<<TXTFile0<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //define the DAC output
  const int nDyPar=2;//Dy2,5,8 mean and sigma (gaus fitted)
  Double_t DACPar[nDyPar][nGID];
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //define the TH1F for pedestal 
  Float_t Cal_DACcode[15]={160,320,480,640,800,960,1120,1280,1440,1600,1760,1920,2080,2240,2400};
  Float_t Cal_Vol[15];
  for(int i=0;i<15;i++){
    Cal_Vol[i]=0.75651*Cal_DACcode[i]+2.64069;
  }
  TH1F *EMC_Cal[nGID];
  TH2F *FEE_G[nGID];
  for (int ng=0;ng<nGID;ng++){
    // GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nch+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
    Int_t idy=ng%Ndy;
    Int_t ibar=((int)(ng/Ndy))%Nch;
    Int_t iside=((int)(ng/Nch/Ndy))%Nside;
    Int_t ilayer=(int)(ng/Nside/Nch/Ndy);
    char cc[50];
    sprintf(cc,"Layer%d_Side%d_Bar%d_Dy%d_cal",ilayer,iside,ibar,(idy*3+2));
    char dd[50];
    sprintf(dd,"Layer%d_Side%d_Bar%d_Dy%d_FEE",ilayer,iside,ibar,(idy*3+2));
    char ee[50];
    sprintf(ee,"Layer%d_Side%d_Bar%d_Dy%d_FEE;mV;ADC channels",ilayer,iside,ibar,(idy*3+2));
    EMC_Cal[ng]=new TH1F(cc,cc,4250,-1000,16000);
    FEE_G[ng]=new TH2F(dd,ee,250,0,2500,1000,-1000,16000);
  }
  TF1 *mylinear=new TF1("mylinear","[0]*x+[1]",100,800);
  TF1 *mygaus=new TF1("mygaus","gaus",-1000,16000);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Loop start
  cout<<"Loop Root File..."<<endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"~~~~~~~~~~~~~~~"<<nentries<<endl;
  Int_t EntryBuf=0;
  for(int t=0;t<15;t++){
    Int_t EntNub=0;
    //single loop
    for (Long64_t jentry=EntryBuf; jentry<nentries;jentry++){
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);   
      if (mode!=2)
      continue;
      EntNub++;
	
      if(EntNub==513){
        EntryBuf=jentry;
	cout<<"the EntryBuf is "<<EntryBuf<<endl;
	cout<<"the t is "<<t<<endl;
	cout<<"the EntNub is "<<EntNub<<endl;
	break;
      }
      if(jentry%1000==0)
      cout<<"Total Events:"<<jentry<<endl;
	
      for(int i=0; i<nHits;i++)
        EMC_Cal[GID[i]]->Fill(ADC[i]);
    }//single event Loop end
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //DAC Fit and Draw
    Double_t GausPar[3];
    TSpectrum *s=new TSpectrum(2,50);
    Int_t nfound =0;
    //printf("Found %d candidate peaks \n",nfound);
    Float_t *xpeaks;
    for(int iGID=0;iGID<nGID;iGID++){
      nfound=s->Search(EMC_Cal[iGID],22,"",0.10);
    // cout<<"the Peaks number is "<<nfound<<endl;
      xpeaks=s->GetPositionX();
      mygaus->SetRange(xpeaks[0]-100,xpeaks[0]+100);
      EMC_Cal[iGID]->Fit(mygaus,"RQ");
      mygaus->GetParameters(GausPar);
      mygaus->SetRange(GausPar[1]-3*GausPar[2],GausPar[1]+3*GausPar[2]);
      EMC_Cal[iGID]->Fit(mygaus,"RQ");
      mygaus->GetParameters(GausPar);
      FEE_G[iGID]->Fill(Cal_Vol[t],GausPar[1]);
      EMC_Cal[iGID]->Reset();
    }
  }
  gStyle->SetOptStat(000);
  TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
  string Star=EPSFile;
  string FroBracket="[";
  Star=Star+FroBracket;
  Start->Print(Star.c_str());
  Double_t par[2];
  TCanvas *EMC_Ped[Nlayer][Nside][Ndy];
  for(int ilayer=0;ilayer<Nlayer;ilayer++)
    for(int iside=0;iside<Nside;iside++){
      char d[3][50];
      for(int idy=0;idy<Ndy;idy++)
	{
	  sprintf(d[idy],"Layer%d_Side%d_Dy%d",ilayer+1,iside,idy*3+2);
	  EMC_Ped[ilayer][iside][idy]=new TCanvas(d[idy],d[idy],0,0,600,800);
	  EMC_Ped[ilayer][iside][idy]->Divide(4,6);
	  for(int ibar=0;ibar<Nch;ibar++){
	    Int_t iGID=((ilayer*Nside+iside)*Nch+ibar)*Ndy+idy;
	    EMC_Ped[ilayer][iside][idy]->cd(ibar+1);
	    mylinear->SetRange(150,1500);
	    gStyle->SetOptFit(1111);
	    FEE_G[iGID]->Draw();
	    FEE_G[iGID]->SetMarkerStyle(29);
	    FEE_G[iGID]->SetMarkerSize(0.5);
	    FEE_G[iGID]->Fit(mylinear,"R");
	    mylinear->GetParameters(par);
	    //extract the DAC parameters
	    DACPar[0][iGID]=par[0];
	    DACPar[1][iGID]=par[1];
	   
	  }
	  EMC_Ped[ilayer][iside][idy]->Print(EPSFile.c_str());	  
	}
    }
  TCanvas *End=new TCanvas("End","End",0,0,600,800);
  string En=EPSFile;
  string BacBracket="]";
  En=En+BacBracket;
  End->Print(En.c_str());
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //cout DACPar
  ofstream PedPar0;
  PedPar0.open(TXTFile0.c_str());
  if (!PedPar0.good())
    {cout<<"Can not open DAC TXTFlie Name File!!!"<<endl;
      exit(0);
    }
  //Double_t DACPar[nTBar][nDyPar];
  for(int ig=0;ig<nGID;ig++)
    {
      PedPar0<<DACPar[0][ig]<<"  "<<DACPar[1][ig]<<" ";
      PedPar0<<"\n";
    }
  PedPar0.close();
 
  //cout DACPar
  ofstream PedPar1;
  PedPar1.open(TXTFile1.c_str());
  if (!PedPar1.good())
    {cout<<"Can not open DAC TXTFlie Name File!!!"<<endl;
      exit(0);
    }
  //Double_t DACPar[nTBar][nDyPar];
  for(int ig=0;ig<nGID;ig++)
    {
      PedPar1<<DACPar[0][ig]<<"  "<<DACPar[1][ig]<<" ";
      PedPar1<<"\n";
    }
  PedPar1.close();

  //filename
  ofstream name;

  name.open("./DAC/filename");
  if(!name.good())
    {cout<<"Can not open namefile !"<<endl;exit(1);}
  name<<RootFileName;
  name.close();

  cout<<"\n~~~~~*****~~~~~\nTHE END"<<endl;
}//Loop end 




//define main()
string filename;
void helpinfo(){
  cout<<"Usage is ./DAC.exe ./Raw2ROOT/<filename>\n";
  cout<<"default output file name is <inputfilename.eps>"<<endl;	
  return;
}
void phrase_command(int argc,char *argv[]){
  if (argc<2){ helpinfo();
    exit(0);
  }else {filename=(string) argv[1]; cout<<"START\n~~~~~*****~~~~~\nInput File : "<<filename<<endl;}
}

int main(int argc,char *argv[])
{

  phrase_command(argc,argv);
  string Tt="CTest";
  //string RootFileName=filename;
  const char *TtreeName=(char *)(Tt.data());
  const Text_t* RawRootFile=(Text_t*)(filename.data());

  char *RawRootFile_c=(char *)(filename.data());
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(RawRootFile);
  if (!f) 
    {
      f = new TFile(RawRootFile_c);
    }
  TTree *tree = (TTree*)gDirectory->Get(TtreeName);
 
  RootConvert xx;
  xx.RootConvert::Init(tree);
  xx.RootConvert::Loop(RawRootFile_c);
  return 1;

}

