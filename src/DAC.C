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
  const int nDyPar=3;//Dy2,5,8 mean and sigma (gaus fitted)
  Double_t DACPar[nDyPar][nGID];
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //define the TH1F for pedestal 
  Float_t Cal_DACcode[15]={160,320,480,640,800,960,1120,1280,1440,1600,1760,1920,2080,2240,2400};
  Float_t Cal_Vol[16][15];
  for(int i=0;i<15;i++){
    Cal_Vol[0][i]=0.76575*Cal_DACcode[i]+5.7013;
    Cal_Vol[1][i]=0.76288*Cal_DACcode[i]+3.28139;
    Cal_Vol[2][i]=0.76615*Cal_DACcode[i]+1.8961;
    Cal_Vol[3][i]=0.75676*Cal_DACcode[i]-0.74459;
    Cal_Vol[4][i]=0.75278*Cal_DACcode[i]-2.12987;
    Cal_Vol[5][i]=0.76327*Cal_DACcode[i]+1.44156;
    Cal_Vol[6][i]=0.76291*Cal_DACcode[i]+0.38528;
    Cal_Vol[7][i]=0.75692*Cal_DACcode[i]+2.47619;
    Cal_Vol[8][i]=0.76002*Cal_DACcode[i]+1.31169;
    Cal_Vol[9][i]=0.75968*Cal_DACcode[i]+1.41991;
    Cal_Vol[10][i]=0.75951*Cal_DACcode[i]-0.26407;
    Cal_Vol[11][i]=0.75845*Cal_DACcode[i]+3.08658;
    Cal_Vol[12][i]=0.76415*Cal_DACcode[i]+1.93074;
    Cal_Vol[13][i]=0.76422*Cal_DACcode[i]-0.58874;
    Cal_Vol[14][i]=0.76073*Cal_DACcode[i]+1.59307;
    Cal_Vol[15][i]=0.75579*Cal_DACcode[i]+3.91775;
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
  TF1 *mylinear=new TF1("mylinear","[0]*x*x+[1]*x+[2]",100,800);
  //TF1 *mygaus=new TF1("mygaus","gaus",-1000,16000);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Loop start
  cout<<"Loop Root File..."<<endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"~~~~~~~~~~~~~~~"<<nentries<<endl;
  Int_t EntryBuf=0;

      int fee[nGID];
      Long64_t ientry = LoadTree(0);
     fChain->GetEntry(0);
      for(int i=0; i<nHits;i++){
        fee[GID[i]]=FEE_ID[i]-16;
    }//single event Loop end
  //for cut last data
  Long64_t tg=time;
  
      for(int t=0;t<15;t++){
    Int_t EntNub=0;
    //single loop
    for (Long64_t jentry=EntryBuf; jentry<nentries;jentry++){
      ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);   
      if (mode!=2)
      continue;
      EntNub++;
      //cut last data
     // if(time-tg>600&&EntNub<2000){
      if(time-tg>600){
        cout<<"Cut last package & Reset Gstart:"<<endl;
        for(int iGID=0;iGID<nGID;iGID++){
        EMC_Cal[iGID]->Reset();
	FEE_G[iGID]->Reset();
      }
        cout<<"Total Events:"<<jentry<<endl;
        EntNub=0;
	t=0;
      }
      tg=time;

	
 //     if(EntNub==600){
      if(EntNub==513){
        EntryBuf=jentry;
	cout<<"the EntryBuf is "<<EntryBuf<<endl;
	cout<<"the t is "<<t<<endl;
	cout<<"the EntNub is "<<EntNub<<endl;
	break;
      }
      if(jentry%1000==0)
      cout<<"Total Events:"<<jentry<<endl;
	
      for(int i=0; i<nHits;i++){
        EMC_Cal[GID[i]]->Fill(ADC[i]);}
    }//single event Loop end
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //DAC Fit and Draw
    Double_t GausPar[3];
 //   TSpectrum *s=new TSpectrum(2,50);
 //   Int_t nfound =0;
    //printf("Found %d candidate peaks \n",nfound);
 //   Float_t *xpeaks;
    for(int iGID=0;iGID<nGID;iGID++){
//      nfound=s->Search(EMC_Cal[iGID],22,"",0.10);
    // cout<<"the Peaks number is "<<nfound<<endl;
//      xpeaks=s->GetPositionX();
//      mygaus->SetRange(xpeaks[0]-100,xpeaks[0]+100);
//      EMC_Cal[iGID]->Fit(mygaus,"RQ");
//      mygaus->GetParameters(GausPar);
//      mygaus->SetRange(GausPar[1]-3*GausPar[2],GausPar[1]+3*GausPar[2]);
//      EMC_Cal[iGID]->Fit(mygaus,"RQ");
//      mygaus->GetParameters(GausPar);
      double mean;
      mean=EMC_Cal[iGID]->GetMean();
      GausPar[1]=mean;
      FEE_G[iGID]->Fill(Cal_Vol[fee[iGID]][t],GausPar[1]);
      EMC_Cal[iGID]->Reset();
    }
  }
  gStyle->SetOptStat(000);
  TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
  string Star=EPSFile;
  string FroBracket="[";
  Star=Star+FroBracket;
  Start->Print(Star.c_str());

  Double_t par[3];

  Double_t error[nDyPar][nGID];

  Double_t Chi[nGID];
  Double_t NDF[nGID];

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
	    mylinear->SetRange(50,1300);
	    gStyle->SetOptFit(1111);
	    FEE_G[iGID]->Draw();
	    FEE_G[iGID]->SetMarkerStyle(29);
	    FEE_G[iGID]->SetMarkerSize(0.5);
	    FEE_G[iGID]->Fit(mylinear,"R");
	    mylinear->GetParameters(par);
	    //extract the DAC parameters
	    DACPar[0][iGID]=par[0];
	    DACPar[1][iGID]=par[1];
        DACPar[2][iGID]=par[2];
	   
        error[0][iGID]=mylinear->GetParError(0);
        error[1][iGID]=mylinear->GetParError(1);
        error[2][iGID]=mylinear->GetParError(2);

        Chi[iGID]=mylinear->GetChisquare();
        NDF[iGID]=mylinear->GetNDF();

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
   PedPar0<<"a*x^2+b*x+c"<<" "<<"\n";
   PedPar0<<"   "<<"a"<<"   "<<"a_error"<<"     "<<"b"<<"    "<<"b_error"<<"    "<<"c"<<"    "<<"c_error"<<"    "<<"Chi"<<"     "<<"NDF"<<""<<"\n";
  for(int ig=0;ig<nGID;ig++)
    {
      PedPar0<<DACPar[0][ig]<<"  "<<error[0][ig]<<"    "<<DACPar[1][ig]<<"    "<<error[1][ig]<<"    "<<DACPar[2][ig]<<"   "<<error[2][ig]<<"    "<<Chi[ig]<<"     "<<NDF[ig]<<"";
      PedPar0<<"\n";
    }

 //   PedPar0<<"\n";
 //   PedPar0<<"chisquare:"<<Chi<<"       "<<"NDF:"<<NDF<<"\n";

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
      PedPar1<<DACPar[0][ig]<<"  "<<DACPar[1][ig]<<" "<<DACPar[2][ig]<<"   ";
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

