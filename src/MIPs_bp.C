#define MIPs_bp_cxx
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

#include "TFile.h"
#include <TF1.h>
#include "TInterpreter.h"
#include <TStyle.h>
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "MIPs_bp.h"
using namespace std;
using namespace cosmictest;

void RootConvert::Loop(char *RootFileName){
  gROOT->ProcessLine(".L /home/zhzhy/DRAW/bes3plotstyle.C");
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //output .root .eps files
  string namesuffix="_MIPs";
  string outputfile=(string)RootFileName;
  outputfile.replace(strlen(RootFileName)-5,5,namesuffix);
  int Lname=sizeof(outputfile.c_str());
  //Lname-sizeof("./Raw2ROOT/")=Lname-11
  string out;
  out.assign(outputfile,11,Lname-11);
  string EPSFile="./MIPs/"+out+".eps";
  string RFile="./Calibration/"+out+".root";
  string RName=out+".root";
  cout<<"Output file 1 "<<EPSFile<<endl;
  cout<<"Output file 2 "<<RFile<<endl;


  //file name
  ofstream name;
  name.open("./Calibration/filenameM");
  name<<RName;
  name.close();
  ifstream name1;
  char PPar[80];
  name1.open("./Calibration/filename");
  name1.getline(PPar,80);
  name1.close();
 string PPPar="./Calibration/"+(string)PPar;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //const int nGbar=616; 14*2*22
  TH1F *EMC_Mips[nGbar];
  for (int nb=0;nb<nGbar;nb++) {
    // GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nbar+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
    Int_t ibar=nb%Nbar;
    Int_t iside=((int)(nb/22))%Nside;
    Int_t ilayer=(int)(nb/Nside/22);
    char cc[30];
    sprintf(cc,"Layer%d_Side%d_Bar%d_Dy%d_MIPs",ilayer,iside,ibar,8);
    if(iside==0||(iside==1&&(ilayer==0||ilayer==13)))
    EMC_Mips[nb]=new TH1F(cc,cc,300,5,9000);
    else
    EMC_Mips[nb]=new TH1F(cc,cc,300,5,4500);
  }
  //TF1 *MIPs=new TF1("MIPs","landau",0,3800);
  TF1 *myMIPs=new TF1("myMIPs",langaufun,0,3800,4);
  myMIPs->SetParNames("Width","MP","Area","GSigma");
  const int nLangaus=4;
  Double_t MIPsPar[nGbar][nLangaus];
  for(int ig=0;ig<nGbar;ig++) {
    memset(MIPsPar[ig],0,sizeof(MIPsPar[ig]));
  }
  gStyle->SetOptStat(111);
  gStyle->SetOptFit(0111);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //LOOP START
  cout<<"Loop root file ..."<<endl;
  CalPar calped;
  CalParIO calpedIO;
  calpedIO.CalParIO::SetPath(PPPar.c_str());
  calpedIO.CalParIO::InitPedPar(&calped); 
  
  CalPar calmip;
  CalParIO calmipIO;
  calmipIO.CalParIO::SetPath(RFile.c_str());
  calmipIO.CalParIO::InitMipPar(&calmip); 
  
  cout<<"Loop Root File..."<<endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"~~~~~~~~~~~~~~~"<<nentries<<endl;
  
  Long64_t ientry;
  Long64_t Gstart=0;
  Long64_t Gstop=0;
  Long64_t EntBuf=0;
  Long64_t timegap=7200;//10mins
  int BadTrack=0;
  int Tevent=0;
  Long64_t nb = 0;
  
  ientry = LoadTree(0);
  fChain->GetEntry(0);
  Gstart=time;
  Long64_t tg=Gstart;

  ientry = LoadTree(nentries-1);
  fChain->GetEntry(nentries-1);
  Gstop=time;

  //get pedestal
  calped.runtime=Gstart;  
  calpedIO.CalParIO::GetPedPar(&calped); 
  cout<<"Get Pedestal Parameters ("<<calped.timestart<<","<<calped.timestop<<")"<<endl;
  bool reget=true;
  if(calped.timestop<calped.runtime)
  reget=false;  


  cout<<"Gstart & Gstop time: "<<Ast2Dat(Gstart)<<"("<<Gstart<<")  to  "<<Ast2Dat(Gstop)<<"("<<Gstop<<")"<<endl;
  //Loop the mip par
  Int_t EntNubChk=0;
  Int_t Nloop=(int)((Gstop-Gstart)/timegap)+1;
  cout<<"Nloop:"<<Nloop-1<<endl;
  for(int t=0;t<Nloop;t++){ 
    ientry = LoadTree(EntBuf);
    fChain->GetEntry(EntBuf);
    calmip.timestart=time;
    EntNubChk=0;
    //A single loop for one event
    for (Long64_t jentry=EntBuf; jentry<nentries;jentry++){
      ientry = LoadTree(jentry);
      if (ientry < 0) {cout<<"ientry="<<ientry<<",break";break;}
      nb = fChain->GetEntry(jentry);   
      if(jentry%10000==0)
      cout<<"Total Events:"<<jentry<<endl;
      
      if (mode!=0&&mode!=1) //select data style
      continue;
      EntNubChk++;  //check Fill number
     
      //cut last data
      if(time-tg>100){
        Gstart=time;
        cout<<"Cut last package & Reset Gstart:"<<endl;
        cout<<"Gstart:"<<Gstart<<endl;
        calmip.timestart=time;
        for(int iGbar=0;iGbar<nGbar;iGbar++){
        EMC_Mips[iGbar]->Reset();
      }
      }
      tg=time;
      //Reset Pedestal Parameters
      if(time>calped.timestop&&reget==true){
        calped.runtime=time;  
        calpedIO.CalParIO::GetPedPar(&calped); 
        cout<<"Reset Pedestal ("<<calped.timestart<<","<<calped.timestop<<")"<<endl;
        if(calped.timestop<calped.runtime)
          reget=false;  
        }
 
      if(time>Gstart+(t+1)*timegap||jentry==nentries-1){
	EntBuf=jentry;
	calmip.timestop=time;
        if(jentry==nentries-1){t=Nloop; break; } 
	if(t<(Nloop-2)&&EntNubChk>=100000)
	  break;
      }
      //fill
      Int_t Lmax[20];
      RawTrack(&calped,Lmax);
      if(Lmax[14]==1){
        for(int i=0;i<nHits;i++)  {
          if(Lmax[Layer[i]]==Bar[i]&&Dy[i]==8&&Bar[i]!=0&&Bar[i]!=23){//both sides
            if(ADC[i]>(calped.PedMean[GID[i]]+3*calped.PedSigma[GID[i]])) {
              int ibar=(Layer[i]*Nside+Side[i])*Nbar+Bar[i]-1;
              EMC_Mips[ibar]->Fill(ADC[i]-calped.PedMean[GID[i]]);
            }
          }
        }
      }
      else{BadTrack++;
      }
      Tevent++;
    }  
    //Loop end
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MIPs Fit and Draw
//    Double_t par[3];
    Double_t mean=0;
    cout<<"Fitting..."<<endl;
    if(t<Nloop-2)
    for(int iGbar=0;iGbar<nGbar;iGbar++){
          mean=EMC_Mips[iGbar]->GetMean();
	  //MIPs->SetRange(mean*0.55,mean*1.2);
	  //MIPs->SetParameters(300,mean*0.8);
	  //EMC_Mips[iGbar]->Fit(MIPs,"R0Q");
	  //MIPs->GetParameters(par);
	  myMIPs->SetRange(mean*0.45,mean*2.0);
	  //myMIPs->SetRange(par[1]-2.8*par[2],par[1]+6.50*par[2]);
	   //myMIPs->SetParameters(par);
	   // Once again, here are the Landau * Gaussian parameters:
	   //   par[0]=Width (scale) parameter of Landau density
	   //   par[1]=Most Probable (MP, location) parameter of Landau density
	   //   par[2]=Total area (integral -inf to inf, normalization constant)
	   //   par[3]=Width (sigma) of convoluted Gaussian function
          myMIPs->SetParLimits(3,0,mean);
          myMIPs->SetParLimits(1,mean*0.6,mean*1.1);
	  myMIPs->SetParameter(1,mean*0.8);
	  myMIPs->SetParameter(3,mean*0.3);
	  myMIPs->SetParameter(0,mean*80);
	  //myMIPs->SetParameter(0,par[2]);
	  //myMIPs->SetParameter(1,par[1]);
	  //myMIPs->SetParameter(3,par[2]);
	  EMC_Mips[iGbar]->Fit(myMIPs,"RQ0");
	  myMIPs->GetParameters(MIPsPar[iGbar]);
          
          if(MIPsPar[iGbar][1]>mean||MIPsPar[iGbar][1]<0.6*mean||MIPsPar[iGbar][3]<0||MIPsPar[iGbar][3]>MIPsPar[iGbar][1]){
	    myMIPs->SetRange(mean*0.50,mean*1.8);
	    //myMIPs->SetParameter(0,80);
	    ///myMIPs->SetParameter(1,mean*0.8);
	    //myMIPs->SetParameter(3,mean*0.3);
	    EMC_Mips[iGbar]->Fit(myMIPs,"RQ0");
	    myMIPs->GetParameters(MIPsPar[iGbar]);
          }
          calmip.MipMPV[iGbar]=MIPsPar[iGbar][1];
          calmip.MipGsigma[iGbar]=MIPsPar[iGbar][3];
    }
    else{
    TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
    string Star=(string)EPSFile;
    string FroBracket="[";
    Star=Star+FroBracket;
    Start->Print(Star.c_str());

    TCanvas *EMCMIPs[Nlayer][Nside];
    for(int il=0;il<Nlayer;il++){
      char d[50];
      for(int is=0;is<Nside;is++){
        sprintf(d,"Layer%d_Side%d",il,is);
        EMCMIPs[il][is]=new TCanvas(d,d,0,0,600,800);
        EMCMIPs[il][is]->Divide(4,6);
        for(int ib=0;ib<Nbar;ib++){
          Int_t iGbar=(il*Nside+is)*Nbar+ib;//idy=2;
	  EMCMIPs[il][is]->cd(ib+1);
          
          mean=EMC_Mips[iGbar]->GetMean();
//	  MIPs->SetRange(mean*0.55,mean*1.2);
//	  MIPs->SetParameters(300,mean*0.8);
//	  EMC_Mips[iGbar]->Fit(MIPs,"R0Q");
//	  MIPs->GetParameters(par);
//	  myMIPs->SetRange(par[1]-2.8*par[2],par[1]+6.00*par[2]);
	  myMIPs->SetRange(mean*0.45,mean*2.0);
	   //myMIPs->SetParameters(par);
	   // Once again, here are the Landau * Gaussian parameters:
	   //   par[0]=Width (scale) parameter of Landau density
	   //   par[1]=Most Probable (MP, location) parameter of Landau density
	   //   par[2]=Total area (integral -inf to inf, normalization constant)
	   //   par[3]=Width (sigma) of convoluted Gaussian function
	    myMIPs->SetParameter(1,mean*0.8);
	    myMIPs->SetParameter(3,mean*0.3);
	    myMIPs->SetParameter(0,80);
//	  myMIPs->SetParameter(0,par[2]);
//	  myMIPs->SetParameter(1,par[1]);
//	  myMIPs->SetParameter(3,par[2]);
          myMIPs->SetParLimits(3,0,mean);
          myMIPs->SetParLimits(1,mean*0.6,mean*1.1);
	  myMIPs->SetLineColor(2);
	  myMIPs->SetLineWidth(2);
	  EMC_Mips[iGbar]->Draw();
	  EMC_Mips[iGbar]->Fit(myMIPs,"R");
	  myMIPs->GetParameters(MIPsPar[iGbar]);
          if(MIPsPar[iGbar][1]>mean||MIPsPar[iGbar][1]<0.4*mean||MIPsPar[iGbar][3]<0||MIPsPar[iGbar][3]>MIPsPar[iGbar][1]){
	    myMIPs->SetRange(mean*0.5,mean*1.8);
	    //myMIPs->SetParameter(0,80);
	    //myMIPs->SetParameter(1,mean*0.8);
	    //myMIPs->SetParameter(3,mean*0.3);
	    EMC_Mips[iGbar]->Fit(myMIPs,"R");
	    myMIPs->GetParameters(MIPsPar[iGbar]);
          }
           calmip.MipMPV[iGbar]=MIPsPar[iGbar][1];
           calmip.MipGsigma[iGbar]=MIPsPar[iGbar][3];
         }
         EMCMIPs[il][is]->Print(EPSFile.c_str());
       } 
    }
    TCanvas *End=new TCanvas("End","End",0,0,600,800);
    string En=(string)EPSFile;
    string BacBracket="]";
    En=En+BacBracket;
    End->Print(En.c_str());
    }
    //Set Pedesrtal Par
    //to txt
    char txtfile[60];
    sprintf(txtfile,"./Calibration/MIPs/M%d-%d",(int)calmip.timestart,(int)calmip.timestop);
    ofstream outtxt;
    outtxt.open(txtfile);
    if(!outtxt.good()){cout<<"Can not open txt file!"<<endl;}
    for(int j=0;j<nGbar;j++)
    outtxt<<calmip.MipMPV[j]<<" "<<calmip.MipGsigma[j]<<endl;
    outtxt.close();

    //to root 
    cout<<"the step time : "<<Ast2Dat(calmip.timestart)<<"("<<calmip.timestart<<")  to  "<<Ast2Dat(calmip.timestop)<<"("<<calmip.timestop<<")"<<endl;
    cout<<"Gstart & Gstop time: "<<Ast2Dat(Gstart)<<"("<<Gstart<<")  to  "<<Ast2Dat(Gstop)<<"("<<Gstop<<")"<<endl;
    cout<<"The Event Number:"<<EntNubChk<<"; Rate: "<<(Double_t)EntNubChk/(Double_t)(calmip.timestop-calmip.timestart)<<"Hz"<<endl;  
    calmipIO.CalParIO::SetMipPar();
  }
    calpedIO.dataIO->Close();
    calmipIO.CalParIO::FinMipPar();
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout<<"Total event: "<<Tevent<<"\nBad Track event: "<<BadTrack<<endl;
}//Loop end 



//define main()
string filename;
void helpinfo(){
  cout<<"Usage is ./MIPs.exe ./Raw2ROOT/<filename>\n";
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
  cout<<"\n~~~~~*****~~~~~\nTHE END"<<endl;
  return 1;
}

