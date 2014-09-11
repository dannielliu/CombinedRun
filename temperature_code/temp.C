int temp(){
char time[20];
char header[4][3];
char temp_flag[3];
char FEE_ID[3];
char line[200];
char tail[3];
Long64_t timeStart;
Long64_t timeN;
Int_t Fee=0;
Int_t resistance[4];
Int_t current[2];
Double_t Te[4];
Double_t Cu[2];
//Int_t FeeC[8]={16,17,18,19,27,28,29,30};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//read map
//S_No, Quadrant FEE_ID chn Layer Bar Position;
ifstream map;
map.open("rst_map.txt");
if(!map.good())
{cout<<"can not open map file!";
exit(1);
}
char mapheader[80];
map.getline(mapheader,80);
//cout<<mapheader<<endl;
//cout<<"~~~~~~~~~~~"<<endl;
Int_t RPM[64][7];//Resistance Position Map;
Int_t Layer[16][4];
Int_t Bar[16][4];
Int_t Position[16][4]; 
for(int i=0;i<64;i++){
for(int j=0;j<7;j++){
map>>RPM[i][j];
}
int l=RPM[i][2]-16;
int c=RPM[i][3]-1;
Layer[l][c]=RPM[i][4];
Bar[l][c]=RPM[i][5];
Position[l][c]=RPM[i][6];
//cout<<" "<<Layer[l][c]<<" "<<Bar[l][c]<<" "<<Position[l][c]<<endl;
}
map.close();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH3F *Tm=new TH3F("Temperature","TemperatureField;X(bars);Y(bars);Z(layers)",22,0,22,22,0,22,14,0,14);
char F_time[20]="2014-03-11-23-57-00";
bool F_flag[16];
for(int f=0;f<16;f++)
F_flag[f]=false;
//
TH2F *Rsc[16][4];
for(int i=0;i<16;i++)
for(int j=0;j<4;j++)
{
char cc[20];
sprintf(cc,"FEE%d_Rsc%d",i+16,j+1);
char dd[20];
sprintf(dd,";time(h);temp",i+16);
Rsc[i][j]=new TH2F(cc,dd,3000,0,2,300,-25,50);
}
    
TH2F *crt[16][2];
for(int i=0;i<16;i++){
	for(int j=0;j<2;j++){
    		char ee[20];
		char ff[20];
		if(j==0){sprintf(ee,"FEE%d",i+16);}
		else {sprintf(ee,"FEE%d_Crt%d",i+16,j);}
		//sprintf(ee,"FEE%d_Crt%d",i+16,j);
    		crt[i][j]=new TH2F(ee,ee,3000,0,25,300,0,400);
	}
}

void ccopy(char* destination, char* source,int start,int num){
for(int i=start;i<start+num;i++)
{destination[i-start]=source[i];}
}

int c2i(char source){
int aa=(int)source;
if (aa>=48&&aa<=57)
aa=aa-48;
if(aa>=97&&aa<=102)
aa=aa-87;
return aa;
}

Long64_t time_c2i(char* source){
Long64_t res=(c2i(source[8])*10+c2i(source[9]))*3600*24+(c2i(source[11])*10+c2i(source[12]))*3600+(c2i(source[14])*10+c2i(source[15]))*60+c2i(source[17])*10+c2i(source[18]);
return res;
}

Double_t TC(Int_t N,char flag){//flag='F':FEE ,flag='B':BGO
Double_t B=0;
Double_t R0=0;
Double_t F=0;
if(flag=='F'){
B=4100.;
R0=17.25;
F=0.9524;
}
if(flag=='B'){
B=3650.;
R0=15;
F=0.9524;
}
Double_t NN=(Double_t)N;
Double_t Tem= 273.15/(1+273.15/B*TMath::Log(10*NN/(4096-F*NN)/R0))-273.15;
return Tem;
}
    
Double_t CC(Int_t C) //Current Converter
{
    Double_t I=(2000*C/4096)-182;
    return I;
}

ifstream readtemp;
readtemp.open("CMD_LOG.txt");
if(!readtemp.good()){
cout<<"Can not open temp file!!"<<endl;
exit(1);
}
readtemp.getline(line,200);
//timeStart=(c2i(line[8])*10+c2i(line[9]))*3600*24+(c2i(line[11])*10+c2i(line[12]))*3600+(c2i(line[14])*10+c2i(line[15]))*60+c2i(line[17])*10+c2i(line[18]);
timeStart=time_c2i(line);
char start[20];
memset(start,0,sizeof(start));
ccopy(start,line,0,19);
cout<<"Start!"<<endl;
cout<<"Time:"<<start<<endl;
cout<<"......"<<endl;
while(!readtemp.eof()){
readtemp.getline(line,200);
if(line[0]=='2')
ccopy(time,line,0,19);
//cout<<time<<endl;
ccopy(header[0],line,21,2);
if(header[0][0]=='5'&&header[0][1]=='5'){
  ccopy(header[1],line,24,2);
  if(header[1][0]=='a'&&header[1][1]=='a'){
    ccopy(header[2],line,27,2);
    if(header[2][0]=='e'&&header[2][1]=='b'){
      ccopy(header[3],line,30,2);
      if(header[3][0]=='9'&&header[3][1]=='0'){
        ccopy(FEE_ID,line,33,2);
        ccopy(temp_flag,line,36,2);
        if(temp_flag[0]=='0'&&temp_flag[1]=='6'){
          ccopy(tail,line,72,2);
          if(tail[0]=='a'&&tail[1]=='5'){
        
            // timeN=(c2i(line[8])*10+c2i(line[9]))*3600*24+(c2i(line[11])*10+c2i(line[12]))*3600+(c2i(line[14])*10+c2i(line[15]))*60+c2i(line[17])*10+c2i(line[18]);
            timeN=time_c2i(line);
            timeN=timeN-timeStart;
            double timeH=(double)timeN/3600;
            Fee=c2i(line[33])*16+c2i(line[34]);
            resistance[0]=(c2i(line[39])*16+c2i(line[40]))*16*16+c2i(line[42])*16+c2i(line[43]);
            resistance[1]=(c2i(line[45])*16+c2i(line[46]))*16*16+c2i(line[48])*16+c2i(line[49]);
            resistance[2]=(c2i(line[51])*16+c2i(line[52]))*16*16+c2i(line[54])*16+c2i(line[55]);
            resistance[3]=(c2i(line[57])*16+c2i(line[58]))*16*16+c2i(line[60])*16+c2i(line[61]);
            Te[0]=TC(resistance[0],'B');
            Te[1]=TC(resistance[1],'B');
            Te[2]=TC(resistance[2],'B');
            Te[3]=TC(resistance[3],'F');

       //   cout<<timeN<<" "<<Fee<<" ";
       //   for(int i=0;i<4;i++)
       //   cout<<resistance[i]<<" ";
       //   cout<<"\n";
            //fill
            int f=Fee%16;
            for(int n=0;n<4;n++)
//          Rsc[f][n]->Fill(timeH,resistance[n]);
            Rsc[f][n]->Fill(timeH,Te[3-n]);//data order is BGO1,BGO2,BGO3,FEE
            //fill temperature field 3D plot
              if(timeN>(time_c2i(F_time)-timeStart)&&F_flag[f]==false){
                for(int n=0;n<4;n++){
                  memset(time,0,sizeof(time));
                  ccopy(time,line,0,19);
                  cout<<time<<":"<<Layer[f][n]<<" "<<Bar[f][n]<<" "<<Position[f][n]<<endl;
                  if(Layer[f][n]*Bar[f][n]!=0){
                  if(Layer[f][n]%2==0&&Fee>=28)
                  Tm->Fill(Bar[f][n],Position[f][n]/60*22,Layer[f][n],Te[3-n]/20);
                  else if(Layer[f][n]%2==0&&Fee>=24)
                  Tm->Fill(Bar[f][n],22-Position[f][n]/60*22,Layer[f][n],Te[3-n]/20);
                  else if(Layer[f][n]%2==1&&Fee>=20)
                  Tm->Fill(Position[f][n]/60*22,Bar[f][n],Layer[f][n],Te[3-n]);
                  else if(Layer[f][n]%2==1&&Fee>=16)
                  Tm->Fill(22-Position[f][n]/60*22,Bar[f][n],Layer[f][n],Te[3-n]);
                   }
                  }
                F_flag[f]=true;
                }
            }
          }
        else if(temp_flag[0]=='0'&&temp_flag[1]=='5'){
            ccopy(tail,line,60,2);
            if(tail[0]=='a'&&tail[1]=='5'){
                timeN=time_c2i(line);
                timeN=timeN-timeStart;
                double timeH=(double)timeN/3600;
		Fee=c2i(line[33])*16+c2i(line[34]);
		int f=Fee%16;		
                current[0]=(c2i(line[39])*16+c2i(line[40]))*16*16+c2i(line[42])*16+c2i(line[43]);
                current[1]=(c2i(line[45])*16+c2i(line[46]))*16*16+c2i(line[48])*16+c2i(line[49]);
                Cu[0]=CC(current[0]);
                Cu[1]=CC(current[1]);
                crt[f][0]->Fill(timeH,Cu[0]);
                crt[f][1]->Fill(timeH,Cu[1]);
                
            }
        }
        }
      }
    }
  }
memset(line,0,sizeof(line));

/*memset(time,0,sizeof(time));
for(int i=0;i<4;i++)
memset(header[i],0,sizeof(header[i]));
memset(FEE_ID,0,sizeof(FEE_ID));
memset(temp_flag,0,sizeof(temp_flag));
memset(blank,0,sizeof(blank));
memset(blankS,0,sizeof(blankS));
*/
}
readtemp.close();

char stop[20];
ccopy(stop,time,0,19);
cout<<"Stop!"<<endl;
cout<<"Time:"<<stop<<endl;
//Draw
//pad for time

//TLatex lt[4];
gStyle->SetOptStat(000);
//TCanvas
TCanvas *A[4];
Int_t Qd[4]={1,3,4,2};
for(int q=0;q<4;q++){
char qC[20];
sprintf(qC,"Quadrant%d",Qd[q]);
A[q]=new TCanvas(qC,qC);
}

TLegend *leg[16];
for(int q=0;q<4;q++){
//draw time pad
char latex[80];
sprintf(latex,"Quadrant%d (Start:%s;Stop:%s)",Qd[q],start,stop);
//pad[q]=new TPad("pad","latex",0.1,0.2,0.9,0.9);
//lt[q]=new TLatex(0.1,0.4,latex);
//lt[q].SetTextSize(0.15);
A[q]->Divide(2,2);

for(int f=0;f<4;f++)
{
int i=q*4+f;
A[q]->cd(f+1);

//lt[q].Draw();

gPad->SetGridy(2);
leg[i]=new TLegend(0.25,0.65,0.88,0.85);
leg[i]->SetTextFont(72);
//leg[i]->SetHeader("Termal-sensitive resistance:");
leg[i]->SetHeader(latex);
leg[i]->SetTextSize(0.04);
leg[i]->SetFillColor(kYellow-9);

Rsc[i][0]->Draw();
Rsc[i][0]->SetMarkerStyle(23);
Rsc[i][0]->SetMarkerSize(0.65);
Rsc[i][0]->SetMarkerColor(1);
for(int j=0;j<4;j++){
Rsc[i][j]->Draw("same");
Rsc[i][j]->SetMarkerStyle(23-j);
Rsc[i][j]->SetMarkerSize(0.65);

Rsc[i][j]->SetMarkerColor(j+1);
char legC[40];
if(Layer[i][j]*Bar[i][j]!=0)
sprintf(legC,"T%d_BGO:Layer%d_Bar%d_%dcm",j+1,Layer[i][j],Bar[i][j],Position[i][j]);
else
sprintf(legC,"T%d_FEE:%d",j+1,i+16);

leg[i]->AddEntry(Rsc[i][j],legC,"pl");
}
leg[i]->Draw();
}
}

//Draw Current
    TCanvas *bilibili=new TCanvas("Current","Current");
    bilibili->Divide(4,4);

    TLegend *lege[16];
    char title[80];
    for(int f=0;f<16;f++){
    sprintf(title,"FEE%d (Start:%s;Stop:%s)",f+16,start,stop);
    bilibili->cd(f+1);
    
    gPad->SetGridy();
    lege[f]=new TLegend(0.30,0.75,0.90,0.90);
    lege[f]->SetTextFont(72);
    lege[f]->SetHeader(title);
    lege[f]->SetTextSize(0.04);
    lege[f]->SetFillColor(kYellow-9);
    crt[f][0]->Draw();
    crt[f][0]->SetMarkerStyle(23);
    crt[f][0]->SetMarkerSize(0.65);
    crt[f][0]->SetMarkerColor(4);

    crt[f][1]->Draw("same");
    crt[f][1]->SetMarkerSize(0.65);
    crt[f][1]->SetMarkerStyle(22);
    crt[f][1]->SetMarkerColor(2);
    lege[f]->AddEntry(crt[f][0],"+2.5V","pl");
    lege[f]->AddEntry(crt[f][1],"-2.5V","pl");
    lege[f]->Draw();
}    

    
//Temperature monitors 3D plots
//gStyle->SetOptStat(111);
//TCanvas *temp3D=new TCanvas("3D temp","3D temp");
//Tm->SetMarkerStyle(20);
//Tm->SetMarkerSize(1);
//Tm->SetMarkerColor(2);
//Tm->Draw("COLZ");
//gPad->SetGridx();
//gPad->SetGridy();
//gPad->SetGridz();
return 222;
}
