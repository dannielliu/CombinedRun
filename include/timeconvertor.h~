#ifndef timeconvertor_h
#define timeconvertor_h
#include <time.h>

char *Ast2Dat(Long64_t Astime){
//  struct tm DatT0={0,0,0,1,0,113,1,0,-1};//{sec,min,hour,mday(1-31),mon(0-11),year(real-1900),wday(0-6),yday(0-365),isdst(+,0,-)}
  //conminder time zone
  struct tm DatT0={0,0,8,1,0,113,1,0,-1};//{sec,min,hour,mday(1-31),mon(0-11),year(real-1900),wday(0-6),yday(0-365),isdst(+,0,-)}
  time_t Astime0=mktime(&DatT0);
  time_t Astime1=Astime0+Astime;
  struct tm *DatT=gmtime(&Astime1);
  Int_t year=DatT->tm_year+1900;
  Int_t month=DatT->tm_mon+1;
  Int_t day=DatT->tm_mday;
  Int_t hour=DatT->tm_hour;
  Int_t min=DatT->tm_min;
  Int_t sec=DatT->tm_sec;
  char* date=(char*)malloc(30);
//  char* date=new char[30];
  sprintf(date,"%d/%02d/%02d %02d:%02d:%02d",year,month,day,hour,min,sec);
  return date;
}
//define char 2 int
Int_t c2i(char ch){
Int_t In=(int)ch-48;
return In;
}

Long64_t Dat2Ast(char *date){
  struct tm DatT0={0,0,0,1,0,113,1,0,-1};//{sec,min,hour,mday(1-31),mon(0-11),year(real-1900),wday(0-6),yday(0-365),isdst(+,0,-)}
//struct tm DatT0={0,0,14,31,11,112,0,365,-1};//{sec,min,hour,mday(1-31),mon(0-11),year(real-1900),wday(0-6),yday(0-365),isdst(+,0,-)}
  time_t Astime0=mktime(&DatT0);
  struct tm DatT;
//time YYYY-MM-DD HH:MM:SS
  Int_t year=c2i(date[0])*1000+c2i(date[1])*100+c2i(date[2])*10+c2i(date[3]);
  DatT.tm_year=year-1900;
  Int_t month=c2i(date[5])*10+c2i(date[6]);
  DatT.tm_mon=month-1;
  Int_t day=c2i(date[8])*10+c2i(date[9]);
  DatT.tm_mday=day;
  Int_t hour=c2i(date[11])*10+c2i(date[12]);
  DatT.tm_hour=hour;
  Int_t min=c2i(date[14])*10+c2i(date[15]);
  DatT.tm_min=min;
  Int_t sec=c2i(date[17])*10+c2i(date[18]);
  DatT.tm_sec=sec;
  time_t Astime1=mktime(&DatT);
  Long64_t Astime=(Long64_t)difftime(Astime1,Astime0);

  return Astime;
}
//char *asctime(struct tm *p)
//char *ctime(long time);
//time_t mktime(struct tm *p)
//stuct mt *gmtime(time_t *t)

#endif
