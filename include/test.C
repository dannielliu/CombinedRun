#include "timeconvertor.h"
int test(){
char *time="2014-06-22-07-05-00";
Long64_t t=Dat2Ast(time);
cout<<"Absolute"<<t<<endl;
cout<<Ast2Dat(t)<<endl;
return;
}
