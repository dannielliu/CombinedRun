#ifndef Constant_h
#define Constant_h

namespace cosmictest {

const int Nlayer=14;
const int Nside=2;
const int Nbar=22;
const int Nch=24;//22+spare ones of one side of single layer
const int Ndy=3;
const int nGID=Nlayer*Nside*Nch*Ndy;
const int nGbar=Nlayer*Nside*Nbar;
const int nGside=Nlayer*Nside;

const int N_FEE=16;
}
#endif
