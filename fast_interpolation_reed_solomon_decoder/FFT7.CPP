//Fedorenko-Trifonov's Direct FFT for GF(2^3) and primitive polynomial 0xb
namespace {
#define MULTIPLY(A,B,C){if (!A||!B) C=0;else {int Log=LogTable[A]+LogTable[B];if (Log>=(FIELDSIZE-1))Log-=(FIELDSIZE-1);C=GF[1+Log];};};
		typedef unsigned GFValue;
	const unsigned EXTENSION = 3;
#define FIELDSIZE (1u<<EXTENSION)
	GFValue GF[8] = { 0,1,2,4,3,6,7,5 };
	int LogTable[8] = { -1,0,1,3,2,6,4,5 };
}
/*Generators of normal bases (LOG notation):
1 0
3 3
*/


void FFT7(GFValue* pSource,GFValue* pDest){
GFValue Products[4];
GFValue LinearizedValues[7];
LinearizedValues[0]=pSource[0];
Products[1]=pSource[2]^pSource[4];
Products[2]=pSource[1]^pSource[2];
Products[3]=pSource[1]^pSource[4];
Products[0]=pSource[1]^Products[1];
MULTIPLY(Products[1],2,Products[1]);
MULTIPLY(Products[2],4,Products[2]);
MULTIPLY(Products[3],6,Products[3]);
LinearizedValues[1]=Products[0];
LinearizedValues[2]=Products[1]^Products[3];
LinearizedValues[3]=Products[2]^Products[3];
Products[1]=pSource[6]^pSource[5];
Products[2]=pSource[3]^pSource[6];
Products[3]=pSource[3]^pSource[5];
Products[0]=pSource[3]^Products[1];
MULTIPLY(Products[1],2,Products[1]);
MULTIPLY(Products[2],4,Products[2]);
MULTIPLY(Products[3],6,Products[3]);
LinearizedValues[4]=Products[0];
LinearizedValues[5]=Products[1]^Products[3];
LinearizedValues[6]=Products[2]^Products[3];
GFValue Stage2[6];
Stage2[0]=LinearizedValues[3]^LinearizedValues[6];
Stage2[2]=LinearizedValues[2]^LinearizedValues[6];
Stage2[3]=LinearizedValues[1]^LinearizedValues[4];
Stage2[4]=LinearizedValues[0]^LinearizedValues[4];
Stage2[5]=LinearizedValues[0]^LinearizedValues[1];
pDest[0]=LinearizedValues[0]^Stage2[3];
pDest[1]=Stage2[0]^Stage2[4];
pDest[5]=Stage2[2]^Stage2[5];
Stage2[1]=LinearizedValues[5]^pDest[1];
pDest[2]=LinearizedValues[2]^Stage2[1];
pDest[3]=Stage2[1]^Stage2[3];
pDest[4]=pDest[2]^Stage2[0];
pDest[6]=pDest[3]^Stage2[2];
};
//6 multiplications, 25 additions
//Construction time 0.390625 s
