#include "cdfvectorholder.h"
#include <iostream>
#include <vector>

using namespace std;

CDFvectorholder::CDFvectorholder() {
	x=0;
    //Initialize CDFvectorcontents
	//vector<vector<double> > CDFvectorcontents;
    //these values are based on simulations of 100,000 trees; agreement values for which there were no samples were arbitrarily assigned a value of 1/100000 (since the inverse of the p value is used, this keeps the value from being 1/0)
	
}

vector<vector<double> > CDFvectorholder::Initialize() {
  vector<vector<double> > CDFvectorcontents;
  CDFvectorcontents.push_back( *(new vector<double>(cdfv03, cdfv03+sizeof(cdfv03)/sizeof(cdfv03[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv04, cdfv04+sizeof(cdfv04)/sizeof(cdfv04[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv05, cdfv05+sizeof(cdfv05)/sizeof(cdfv05[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv06, cdfv06+sizeof(cdfv06)/sizeof(cdfv06[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv07, cdfv07+sizeof(cdfv07)/sizeof(cdfv07[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv08, cdfv08+sizeof(cdfv08)/sizeof(cdfv08[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv09, cdfv09+sizeof(cdfv09)/sizeof(cdfv09[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv10, cdfv10+sizeof(cdfv10)/sizeof(cdfv10[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv11, cdfv11+sizeof(cdfv11)/sizeof(cdfv11[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv12, cdfv12+sizeof(cdfv12)/sizeof(cdfv12[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv13, cdfv13+sizeof(cdfv13)/sizeof(cdfv13[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv14, cdfv14+sizeof(cdfv14)/sizeof(cdfv14[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv15, cdfv15+sizeof(cdfv15)/sizeof(cdfv15[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv16, cdfv16+sizeof(cdfv16)/sizeof(cdfv16[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv17, cdfv17+sizeof(cdfv17)/sizeof(cdfv17[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv18, cdfv18+sizeof(cdfv18)/sizeof(cdfv18[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv19, cdfv19+sizeof(cdfv19)/sizeof(cdfv19[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv20, cdfv20+sizeof(cdfv20)/sizeof(cdfv20[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv21, cdfv21+sizeof(cdfv21)/sizeof(cdfv21[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv22, cdfv22+sizeof(cdfv22)/sizeof(cdfv22[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv23, cdfv23+sizeof(cdfv23)/sizeof(cdfv23[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv24, cdfv24+sizeof(cdfv24)/sizeof(cdfv24[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv25, cdfv25+sizeof(cdfv25)/sizeof(cdfv25[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv26, cdfv26+sizeof(cdfv26)/sizeof(cdfv26[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv27, cdfv27+sizeof(cdfv27)/sizeof(cdfv27[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv28, cdfv28+sizeof(cdfv28)/sizeof(cdfv28[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv29, cdfv29+sizeof(cdfv29)/sizeof(cdfv29[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv30, cdfv30+sizeof(cdfv30)/sizeof(cdfv30[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv31, cdfv31+sizeof(cdfv31)/sizeof(cdfv31[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv32, cdfv32+sizeof(cdfv32)/sizeof(cdfv32[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv33, cdfv33+sizeof(cdfv33)/sizeof(cdfv33[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv34, cdfv34+sizeof(cdfv34)/sizeof(cdfv34[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv35, cdfv35+sizeof(cdfv35)/sizeof(cdfv35[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv36, cdfv36+sizeof(cdfv36)/sizeof(cdfv36[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv37, cdfv37+sizeof(cdfv37)/sizeof(cdfv37[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv38, cdfv38+sizeof(cdfv38)/sizeof(cdfv38[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv39, cdfv39+sizeof(cdfv39)/sizeof(cdfv39[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv40, cdfv40+sizeof(cdfv40)/sizeof(cdfv40[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv41, cdfv41+sizeof(cdfv41)/sizeof(cdfv41[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv42, cdfv42+sizeof(cdfv42)/sizeof(cdfv42[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv43, cdfv43+sizeof(cdfv43)/sizeof(cdfv43[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv44, cdfv44+sizeof(cdfv44)/sizeof(cdfv44[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv45, cdfv45+sizeof(cdfv45)/sizeof(cdfv45[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv46, cdfv46+sizeof(cdfv46)/sizeof(cdfv46[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv47, cdfv47+sizeof(cdfv47)/sizeof(cdfv47[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv48, cdfv48+sizeof(cdfv48)/sizeof(cdfv48[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv49, cdfv49+sizeof(cdfv49)/sizeof(cdfv49[0]))) );
  CDFvectorcontents.push_back( *(new vector<double>(cdfv50, cdfv50+sizeof(cdfv50)/sizeof(cdfv50[0]))) );

  return CDFvectorcontents;
}
