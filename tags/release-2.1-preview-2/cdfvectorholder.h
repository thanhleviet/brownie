#ifndef CDFVECTORHOLDER_H
#define CDFVECTORHOLDER_H
#include <iostream>
#include <vector>

using namespace std;

class CDFvectorholder
{
	//friend class BROWNIE;
public:
	int x;
	//vector<vector<double> > CDFvectorcontents;
public:
	CDFvectorholder();
	vector<vector<double> > Initialize();
//	~CDFvectorholder();
	//virtual vector<vector<double> > returnCDFvectorcontents();
};

#endif