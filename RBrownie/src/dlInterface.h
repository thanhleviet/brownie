#ifndef _RBrownie_DLINTERFACE_H
#define _RBrownie_DLINTERFACE_H

#include "dl.h"
#include "browniecmds.h"

class dlInterface
{
	private:
	BROWNIE brownie;
	
	
	public:
	dlInterface();
	~dlInterface();
	
	// Turn on/off Rcpp streaming 
	// (note: turn off for very large files)
	void disableReturnStream();
	
	// lower level functions:
	void pipe( std::string );
	void execute(std::string);
	
	// RETURN ITEM COUNTS:
	int getNumLoadedTrees();
	int getNumTaxa();
	int getNumChars(); // mainly for debugging
	int getNumDiscreteChars();
	int getNumContinuousChars();
	int getNumRetTrees();
	
	// CHARACTERS
	std::vector<std::string> getCharLabels(bool cont = true);
	std::vector<char> getDiscreteChar(int colindex=0);
	std::vector<float> getContChar(int colindex=0);
	
	// TAXA
	int getNumTaxaSets();
	std::vector<std::string> getTaxaSetNames();
	std::vector< std::vector<std::string> > getTaxaSets();
	
	// TREES
	std::string getTree(int,bool=true);  // these two might not return the same trees.
	bool getTree(int, std::ostream &f, bool=true);
	float getTreeWeight(int);
	bool writeTrees(std::string);
	bool hasRetTrees();
	std::string getRetTree(int index=0);
	
	// RETURNED STRINGS (most important):
	std::vector<std::string> getReturnStrings();
	std::vector<std::string> getReturnTrees();
};



#endif
