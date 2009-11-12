/*
 *  treepiper.h
 *  
 *
 *  Created by Brian O'Meara on Fri Jul 07 2006.
 *  Copyright (c) 2006 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef TREEPIPERH
#define TREEPIPERH
#include "treewriter.h"
#include "taxablock.h"
#include "treesblock.h"

class TreePiper : public NewickTreeWriter
{
	friend class BROWNIE;
public:
	//BROWNIE brownie;
    TreePiper (Tree *tree, TaxaBlock *taxa, TreesBlock *trees, vector<int> convertsamplestospecies, int chosentree): NewickTreeWriter (tree) {t=tree; writeEdgeLengths = false;};
	virtual ~TreePiper () {};
	vector<int> convertsamplestospecies;
	int chosentree;
	TaxaBlock *taxa;
	TreesBlock *trees;
	virtual nxsstring PipeLeafGTP ();
	virtual nxsstring PipeLeafSpeciesTree ();
	virtual nxsstring PipeGTP ();
	virtual nxsstring PipeSpeciesTree ();
	virtual nxsstring PipeEndOfTreeGTP ();
	virtual nxsstring PipeEndOfTreeSpeciesTree ();	
	virtual nxsstring PipeLeftParenthesis();
	virtual nxsstring PipeRightParenthesis();
	virtual nxsstring PipeSiblingSymbol();
	virtual nxsstring PipeInternal();
};

#endif