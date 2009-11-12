/*
 *  treepiper.cpp
 *  
 *
 *  Created by Brian O'Meara on Fri Jul 07 2006.
 *  Copyright (c) 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include <strstream>
#include <fstream.h>
#include <iomanip.h>
#include <unistd.h>
#include <stdio.h>
#include <set>
#include <math.h>

#include "nexusdefs.h"
#include "xnexus.h"
#include "nexustoken.h"
#include "nexus.h"
#include "taxablock.h"
#include "assumptionsblock.h"
#include "treesblock.h"
#include "discretedatum.h"
#include "discretematrix.h"
#include "charactersblock.h"
#include "gport.h"
#include "profile.h"
#include "ntree.h"
#include "nodeiterator.h"
#include "treeorder.h"
#include "treedrawer.h"
#include "TreeLib.h"
#include "gtree.h"
#include "treereader.h"
#include "treewriter.h"
#include <time.h>
#include <map>
#include "brownie.h"
#include "optimizationfn.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_exp.h>
#include "treepiper.h"



nxsstring TreePiper::PipeLeafGTP ()
{
	int ntax=taxa->GetNumTaxonLabels();
	cout<<"Num taxon labels in PipeLeafGTP="<<ntax<<endl;
	nxsstring output="";
	nxsstring TaxonLabel="";
	TaxonLabel+=cur->GetLabel();
	int TaxonNumber=taxa->FindTaxon(TaxonLabel);
	nxsstring NewLabel="taxon";
	NewLabel+=convertsamplestospecies[TaxonNumber];
	output=NEXUSString (NewLabel);
	return output;
}

nxsstring TreePiper::PipeLeafSpeciesTree ()
{
	nxsstring output="";
	output+=(cur->GetLabel());
	return output;
}

nxsstring TreePiper::PipeGTP ()
{
	nxsstring TreeDescription="";
	cur = t->GetRoot();
	
    while (cur)
    {
        if (cur->GetChild())
        {
            TreeDescription+=PipeLeftParenthesis ();
            stk.push (cur);
            cur = cur->GetChild();
        }
        else
        {
            TreeDescription+=PipeLeafGTP ();
            while (!stk.empty() && (cur->GetSibling() == NULL))
            {
                TreeDescription+=PipeRightParenthesis ();
                cur = stk.top();
                PipeInternal ();
                stk.pop();
            }
            if (stk.empty())
                cur = NULL;
            else
            {
                TreeDescription+=PipeSiblingSymbol ();
                cur = cur->GetSibling();
            }
        }
    }
    TreeDescription+=PipeEndOfTreeGTP ();
	return TreeDescription;
}


nxsstring TreePiper::PipeSpeciesTree ()
{
	nxsstring TreeDescription="";
	cur = t->GetRoot();
	
    while (cur)
    {
        if (cur->GetChild())
        {
            TreeDescription+=PipeLeftParenthesis ();
            stk.push (cur);
            cur = cur->GetChild();
        }
        else
        {
            TreeDescription+=PipeLeafSpeciesTree ();
            while (!stk.empty() && (cur->GetSibling() == NULL))
            {
                TreeDescription+=PipeRightParenthesis ();
                cur = stk.top();
                TreeDescription+=PipeInternal ();
                stk.pop();
            }
            if (stk.empty())
                cur = NULL;
            else
            {
                TreeDescription+=PipeSiblingSymbol ();
                cur = cur->GetSibling();
            }
        }
    }
    TreeDescription+=PipeEndOfTreeSpeciesTree ();
	return TreeDescription;
}



nxsstring TreePiper::PipeEndOfTreeGTP ()
{
	nxsstring output="";
	output+=':';
	output+=trees->GetTreeWeight(chosentree);
	output+=';';
	return output;
}

nxsstring TreePiper::PipeEndOfTreeSpeciesTree ()
{
	nxsstring output="";
	output+=';';
	return output;
}


//------------------------------------------------------------------------------
nxsstring TreePiper::PipeLeftParenthesis ()
{
	nxsstring output="";
	output+='(';
	return output;
}

//------------------------------------------------------------------------------
nxsstring TreePiper::PipeRightParenthesis ()
{
	nxsstring output="";
	output+=')';
	return output;
}

//------------------------------------------------------------------------------
nxsstring TreePiper::PipeSiblingSymbol ()
{
	nxsstring output="";
	output+=',';
	return output;
}


//------------------------------------------------------------------------------
nxsstring TreePiper::PipeInternal ()
{
	nxsstring output="";
	if (cur->GetLabel() != "")
		output=NEXUSString (cur->GetLabel());
	if (t->GetHasEdgeLengths () && writeEdgeLengths)
	{
		output+=":";
	output+=cur->GetEdgeLength();
	}
	return output;
}
