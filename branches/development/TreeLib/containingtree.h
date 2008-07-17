#ifndef CONTAININGTREEH
#define CONTAININGTREEH

#include "stree.h"

#include <string>
#include <stack>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "nexusdefs.h"
#include "xnexus.h"
#include "nexustoken.h"
#include "nexus.h"
extern gsl_rng * r;

/*
 *  containingtree.h
 *
 *
 *  Created by Brian O'Meara on Tue Jul 04 2006.
 *  Copyright (c) 2006 __MyCompanyName__. All rights reserved.
 *
 */

/**
* @class ContainingTree
 * ContainingTree is a tree class to support trees that contain input trees (i.e., species trees for gene trees.
                                                                             * It extends STree to add a stack for keeping track of the growing supertree.
                                                                             *
                                                                             */
class ContainingTree : public STree
{
public:
    ContainingTree () { stk.empty(); };
	
	//destructor
	virtual ~ContainingTree();
    
    //Create a random n-taxon tree
    virtual void RandomTree (int ntax);
    
	//Speciate at a selected node
    void SpeciateAtLeaf (NodePtr Node);
    
    //Do an NNI move
    virtual void NNI ();
	
	//Take a given node and do an NNI
	void NonRandomNNIAtNode (int chosennodenumber, int resolution);
    
    //Do an SPR move
    virtual void SPR ();

    //Do next possible SPR move; if this is the last possible move, returns false
    virtual bool NextSPR();

    //Generates break vectors for the possible moves
    virtual void ResetBreakVector();

    //Generates break vectors for the possible moves
    virtual void ResetAttachVector();

    virtual void AttachTraverse(NodePtr p);

    NodePtr NodeToBreak;
    NodePtr NodeToAttach;
        //Vector containing a list of nodes that can be broken on
    vector<int> BreakVector;
    //Vector containing a list of nodes that can be reattached on (reinitialized for each break point)
    vector<int> AttachVector;

    vector<int> GetBreakVector() {return BreakVector; };
    vector<int> GetAttachVector() {return AttachVector; };

    virtual void SetBreakVector(vector<int> inputvector) {BreakVector=inputvector; };
    virtual void SetAttachVector(vector<int> inputvector) {AttachVector=inputvector; };
    //Vector containining a list of nodes to reroot on (TBR only)
    vector<int> ReRootVector;

    virtual void GetNodeToBreak();

    virtual void GetNodeToAttach();

    void FindAttachTraverse (NodePtr p);

    int NodeToBreakInt;
    
    int NodeToAttachInt;

    int nodecount;

    NodePtr GetLeafWithNumber (int i);

    vector<nxsstring> GetLeafLabelVector();

    vector<int> GetLeafNumberVector();

    std::map<Node *, int, std::less<Node *> > depth;

    vector<int> GetLCADepthVector(nxsstring a, nxsstring b, nxsstring c);

    int LCADepthQuery(nxsstring a, nxsstring b);
    
    //Do a TBR move
    virtual void TBR ();
    
    //Split a leaf into a cherry (so, split one taxon into two). The first number is the number of the taxon that was split, the second number is the number of the new taxon
    virtual vector<int> SplitLeaf(int LeafToSplit);
    
    //Collapse a cherry into a leaf (combine two neighboring taxa into one). Note that I used a relaxed cherry definition, but only two taxa are collapsed at this step. The first number of the vector is the taxon number of the taxon the pair has been joined into, the second number is the deleted taxon.
    virtual vector<int> CollapseCherry(int chosencherrynum);

    int LeavesCreated;

    vector <NodePtr> CherryNodes; //Contains pointers to all the cherries (needed for cherry contraction)

    virtual int GetNumCherries() { return CherryNodes.size(); };

	//Convert names on tree to sequential numbers
    virtual void ConvertTaxonNamesToOrderedTaxonNumbers();

    //Convert taxon names to random taxon numbers
    virtual void ConvertTaxonNamesToRandomTaxonNumbers();
    
    //Delete the edge below a node. Also eliminates any nodes with just one descendant as a result.
    virtual void DeleteSubtendingEdge(NodePtr SelectedNode);

    //Reroot a tree on the chosen node and return a pointer to the new node
    virtual NodePtr ReRootTree(NodePtr SelectedNode);

    //Find and then suppress internal nodes with one descendant (generated during rerooting or pruning)
    virtual void SuppressInternalNodesWithOneDescendant();

    virtual void SuppressInternalNode(NodePtr SelectedNode);

    virtual void TestRerooting();
	
	virtual void SetLeafNumbers();

        virtual void FixLeafNamesAndNumbers(int deletedleafnumber);

    //Add a node to the branch subtending the given node. This node will have just one descendant.
    virtual NodePtr AddNodeToSubtendingEdge (NodePtr Node);

    virtual int TraverseToGetNodeCount (NodePtr p);

    virtual void UpdateCherries();

    virtual void FindAndSetRoot();

    virtual void ReportTreeHealth();

    virtual vector<ContainingTree> SplitOnTaxon (vector<nxsstring> taxatoexclude);

    virtual vector<int> GetPotentialNewRoots();

    virtual NodePtr SelectNodeToReRootOn(int chosennodenum);

    virtual void ClearInternalLabels();
	
	virtual void InitializeMissingBranchLengths();

	virtual void RandomlyModifySingleBranchLength(double markedmultiplier, double brlensigma);
	
	virtual void NodeSlideBranchLength(double markedmultiplier);
	
	virtual void ModifyTotalBranchLength(double brlensigma);
};



#endif
