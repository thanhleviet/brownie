/*
 *  containingtree.cpp
 *
 *
 *  Created by Brian O'Meara on Tue Jul 04 2006.
 *  Copyright (c) 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "containingtree.h"
#include "nexusdefs.h"
#include "xnexus.h"
#include "nexustoken.h"
#include "nexus.h"
#include "nodeiterator.h"
#include "treedrawer.h"
#include "profile.h"
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
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_exp.h>
#include <math.h>
#include <string>
#include <stack>
#include "brownie.h"


ContainingTree::~ContainingTree() {
	//deletetraverse(Root);
	//delete [] Nodes;
	//Nodes = NULL;
}

void ContainingTree::ResetBreakVector()
{
    FindAndSetRoot();
    if (GetNumLeaves()>1) {
        vector<int> TempBreakVector;
        BreakVector.clear();
        NodeIterator <Node> n (GetRoot());
        NodePtr currentnode = n.begin();
        nodecount=0;
        while (currentnode)
        {
            if (currentnode!=GetRoot()) {
                TempBreakVector.push_back(nodecount);
            }
            nodecount++;
            currentnode = n.next();
        }
        gsl_permutation * p = gsl_permutation_alloc (TempBreakVector.size()); //all this is so that the nodes are examined in random order
        gsl_permutation_init (p);
        gsl_ran_shuffle (r, p->data, TempBreakVector.size(), sizeof(size_t));
        for (int i=0; i<TempBreakVector.size(); i++) {
            BreakVector.push_back(TempBreakVector[gsl_permutation_get (p,i)]); //so, assigning nodes in random order to the BreakVector
			//cout<<"TempBreakVector[gsl_permutation_get (p,i)] = "<<TempBreakVector[gsl_permutation_get (p,i)]<<endl;
        }
        //cout<<"BreakV =";
        //for (int j=0; j<BreakVector.size(); j++) {
        //    cout<<" "<<BreakVector[j];
        //}
        //cout<<"BreakVector.back() = "<<BreakVector.back()<<endl;
        //cout<<"BreakVector.back() = "<<BreakVector.back()<<endl;
        //cout<<endl;
        gsl_permutation_free(p);
        NodeToBreakInt=BreakVector.back();
        GetNodeToBreak();
        ResetAttachVector();
    }
}

void ContainingTree::ResetAttachVector()
{
    FindAndSetRoot();
    if(GetNumLeaves()>1) {
        vector<int> TempAttachVector;
        AttachVector.clear();
        int OriginalNodeToBreakInt=NodeToBreakInt;
        NodeToBreakInt=BreakVector.back();
        //cout<<"ResetAttachVector BreakVector.back() = "<<BreakVector.back()<<endl;
        //cout<<"NodeToBreakInt = "<<NodeToBreakInt<<endl;
        GetNodeToBreak();
        nodecount=0;
        AttachTraverse(Root); //Fills AttachVector with nodes in order (without needing to actually delete the subtree)
        TempAttachVector.swap(AttachVector); //Now AttachVector is empty again and TempAttachVector has nodes in order
        gsl_permutation * p = gsl_permutation_alloc (TempAttachVector.size()); //all this is so that the nodes are examined in random order
        gsl_permutation_init (p);
        gsl_ran_shuffle (r, p->data, TempAttachVector.size(), sizeof(size_t));
        //cout<<"AttachV = ";
        for (int i=0; i<TempAttachVector.size(); i++) {
            AttachVector.push_back(TempAttachVector[gsl_permutation_get (p,i)]); //so, assigning nodes in random order to the AttachVector
			//cout<<TempAttachVector[gsl_permutation_get (p,i)]<<" ";
        }
        //cout<<endl;
        gsl_permutation_free(p);
        NodeToBreakInt=OriginalNodeToBreakInt; //reset to initial values
        GetNodeToBreak();
    }
}

void ContainingTree::AttachTraverse (NodePtr p)
{
    if (p)
    {
        // cout<<"AttachTraverse on node "<<p<<endl;
        if (p!=NodeToBreak) { //so we don't grab the node to break on or its descendants
            AttachTraverse (p->GetChild());
            if (p!=Root && (p!=(NodeToBreak->GetAnc()))) { //Since we're basically dealing with unrooted species trees, don't want to bother adding things below the root AND we also don't want to bother breaking below a node and then reattaching in the same place
                AttachVector.push_back(nodecount);
            }
            nodecount++;
        }
        AttachTraverse (p->GetSibling());
    }
}

void ContainingTree::FindAttachTraverse (NodePtr p)
{
    if (p)
    {
        //cout<<"FindAttachTraverse on node "<<p<<endl;
        if (p!=NodeToBreak) { //so we don't grab the node to break on or its descendants
            FindAttachTraverse (p->GetChild());
            if (nodecount==NodeToAttachInt) {
                NodeToAttach=p;
                //cout<<"Found NodeToAttach = "<<NodeToAttach<<" "<<nodecount<<endl;
            }
            nodecount++;
        }
        FindAttachTraverse (p->GetSibling());
    }
}

bool ContainingTree::NextSPR()
{
    FindAndSetRoot();
    Update();
    //cout<<"BreakVector now: ";
    //for (int j=0; j<BreakVector.size(); j++) {
    //    cout<<" "<<BreakVector[j];
    //}
    //cout<<endl;
	bool ThereAreMoreMoves=true;
	NodeToBreakInt=BreakVector.back();
	//cout<<"BreakVector.back() = "<<BreakVector.back()<<endl;
	//cout<<"NodeToBreakInt = "<<NodeToBreakInt<<endl;
	NodeToAttachInt=AttachVector.back();
	AttachVector.pop_back();
	if (AttachVector.empty()) { //We've run out of attachment points with that breakpoint
		BreakVector.pop_back(); //So delete that breakpoint
		//cout<<"after popback, BreakVector now: ";
		//for (int j=0; j<BreakVector.size(); j++) {
		//    cout<<" "<<BreakVector[j];
		//}
		//cout<<endl<<endl;
		if (BreakVector.empty()) { //We've also run out of break points
			ThereAreMoreMoves=false;
			//  cout<<"\n\n\n*************************Done with that initial tree**********************\n\n\n";
		}
		else { //There are still break points, but now need a new AttachVector
			ResetAttachVector();
		}
	}
	//Now we can actually do this SPR
	//cout<<"NodeToBreakInt = "<<NodeToBreakInt<<endl;
	//cout<<"NodeToAttachInt = "<<NodeToAttachInt<<endl;
	GetNodeToBreak();
	GetNodeToAttach();
	//cout<<"NodeToBreak = "<<NodeToBreak<<endl;
	//cout<<"NodeToAttach = "<<NodeToAttach<<endl;
	assert(NodeToAttach!=NULL);
	assert(NodeToBreak!=NULL);
	assert(NodeToAttach!=NodeToBreak->GetAnc());
	DeleteSubtendingEdge(NodeToBreak);
	//cout<<endl<<"Just deleted subtending edge"<<endl;
	//ReportTreeHealth();
	FindAndSetRoot();
	NodePtr NewPosition=AddNodeToSubtendingEdge(NodeToAttach);
	(NewPosition->GetChild())->SetSibling(NodeToBreak); //We know NewPosition has one child with no siblings as a result of how we created New Position
	NodeToBreak->SetAnc(NewPosition);
	NewPosition->IncrementDegree();
	//cout<<endl<<endl;
	//ReportTreeHealth();
	//  cout<<endl<<endl;
	FindAndSetRoot();
	Update();
	return ThereAreMoreMoves;
}

void ContainingTree::GetNodeToBreak() //We need this rigamarole because after doing SPR, the tree has changed, so if we give nextSPR the original tree again, the node ptrs may have changed, but not the order.
{
    //ReportTreeHealth();
    int NodesTouched=0;
    FindAndSetRoot();
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    bool selectedanode=false;
    while (currentnode)
    {
        // cout<<"GetNodeToBreak\t"<<NodesTouched<<"\t"<<NodeToBreakInt<<endl;
        if (NodeToBreakInt==NodesTouched) { // We've found the node to split on
            NodeToBreak=currentnode;
            selectedanode=true;
            // cout<<"Selected node to break is "<<NodeToBreak<<endl;
            break;
        }
        NodesTouched++;
        currentnode = n.next();
    }
    if (selectedanode==false) {
        cout<<"Error -- didn't select a node for a break"<<endl<<endl;
        cout<<"NodeToBreakInt = "<<NodeToBreakInt<<endl;
        cout<<"NodeToAttachInt = "<<NodeToAttachInt<<endl;
        cout<<"NodeToBreak = "<<NodeToBreak<<endl;
        cout<<"NodeToAttach = "<<NodeToAttach<<endl;
        ReportTreeHealth();
        cout<<endl<<endl;
    }
    assert(selectedanode);
}

void ContainingTree::GetNodeToAttach()
{
    nodecount=0;
    FindAndSetRoot();
    FindAttachTraverse(Root);
}

void ContainingTree::ReportTreeHealth()
{
    cout<<"NumLeaves = "<<GetNumLeaves()<<endl;
    cout<<"NumInternals = "<<GetNumInternals()<<endl;
    PreorderIterator <Node> n (Root);
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if (currentnode->IsLeaf())
        {
            cout<<currentnode<<"=Leaf "<<currentnode->GetLabel()<<"\tAnc="<<currentnode->GetAnc();
            if ((currentnode->GetAnc())!=NULL) {
                cout<<" "<<(currentnode->GetAnc())->GetLabel();
            }
            if ((currentnode->GetSibling())!=NULL) {
                cout<<"\tSib="<<currentnode->GetSibling()<<" "<<(currentnode->GetSibling())->GetLabel();
            }
            else {
                cout<<"\tSib=NULL";
            }
            if ((currentnode->GetChild())!=NULL) {
                cout<<"\tChild="<<currentnode->GetChild()<<" "<<(currentnode->GetChild())->GetLabel();
            }
            else {
                cout<<"\tChild=NULL";
            }
            cout<<"\tDepth="<<currentnode->GetDepth()<<"\tHeight="<<currentnode->GetHeight();
        }
        else if (currentnode==Root)
        {
            cout<<currentnode<<"=Root "<<currentnode->GetLabel();
            if ((currentnode->GetAnc())!=NULL) {
                cout<<"\tAnc="<<currentnode->GetAnc();
            }
            else {
                cout<<"\tAnc=NULL";
            }
            if ((currentnode->GetSibling())!=NULL) {
                cout<<"\tSib="<<currentnode->GetSibling()<<" "<<(currentnode->GetSibling())->GetLabel();
            }
            else {
                cout<<"\tSib=NULL";
            }
            if ((currentnode->GetChild())!=NULL) {
                cout<<"\tChild="<<currentnode->GetChild()<<" "<<(currentnode->GetChild())->GetLabel();
            }
            else {
                cout<<"\tChild=NULL";
            }
            cout<<"\tDepth="<<currentnode->GetDepth()<<"\tHeight="<<currentnode->GetHeight();
        }
        else
        {
            cout<<currentnode<<"=Internal "<<currentnode->GetLabel()<<"\tAnc="<<currentnode->GetAnc();
            if ((currentnode->GetAnc())!=NULL) {
                cout<<" "<<(currentnode->GetAnc())->GetLabel();
            }
            if ((currentnode->GetSibling())!=NULL) {
                cout<<"\tSib="<<currentnode->GetSibling()<<" "<<(currentnode->GetSibling())->GetLabel();
            }
            else {
                cout<<"\tSib=NULL";
            }
            if ((currentnode->GetChild())!=NULL) {
                cout<<"\tChild="<<currentnode->GetChild()<<" "<<(currentnode->GetChild())->GetLabel();
            }
            else {
                cout<<"\tChild=NULL";
            }
            cout<<"\tDepth="<<currentnode->GetDepth()<<"\tHeight="<<currentnode->GetHeight();
        }
		if (GetHasEdgeLengths()) {
			cout<<" Edge length = "<<currentnode->GetEdgeLength();
			if (((currentnode->GetEdgeLength())!=(currentnode->GetEdgeLength())) || gsl_isnan(currentnode->GetEdgeLength()) || gsl_isinf(currentnode->GetEdgeLength())) {
				cout<<" <-is nan";
			}
		}
		cout<<endl;
        currentnode = n.next();
    }
	
}

void ContainingTree::RandomTree (int ntax) {
    LeavesCreated=0;
    CurNode   = NewNode();
    Root 	  = CurNode;
    //cout<<"CurNode = "<<CurNode<<endl<<endl;
    MakeChild();
    LeavesCreated++;
	
    nxsstring buf;
    buf="t";
    buf+=LeavesCreated;
    AddLeaf(buf);
    if (ntax>1) {
        MakeSibling();
        LeavesCreated++;
        buf="t";
        buf+=LeavesCreated;
        AddLeaf(buf);
    }
    for (int nodenumber=0;nodenumber<(ntax-2);nodenumber++) { //A rooted tree with ntax has ntax-1 internal nodes; we've already created 1
		//Here, we'll randomly select a leaf and cause it to speciate (change it to an internal node).
		//cout<<endl<<endl<<endl<<"Number of leaves = "<<GetNumLeaves()<<endl;
        int chosenleafnumber=int(ceil(gsl_ran_flat (r,0,GetNumLeaves())));
        //cout<<chosenleafnumber<<"=chosenleafnumber"<<endl;
        int currentleafnumber=0;
        NodeIterator <Node> n (GetRoot());
        NodePtr NodeToSpeciate = n.begin();
        //cout<<endl<<"Starting NodeToSpeciate="<<NodeToSpeciate<<endl;
        bool successfulselection=false;
        while (NodeToSpeciate)
        {
            if (NodeToSpeciate->IsLeaf())
            {
                currentleafnumber++;
                //cout<<endl<<currentleafnumber<<": Now examining "<<NodeToSpeciate;
                if (currentleafnumber==chosenleafnumber) {
                    //cout<<" CHOSEN NODE";
                    successfulselection=true;
                    break;
                }
            }
            NodeToSpeciate = n.next();
        }
        if (successfulselection) {
            SpeciateAtLeaf(NodeToSpeciate);
        }
        else {
            //cout<<endl<<"Error: failed to select a node, though NodeToSpeciate="<<NodeToSpeciate<<endl;
        }
    }
}

//Involves deleting the leaf, having the ancestor create a new node,
//and having that new node have two children that each becomes a leaf.
void ContainingTree::SpeciateAtLeaf (NodePtr Node)
{
    assert(Node->IsLeaf());
    NodePtr Ancestor = Node->GetAnc();
    int startingleafnumber= Node->GetLeafNumber();
    NodePtr q;
    if (Node->IsTheChild())
    {
        Ancestor->SetChild (Node->GetSibling());
        q = Node->GetSibling ();
    }
    else
    {
        q = Node->LeftSiblingOf ();
        q->SetSibling (Node->GetSibling ());
    }
    Node->SetSibling (NULL);
    Node->SetAnc (NULL); //So we've deleted the node that speciates
    Ancestor->SetDegree (Ancestor->GetDegree() - 1);
    Leaves--;
    CurNode=Ancestor->GetChild();
    MakeSibling(); //We replace the leaf that we've just deleted with an internal node (that's sister to its ancestor's daughter).
	//Now we create names for the new taxa
    nxsstring NewLeaf1="t";
    LeavesCreated++;
    NewLeaf1+=LeavesCreated;
    nxsstring NewLeaf2="t";
    LeavesCreated++;
    NewLeaf2+=LeavesCreated;
    AddCherry(NewLeaf1,NewLeaf2);
}

void ContainingTree::ConvertTaxonNamesToOrderedTaxonNumbers()
{
    int LeavesTouched=0;
    nxsstring NewLeafLabel;
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if (currentnode->IsLeaf())
        {
            LeavesTouched++;
            NewLeafLabel="taxon";
            NewLeafLabel+=LeavesTouched;
            currentnode->SetLabel(NewLeafLabel);
            currentnode->SetLeafNumber(LeavesTouched);
        }
        currentnode = n.next();
    }
	
}

void ContainingTree::ConvertTaxonNamesToRandomTaxonNumbers()
{
    const size_t N = GetNumLeaves();
    gsl_permutation * p = gsl_permutation_alloc (N);
    gsl_permutation_init (p);
    gsl_ran_shuffle (r, p->data, N, sizeof(size_t));
    int LeavesTouched=0;
    nxsstring NewLeafLabel;
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if (currentnode->IsLeaf())
        {
            LeavesTouched++;
            NewLeafLabel="taxon";
            NewLeafLabel+=int(1+gsl_permutation_get(p,LeavesTouched-1));
            currentnode->SetLabel(NewLeafLabel);
            currentnode->SetLeafNumber(LeavesTouched);
        }
        currentnode = n.next();
    }
    gsl_permutation_free(p);
}


void ContainingTree::NNI () {
}

//This function does an NNI over the edge below the given node. This node should not be the root node, a child of the root node, or a leaf. There are two possible NNIs for an edge that has each end node of degree 2, resolution=1 or resolution=2 selects which to return
void ContainingTree::NonRandomNNIAtNode (int chosennodenumber, int resolution) { 
	
	FindAndSetRoot();
	NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    NodePtr SelectedEdgeTipwardNode;
	int NodesTouched=0;
    bool selectedanode=false;
    while (currentnode)
    {
        NodesTouched++;
        //  cout<<chosennodenumber<<"\t"<<NodesTouched<<"\t"<<GetNumNodes()<<endl;
        if (chosennodenumber==NodesTouched) { // We've found the node to split on
            if (currentnode!=GetRoot()) {
                SelectedEdgeTipwardNode=currentnode;
                selectedanode=true;
                break;
            }
            else {
                chosennodenumber++; //Don't want to split on the root node
				//  cout<<"Chosen node was root; choosing next instead"<<endl;
            }
        }
        currentnode = n.next();
    }
    if (selectedanode==false) {
        cout<<"Error -- didn't select a node for a break"<<endl<<endl;
    }
	
	if (resolution==1) {
		NodeToBreak=SelectedEdgeTipwardNode->GetChild();
	}
	else {
		NodeToBreak=(SelectedEdgeTipwardNode->GetChild())->GetSibling();
	}
	NodePtr NodeToAdd;
	if((SelectedEdgeTipwardNode->GetAnc()->GetChild())==SelectedEdgeTipwardNode) {
		NodeToAdd=(SelectedEdgeTipwardNode->GetAnc()->GetChild())->GetSibling();
	}
	else {
		NodeToAdd=(SelectedEdgeTipwardNode->GetAnc()->GetChild());
	}
	if( (SelectedEdgeTipwardNode->GetDegree() > 2) || ((SelectedEdgeTipwardNode->GetAnc())->GetDegree() > 2) ) {
		cout<<"Error: this edge has a multifurcation at at least one end";
		DeleteSubtendingEdge(NodeToBreak); // just to break the tree
	}
	else if (NodeToBreak == NULL || NodeToAdd == NULL) {
		cout<<"Error: NodeToBreak "<<NodeToBreak<<" or NodeToAdd "<<NodeToAdd<<" is NULL"<<endl;
		ReportTreeHealth();
	}
	else {
		DeleteSubtendingEdge(NodeToBreak);
		NodePtr NewPosition=AddNodeToSubtendingEdge(NodeToAdd);
		(NewPosition->GetChild())->SetSibling(NodeToBreak); //We know NewPosition has one child with no siblings as a result of how we created New Position
		NodeToBreak->SetAnc(NewPosition);
		NewPosition->IncrementDegree();
		FindAndSetRoot();
		Update();	
	}
}

void ContainingTree::SPR () {
    //First, select a node to break under
    //cout<<endl<<"Now starting SPR"<<endl;
    int chosennodenumber=int(ceil(gsl_ran_flat (r,0,GetNumNodes()-1))); //Since we don't want the root node
    int NodesTouched=0;
    FindAndSetRoot();
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    NodeToBreak;
    bool selectedanode=false;
    while (currentnode)
    {
        NodesTouched++;
        //  cout<<chosennodenumber<<"\t"<<NodesTouched<<"\t"<<GetNumNodes()<<endl;
        if (chosennodenumber==NodesTouched) { // We've found the node to split on
            if (currentnode!=GetRoot()) {
                NodeToBreak=currentnode;
                selectedanode=true;
                break;
            }
            else {
                chosennodenumber++; //Don't want to split on the root node
				//  cout<<"Chosen node was root; choosing next instead"<<endl;
            }
        }
        currentnode = n.next();
    }
    if (selectedanode==false) {
        cout<<"Error -- didn't select a node for a break"<<endl<<endl;
    }
	
    //cout<<"NodeToBreak = "<<NodeToBreak<<" "<<NodeToBreak->GetLabel()<<endl;
    //Tree SubtreeCopy;
    //SubtreeCopy.SetRoot(CopyOfSubtree(NodeToBreak));
    //SubtreeCopy.Update();
    //cout<<"Subtree: "<<endl;
    //SubtreeCopy.Draw(cout);
    //cout<<endl<<endl;
	
    DeleteSubtendingEdge(NodeToBreak);
	
    //Draw(cout);
    //int numberremainingnodes=GetNumNodes()-TraverseToGetNodeCount(NodeToBreak);
    int numberremainingnodes=GetNumNodes();
    int chosennodenumberforinsertion=int(ceil(gsl_ran_flat (r,0,numberremainingnodes)));
    NodesTouched=0;
    FindAndSetRoot();
    NodeIterator <Node> m (GetRoot());
    NodePtr currentnodeforinsertion = m.begin();
    NodePtr NodeToAdd;
    //selectedanode=false;
    while (currentnodeforinsertion)
    {
        NodesTouched++;
        //cout<<chosennodenumberforinsertion<<"\t"<<NodesTouched<<"\t"<<GetNumNodes()<<"\t"<<numberremainingnodes<<"\t"<<TraverseToGetNodeCount(NodeToBreak)<<endl;
        //cout<<"Now examining "<<currentnodeforinsertion<<" "<<currentnodeforinsertion->GetLabel()<<endl;
        if (chosennodenumberforinsertion==NodesTouched) { // We've found the node to split on
			//if (currentnode!=Root) {
            NodeToAdd=currentnodeforinsertion;
            //      selectedanode=true;
            break;
            //  }
            //  else {
            //     chosennodenumber++; //Don't want to insert on the root node
            //     cout<<"Chosen node was root; choosing next instead"<<endl;
            //  }
        }
		
        currentnodeforinsertion = m.next();
    }
    //if (selectedanode==false) {
    //  cout<<"Error -- didn't select a node for addition"<<endl<<endl;
    // }
    // cout<<"NodeToAdd = "<<NodeToAdd<<" "<<NodeToAdd->GetLabel()<<endl;
    //  cout<<"\nDraw Line 209\n";
	
    //    Draw(cout);
    NodePtr NewPosition=AddNodeToSubtendingEdge(NodeToAdd);
    //cout<<"\nDraw Line 212\n";
    //Draw(cout);
    (NewPosition->GetChild())->SetSibling(NodeToBreak); //We know NewPosition has one child with no siblings as a result of how we created New Position
    NodeToBreak->SetAnc(NewPosition);
    NewPosition->IncrementDegree();
    FindAndSetRoot();
    Update();
    //cout<<"Draw Line 218\n";
    //Draw(cout);
}

void ContainingTree::TBR () {
}

vector<int> ContainingTree::SplitLeaf (int LeafToSplit) {
	// cout<<"\nStarting tree for split leaf"<<endl;
	// ReportTreeHealth();
    vector<int> output(2,0);
    int startingnumberofleaves=GetNumLeaves();
    output[0]=LeafToSplit;
    SetLeafNumbers();
    NodeIterator <Node> m (GetRoot());
    NodePtr currentnode = m.begin();
    NodePtr SelectedLeafPtr;
    bool success=false;
    while (currentnode)
    {
        if (currentnode->IsLeaf()) {
            // cout<<"currentnode->GetLeafNumber() = "<<currentnode->GetLeafNumber()<<" LeafToSplit = "<<LeafToSplit<<endl;
            if (currentnode->GetLeafNumber()==LeafToSplit) {
                SelectedLeafPtr=currentnode;
                success=true;
                break;
            }
        }
        currentnode = m.next();
    }
    //  if (!success) {
    //      cout<<"Tried to match "<<LeafToSplit<<" but failed."<<endl;
    //      ReportTreeHealth();
    //  }
    assert(success);
    NodePtr q;
    if (GetNumLeaves()==1) {
        SuppressInternalNodesWithOneDescendant();
        //cout<<"\n\n{{{{{{{{{{{{{{{   NOW WORKING WITH SINGLE LEAF TREE }}}}}}}}}}}}}}}"<<endl;
        SelectedLeafPtr->SetAnc(NULL);
        SelectedLeafPtr->SetChild(NULL);
        SelectedLeafPtr->SetLeaf(false);
        nxsstring NewLabel="";
        SelectedLeafPtr->SetLabel(NewLabel);
        CurNode=SelectedLeafPtr;
    }
    else {
		// ReportTreeHealth();
		// cout<<endl;
        NodePtr Ancestor = SelectedLeafPtr->GetAnc();
        if (SelectedLeafPtr->IsTheChild())
        {
            Ancestor->SetChild (SelectedLeafPtr->GetSibling());
            q = SelectedLeafPtr->GetSibling ();
        }
        else
        {
            q = SelectedLeafPtr->LeftSiblingOf ();
            q->SetSibling (SelectedLeafPtr->GetSibling ());
        }
        SelectedLeafPtr->SetSibling (NULL);
        SelectedLeafPtr->SetAnc (NULL); //So we've deleted the node that speciates
        Ancestor->SetDegree (Ancestor->GetDegree() - 1);
        CurNode=Ancestor->GetChild();
		// ReportTreeHealth();
        MakeSibling(); //We replace the leaf that we've just deleted with an internal node (that's sister to its ancestor's daughter).
    }
    nxsstring NewLeaf1="taxon";
    NewLeaf1+=LeafToSplit;
    nxsstring NewLeaf2="taxon";
    NewLeaf2+=startingnumberofleaves+1;
    AddCherry(NewLeaf1,NewLeaf2);
    output[1]=startingnumberofleaves+1;
	//  cout<<"\nInitial health"<<endl;
	//   ReportTreeHealth();
    Update();
    SetLeafNumbers();
	// cout<<"\nFinal health"<<endl;
	// ReportTreeHealth();
    return output;
}

vector<int> ContainingTree::CollapseCherry (int chosencherrynum) {
    //this will involve changing the vector assigning gene samples to species
    //Note that initially, this will probably make the score much worse
    vector<int> output(2,0);
    SetLeafNumbers();
    UpdateCherries();
    NodePtr ChosenCherry=CherryNodes[chosencherrynum];
    NodePtr NodeA=ChosenCherry->GetChild();
    NodePtr NodeB=NodeA->GetSibling();
    //Draw(cout);
    // cout<<"\nNodeA = "<<NodeA->GetLabel()<<"\nNodeB = "<<NodeB->GetLabel()<<endl;
    int NodeANum=NodeA->GetLeafNumber();
    int NodeBNum=NodeB->GetLeafNumber();
    if (NodeANum>NodeBNum) { //We kill the one with the higher taxon number
        output[0]=NodeBNum;
        output[1]=NodeANum;
        DeleteSubtendingEdge(NodeA);
        delete NodeA;
        FixLeafNamesAndNumbers(NodeANum);
    }
    else {
        output[0]=NodeANum;
        output[1]=NodeBNum;
        DeleteSubtendingEdge(NodeB);
        delete NodeB;
        FixLeafNamesAndNumbers(NodeBNum);
    }
    FindAndSetRoot();
    Update();
    SetLeafNumbers();
    UpdateCherries();
    return output;
}

void  ContainingTree::DeleteSubtendingEdge(NodePtr SelectedNode) {
    //cout<<"Now in DeleteSubtendingEdge"<<endl;
    if (SelectedNode->GetAnc()!=NULL) {
        NodePtr AncestralNode = SelectedNode->GetAnc();
        if (AncestralNode->GetChild()==SelectedNode) { //Selected Node is a child
            NodePtr SiblingOfSelectedNode=SelectedNode->GetSibling();
            SelectedNode->SetSibling(NULL);
            SelectedNode->SetAnc(NULL);
            AncestralNode->SetChild(SiblingOfSelectedNode);
            AncestralNode->SetDegree(-1+(AncestralNode->GetDegree()));
        }
        else { //Selected node is a sibling of something else
            NodePtr LinkingNode=AncestralNode->GetChild();
            while (LinkingNode->GetSibling()!=SelectedNode) {
                LinkingNode=LinkingNode->GetSibling();
            }
            LinkingNode->SetSibling(SelectedNode->GetSibling()); //So, move the sibiling link over one (even if it's null)
            SelectedNode->SetSibling(NULL);
            SelectedNode->SetAnc(NULL);
            AncestralNode->SetDegree(-1+(AncestralNode->GetDegree()));
        }
		//cout<<"deleted subtending, but not suppressed internal nodes yet\n\n";
		//ReportTreeHealth();
        SuppressInternalNodesWithOneDescendant();
		//cout<<"\nNow suppression has happened\n\n";
    }
}


void ContainingTree::TestRerooting()
{
	
    cout<<"Original Tree"<<endl;
    NodeIterator <Node> m (GetRoot());
    NodePtr currentnodem = m.begin();
    int currentcount=0;
    while (currentnodem)
    {
        currentcount++;
        nxsstring newlabel="";
        newlabel+=currentcount;
		
        currentnodem->SetLabel(newlabel);
        currentnodem = m.next();
    }
	
    Draw(cout);
    cout<<endl;
    //   Write(cout);
    //   cout<<endl;
    //   for (int SPRnum=1;SPRnum<50;SPRnum++) {
    //        cout<<"After "<<SPRnum<<" SPR"<<endl;
    //        SPR();
    //        ReportTreeHealth();
    //        cout<<endl;
    //       Draw(cout);
    //       cout<<endl;
    //       Write(cout);
    //      cout<<endl;
    //      cout<<"Tree has "<<GetNumLeaves()<<" leaves and "<<GetNumInternals()<<" internals"<<endl;
    //  }
    int NodesTouched=0;
    nxsstring NewLeafLabel;
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
		
        NodesTouched++;
        if (NodesTouched==7) {
            NodePtr NewRoot= ReRootTree(currentnode);
            break;
        }
        currentnode = n.next();
    }
}



//This will reroot a tree on the branch subtending the current node and return a pointer to the new root
NodePtr ContainingTree::ReRootTree(NodePtr SelectedNode)
{
    //ReportTreeHealth();
    NodePtr OldRoot = GetRoot();
    NodePtr ThisNode;
    // cout<<"OldRoot ="<<OldRoot->GetLabel()<<endl;
    //  cout<<"SelectedNode = "<<SelectedNode->GetLabel()<<endl;
	
    //Need to change things on path between new root and old root, change sibling of selectednode if it has one, change thing whose sibling is selected node if it exists, and delete old root node
    NodePtr NewRoot=AddNodeToSubtendingEdge(SelectedNode); //adds a degree 2 node below the selected node
	//nxsstring Newlabel="NEWROOT";
	//NewRoot->SetLabel(Newlabel);
	//  cout<<"NewRoot = "<<NewRoot->GetLabel()<<endl;
	// NodePtr TempNewNode=NewNode();
	// TempNewNode->SetAnc(NewRoot);
	//  (NewRoot->GetChild())->SetSibling(TempNewNode);
	// TempNewNode->SetLeaf(true);
	// TempNewNode->SetLabel("newrootpos");
	// Update();
	// Draw(cout);
    (NewRoot->GetChild())->SetSibling(NewRoot->GetAnc());
    NodePtr TaxonToMakeSibNULLLater=Root->GetChild();
    vector<NodePtr> PathToOldRoot;
    NodePtr CurrentNode=NewRoot->GetAnc();
    PathToOldRoot.push_back(NewRoot);
    while (CurrentNode!=NULL) {
        //cout<<"CurrentNode = "<<CurrentNode->GetLabel()<<endl;
        //ThisNode=CurrentNode;
        //NodePtr PreviousNode=PathToOldRoot.back();
        //cout<<"PreviousNode = "<<PreviousNode->GetLabel()<<endl;
        PathToOldRoot.push_back(CurrentNode);
        CurrentNode=CurrentNode->GetAnc();
        //ThisNode->SetAnc(PreviousNode);
        //(PreviousNode->GetChild())->SetSibling(ThisNode);
        //ThisNode->SetChild(CurrentNode);
    }
    // cout<<"finished while loop"<<endl;
    for (int i=1;i<PathToOldRoot.size();i++) {
        //    cout<<"On node "<<PathToOldRoot[i]->GetLabel()<<endl;
        (PathToOldRoot[i])->SetAnc(PathToOldRoot[i-1]);
        // cout<<"  to "<<((PathToOldRoot[i])->GetAnc())->GetLabel()<<endl;
        //  cout<<"\tNode "<<((PathToOldRoot[i-1])->GetChild())->GetLabel()<<" sib changed from "<<(((PathToOldRoot[i-1])->GetChild())->GetSibling())->GetLabel();
        ((PathToOldRoot[i-1])->GetChild())->SetSibling((PathToOldRoot[i]));
        //  cout<<" to "<<(((PathToOldRoot[i-1])->GetChild())->GetSibling())->GetLabel()<<endl;
        NodePtr NewChild=( (PathToOldRoot[i])->GetChild());
        if (NewChild==(PathToOldRoot[i-1])) {
            NewChild=NewChild->GetSibling();
        }
        //  cout<<"\tChild changed from "<<((PathToOldRoot[i])->GetChild())->GetLabel();
        (PathToOldRoot[i])->SetChild(NewChild);
    }
    NewRoot->SetSibling(NULL);
    NewRoot->SetAnc(NULL);
    TaxonToMakeSibNULLLater->SetSibling(NULL);
    Root=NewRoot;
    SetRoot(NewRoot);
    //  cout<<endl;
    //  ReportTreeHealth();
    //   cout<<endl;
    // ReportTreeHealth();
    NodePtr NodeToCheck=NewRoot;
    while (NodeToCheck!=NULL) {
        NodeToCheck->SetSibling(NULL);
        if ((NodeToCheck->GetChild())!=NULL) {
            NodeToCheck=(NodeToCheck->GetChild())->GetSibling();//a sibling should have no siblings with binary trees
        }
        else {
            NodeToCheck=NULL;
        }
    }
    SuppressInternalNodesWithOneDescendant();
    // ReportTreeHealth();
    //Draw(cout);
    FindAndSetRoot();
    Update();
	
    GetNodeDepths();
    //Draw(cout);
    // ReportTreeHealth();
    return NewRoot;
	
}

NodePtr ContainingTree::AddNodeToSubtendingEdge(NodePtr NodeWhereAttach) {
    assert(NodeWhereAttach!=NULL);
    NodePtr AddedNode=NewNode();
    NodePtr AncestralNode=NodeWhereAttach->GetAnc();
    NodePtr DescendantNode=NodeWhereAttach;
    if (AncestralNode==NULL) { //Means we're inserting below the root node
        AddedNode->SetAnc(NULL);
        SetRoot(AddedNode);
        AddedNode->SetSibling(DescendantNode->GetSibling());
        DescendantNode->SetSibling(NULL);
        DescendantNode->SetAnc(AddedNode);
        AddedNode->SetChild(DescendantNode);
        AddedNode->SetDegree(1);
        AddedNode->SetDepth(DescendantNode->GetDepth());
        DescendantNode->SetDepth(1+(DescendantNode->GetDepth()));
    }
    else {
        if (AncestralNode->GetChild()==DescendantNode) { //So, node to substitute below is a child node
            AddedNode->SetAnc(AncestralNode);
            AncestralNode->SetChild(AddedNode);
            AddedNode->SetSibling(DescendantNode->GetSibling());
            DescendantNode->SetSibling(NULL);
            DescendantNode->SetAnc(AddedNode);
            AddedNode->SetChild(DescendantNode);
            AddedNode->SetDegree(1);
            AddedNode->SetDepth(DescendantNode->GetDepth());
            DescendantNode->SetDepth(1+(DescendantNode->GetDepth()));
        }
        else { //Node to substitute below is a sibling of some other node
            NodePtr ChildOfAncestor=AncestralNode->GetChild();
            //Now must move over until we get the sibling closest to the node of interest
            while ((ChildOfAncestor->GetSibling())!=DescendantNode) {
                ChildOfAncestor=ChildOfAncestor->GetSibling();
            }
            ChildOfAncestor->SetSibling(AddedNode);
            AddedNode->SetAnc(AncestralNode);
            AddedNode->SetSibling(DescendantNode->GetSibling());
            DescendantNode->SetSibling(NULL);
            DescendantNode->SetAnc(AddedNode);
            AddedNode->SetChild(DescendantNode);
            AddedNode->SetDegree(1);
            AddedNode->SetDepth(DescendantNode->GetDepth());
            DescendantNode->SetDepth(1+(DescendantNode->GetDepth()));
        }
    }
    return AddedNode;
}

void ContainingTree::SuppressInternalNodesWithOneDescendant()
{
    //cout<<"Now in SuppressInternalNodesWithOneDescendant"<<endl;
    Update();
    vector <NodePtr> NodesToDelete;
    PreorderIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    NodePtr SelectedNode;
    //  NodePtr AncestralNode;
    //  NodePtr SiblingOfSelectedNode;
    //  NodePtr DescendantNode;
    while (currentnode)
    {
        //cout<<"currentnode "<<currentnode<<" degree "<<currentnode->GetDegree();
        if ((currentnode->GetDegree()==1))
        {
            NodesToDelete.push_back(currentnode);
			// cout<<" <----Want to delete "<<currentnode<<" "<<currentnode->GetLabel();
        }
		// cout<<endl;
        currentnode = n.next();
    }
    while (!NodesToDelete.empty()) {
        SelectedNode=NodesToDelete.back(); //Pulls a node to delete
        NodesToDelete.pop_back(); //And deletes it from the stack
		//cout<<"Now deleting node "<<SelectedNode<<" within Suppression setp"<<endl;
        if (SelectedNode!=NULL) {
            SuppressInternalNode(SelectedNode);
        }
    }
    FindAndSetRoot();
    Update();
}

void ContainingTree::SuppressInternalNode(NodePtr SelectedNode)
{
    //ReportTreeHealth();
    //cout<<"\nSelectedNode is "<<SelectedNode<<endl;
    //cout<<"Now in SuppressInternalNode"<<endl;
    assert(!(SelectedNode->IsLeaf()));
	
    NodePtr AncestralNode = SelectedNode->GetAnc();
    if (AncestralNode!=NULL) {
        NodePtr DescendantNode = SelectedNode->GetChild();
        if (AncestralNode->GetChild()==SelectedNode) { //Selected Node is a child
            AncestralNode->SetChild(DescendantNode);
            DescendantNode->SetAnc(AncestralNode);
			NodePtr SibNode=DescendantNode->GetSibling();
			if (SibNode==NULL) { //covers case where we're suppressing a node with one descendant, so the descendant has no sibs
				DescendantNode->SetSibling(SelectedNode->GetSibling());
			}
			while (SibNode!=NULL) {
				SibNode->SetAnc(AncestralNode);
				if((SibNode->GetSibling())==NULL) {
					SibNode->SetSibling(SelectedNode->GetSibling());
					break;
				}
				SibNode=SibNode->GetSibling();
			}
/*			if ((SelectedNode->GetSibling())!=NULL) {
				(SelectedNode->GetSibling())->SetAnc(AncestralNode);
				(SelectedNode->GetSibling())->SetSibling(DescendantNode->GetSibling());
				DescendantNode->SetSibling(SelectedNode->GetSibling());
			}			
			NodePtr SibNode=DescendantNode->GetSibling();
			while (SibNode!=NULL) {
				SibNode->SetAnc(AncestralNode);
				SibNode=SibNode->GetSibling();
			} */
            SelectedNode->SetSibling(NULL);
            SelectedNode->SetAnc(NULL);
        }
        else { //Selected node is a sibling of something else
            NodePtr LinkingNode=AncestralNode->GetChild();
			while (LinkingNode->GetSibling()!=SelectedNode) { //only necessary if dealing with non-binary trees
                LinkingNode=LinkingNode->GetSibling();
			}
            LinkingNode->SetSibling(DescendantNode); //So, move the sibiling link over one (even if it's null)
			
            DescendantNode->SetAnc(AncestralNode);
			NodePtr SibNode=DescendantNode->GetSibling();
			if (SibNode==NULL) { //covers case where we're suppressing a node with one descendant, so the descendant has no sibs
				DescendantNode->SetSibling(SelectedNode->GetSibling());
			}			
			while (SibNode!=NULL) {
				SibNode->SetAnc(AncestralNode);
				if((SibNode->GetSibling())==NULL) {
					SibNode->SetSibling(SelectedNode->GetSibling());
					break;
				}
				SibNode=SibNode->GetSibling();
			}
			
		/*	if ((SelectedNode->GetSibling())!=NULL) {
				(SelectedNode->GetSibling())->SetAnc(AncestralNode);
				(SelectedNode->GetSibling())->SetSibling(DescendantNode->GetSibling());
				DescendantNode->SetSibling(SelectedNode->GetSibling());
			}			
			NodePtr SibNode=DescendantNode->GetSibling();
			while (SibNode!=NULL) {
				SibNode->SetAnc(AncestralNode);
				SibNode=SibNode->GetSibling();
			} */
            SelectedNode->SetSibling(NULL);
            SelectedNode->SetAnc(NULL);
			SelectedNode->SetChild(NULL);
        }
        delete SelectedNode;
        FindAndSetRoot();
        Update();
    }
    else { //So all that's left on the tree is one taxon and the root. Kill the root
        NodePtr DescendantNode = SelectedNode->GetChild();
        DescendantNode->SetAnc(NULL); //because it's now the root
        SetRoot(DescendantNode);
        SelectedNode->SetSibling(NULL);
        SelectedNode->SetAnc(NULL);
        SelectedNode->SetChild(NULL);
        delete SelectedNode;
        Update();
    }
}

int ContainingTree::TraverseToGetNodeCount (NodePtr p)
{
    nodecount=0;
    if (p)
    {
        traverse (p->GetChild());
        nodecount++;
        traverse (p->GetSibling());
        nodecount++;
    }
    return nodecount;
}

//Note that I use a relaxed cherry definition: A cherry is a node that is directly connected to two OR MORE leaves (i.e., a polytomy of three tip taxa is a "cherry"). For trees of 2 or more taxa, there is always at least one, and perhaps as many as ntax/2, cherries.
void ContainingTree::UpdateCherries()
{
    CherryNodes.clear();
    if(GetNumLeaves()>1) {
        // ReportTreeHealth();
        PreorderIterator <Node> n (GetRoot());
        bool IsCherry=false;
        NodePtr currentnode = n.begin();
        //  cout<<"currentnode is "<<currentnode;
        // cout<<" "<<currentnode->GetLabel()<<endl;
        NodePtr sibnode;
        while (currentnode)
        {
            if ((currentnode->GetChild())!=NULL) {
                if ((currentnode->GetChild())->IsLeaf())
                {
                    IsCherry=true;
                    sibnode=(currentnode->GetChild())->GetSibling();
                    if (sibnode->IsLeaf()) {
                        while ((sibnode->GetSibling())!=NULL) {
                            sibnode=sibnode->GetSibling();
                            if (!(sibnode->IsLeaf())) {
                                IsCherry=false;
                            }
							
                        }
                    }
                    else {
                        IsCherry=false;
                    }
                    if (IsCherry) {
                        CherryNodes.push_back(currentnode);
                        //  cout<<"currentnode "<<currentnode<<" is a cherry\n";
                    }
                }
            }
            currentnode = n.next();
        }
    }
}

void ContainingTree::FindAndSetRoot()
{
    NodePtr currentnode = GetRoot();
    //cout<<"currentnode is "<<currentnode<<" with label "<<currentnode->GetLabel()<<endl;
    while ((currentnode->GetAnc())!=NULL) {
        //cout<<"FindAndSetRoot: currentnode="<<currentnode<<endl;
        assert(currentnode!=(currentnode->GetAnc()));
        currentnode=currentnode->GetAnc();
    }
    SetRoot(currentnode);
    Root=currentnode;
    // cout<<"Final root = "<<Root<<endl;
}

vector<ContainingTree> ContainingTree::SplitOnTaxon (vector<nxsstring> taxatoexclude)
{
	
    //ignore the following comment -- since we now assume rooted gene trees, don't have to reroot
    //So, we want to split whereever there's a connection between a subtree consisting only of the chosen taxa and a subtree consisting only of the excluded taxa. problem is rooting: what if tree were ((C1(C2(C3(E1(E2(C4(C5))))))),C6), with the root between C6 and the rest of the tree? We want ((C1,C2),C3) and ((C6,C5),C4), rooted on the node that connects each subtree to the rest of the tree. Solution: first reroot somewhere between the E taxa.
    vector<ContainingTree> SplitTreeVector;
    //cout<<"Before splitting, gene tree is "<<endl;
    Update();
    FindAndSetRoot();
    GetNodeDepths();
    //Draw(cout);
    //cout<<endl;
    MarkNodes(false);
    vector<NodePtr> FoundLeafVector;
    MakeNodeList(); //needed before finding leaves
    for (int i=0; i<taxatoexclude.size(); i++) {
        nxsstring taxontofind=taxatoexclude[i];
        // cout<<"Looking for taxon "<<taxontofind;
        NodePtr FoundLeaf=GetLeafWithLabel(taxontofind);
        if (FoundLeaf!=NULL) {
            FoundLeafVector.push_back(FoundLeaf);
            //  cout<<" ... found it";
        }
        //  cout<<endl;
    }
    if (FoundLeafVector.size()>0) {
        //cout<<"Now getting subtrees...\n";
        //NodePtr Placeholder=ReRootTree(FoundLeafVector[0]); //just reroot it on one of the excluded leaves
        //since we assume rooted trees now, we shouldn't reroot
        //ReportTreeHealth();
        for (int i=0; i<FoundLeafVector.size(); i++) {
            NodePtr FoundLeaf=FoundLeafVector[i];
            while ((FoundLeaf!=NULL) && ((FoundLeaf->IsMarked())==false)) { //if it's already marked, all its descendants will be, too
				// cout<<"Now examining "<<FoundLeaf<<" "<<FoundLeaf->GetLabel();
                FoundLeaf->SetMarked(true);
                FoundLeaf=FoundLeaf->GetAnc();
                // cout<<" moving on to "<<FoundLeaf<<endl;
            }
        }
        //so, all the nodes connected to deleted taxa are marked true
        NodeIterator <Node> n (GetRoot());
        NodePtr currentnode = n.begin();
        while (currentnode)
        {
            if (currentnode!=Root) {
                //   if ((currentnode->IsMarked())==false) {
                //     cout<<"Node "<<currentnode<<": "<<currentnode->GetLabel()<<" is not marked"<<endl;
                //  }
                if (((currentnode->IsMarked())==false) && (((currentnode->GetAnc())->IsMarked())==true)) {
                    ContainingTree ChosenSubtree;
                    ChosenSubtree.SetRoot(CopyOfSubtree(currentnode));
                    // cout<<"\nThis subtree now"<<endl;
                    ChosenSubtree.FindAndSetRoot();
                    ChosenSubtree.Update();
                    ChosenSubtree.GetNodeDepths();
                    // ChosenSubtree.Draw(cout);
                    ChosenSubtree.MarkNodes(false);
                    SplitTreeVector.push_back(ChosenSubtree);
                }
            }
            currentnode = n.next();
        }
    }
    else { //there was nothing to exclude
        ContainingTree FinalTree;
        FinalTree.SetRoot(CopyOfSubtree(Root));
        SplitTreeVector.push_back(FinalTree);
    }
    return SplitTreeVector;
}

void ContainingTree::FixLeafNamesAndNumbers(int deletedleafnumber)
{
    //this only works if taxa are named taxon1, taxon17, etc.
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if(currentnode->IsLeaf()) {
            int speciesnumber;
            nxsstring specieslabel=currentnode->GetLabel();
            string speciesstring=specieslabel.c_str();
            size_t index = speciesstring.find("taxon",0);
            if (index!=string::npos) {
                speciesstring.erase(index,5); //erase "taxon"
                speciesnumber=atoi(speciesstring.c_str());
                if (speciesnumber>deletedleafnumber) {
                    //cout<<"Need to rename taxon"<<speciesnumber<<" to taxon"<<speciesnumber-1<<endl;
                    currentnode->SetLeafNumber(speciesnumber-1);
                    currentnode->SetLabelNumber(speciesnumber-1);
                    nxsstring newspecieslabel="taxon";
                    newspecieslabel+=(speciesnumber-1);
                    currentnode->SetLabel(newspecieslabel);
                }
            }
        }
        currentnode = n.next();
    }
    //cout<<"After FixLeafNamesAndNumbers\n";
    //Draw(cout);
}

void ContainingTree::SetLeafNumbers()
{
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    int leafnodecount=0;
    while (currentnode)
    {
        if(currentnode->IsLeaf()) {
            int speciesnumber;
            leafnodecount++;
            nxsstring specieslabel=currentnode->GetLabel();
            string speciesstring=specieslabel.c_str();
            size_t index = speciesstring.find("taxon",0);
            if (index!=string::npos) {
                speciesstring.erase(index,5); //erase "taxon"
                speciesnumber=atoi(speciesstring.c_str());
                currentnode->SetLeafNumber(speciesnumber);
                currentnode->SetLabelNumber(speciesnumber);
            }
            else {
                currentnode->SetLeafNumber(leafnodecount);
                currentnode->SetLabelNumber(leafnodecount);
            }
        }
        currentnode = n.next();
    }
}

vector<int> ContainingTree::GetLeafNumberVector()
{
    vector<int> LeafNumberVect;
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if(currentnode->IsLeaf()) {
            LeafNumberVect.push_back(currentnode->GetLeafNumber());
        }
        currentnode = n.next();
    }
    return LeafNumberVect;
}

vector<nxsstring> ContainingTree::GetLeafLabelVector()
{
    vector<nxsstring> LeafLabelVector;
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if(currentnode->IsLeaf()) {
            LeafLabelVector.push_back(currentnode->GetLabel());
        }
        currentnode = n.next();
    }
    return LeafLabelVector;
}

vector<int> ContainingTree::GetLCADepthVector(nxsstring a, nxsstring b, nxsstring c)
{
    //first, initialize depth vector
    PreorderIterator <Node> n (GetRoot());
    int count = 0;
    Node *q = n.begin();
    while (q)
    {
        depth[q] = count++;
        q = n.next();
    }
	
    //Now get depths for each pair
    vector<int> LCADepthVector;
    //MakeNodeList(); //Note: you must call MakeNodeList before using this
    LCADepthVector.push_back(LCADepthQuery(a,b));
    LCADepthVector.push_back(LCADepthQuery(b,c));
    LCADepthVector.push_back(LCADepthQuery(a,c));
    return LCADepthVector;
}

int ContainingTree::LCADepthQuery(nxsstring a, nxsstring b) {
    //make sure you've called MakeNodeList before using this
    NodePtr p = GetLeafWithLabel(a);
    NodePtr q = GetLeafWithLabel(b);
	//  cout<<"Comparing "<<p->GetLabel()<<" and "<<q->GetLabel();
    while (depth[p] != depth[q])
    {
        if (depth[p] < depth[q])
            q = q->GetAnc();
        else
            p = p->GetAnc();
    }
	//   cout<<" depth is "<<depth[p]<<endl;
    return depth[p];
}

NodePtr ContainingTree::GetLeafWithNumber (int i)
{
    //make sure you've called MakeNodeList before using this
    NodePtr result = NULL;
    // cout<<"result is "<<result;
    nxsstring LeafLabel="taxon";
    LeafLabel+=i;
    result=GetLeafWithLabel(LeafLabel.c_str());
    // cout<<" now "<<result;
    if (result==NULL) { //means we didn't find it using the default name
        NodeIterator <Node> n (GetRoot());
        NodePtr currentnode = n.begin();
        while (currentnode) {
            // cout<<"current node is "<<currentnode;
            if(currentnode->IsLeaf()) {
                //    cout<<" "<<currentnode->GetLabel();
                if (currentnode->GetLeafNumber()==i) {
                    result=currentnode;
                    break;
                }
            }
            //  cout<<endl;
            currentnode = n.next();
        }
		
    }
    // cout<<" finally "<<result<<endl;
	
    return result;
}

vector<int> ContainingTree::GetPotentialNewRoots()
{
    //ReportTreeHealth();
    vector<int> TempNodesToReRootOn;
    vector<int> NodesToReRootOn;
    int currentnodenum=0;
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if(currentnode!=Root) {
            if((currentnode->GetAnc())!=Root) {
                TempNodesToReRootOn.push_back(currentnodenum);
            }
        }
        currentnode = n.next();
        currentnodenum++;
    }
    if (TempNodesToReRootOn.size()>0) { //since we need at least three taxa in the species tree to have an interesting rerooting, there may be no nodes to reroot on
        gsl_permutation * p = gsl_permutation_alloc (TempNodesToReRootOn.size());
        gsl_permutation_init (p);
        gsl_ran_shuffle (r, p->data, TempNodesToReRootOn.size(), sizeof(size_t));
        for (int i=0; i<TempNodesToReRootOn.size(); i++) {
            NodesToReRootOn.push_back(TempNodesToReRootOn[gsl_permutation_get (p,i)]);
        }
        gsl_permutation_free(p);
    }
    return NodesToReRootOn;
}

NodePtr ContainingTree::SelectNodeToReRootOn(int chosennodenum)
{
    int currentnodenum=0;
    NodePtr NodeToReRootOn=Root;
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if (currentnodenum==chosennodenum) {
            NodeToReRootOn=currentnode;
            break;
        }
        currentnode = n.next();
        currentnodenum++;
    }
    return NodeToReRootOn;
}

void ContainingTree::ClearInternalLabels()
{
    NodeIterator <Node> n (GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if (!(currentnode->IsLeaf())) {
            currentnode->SetLabel("");
        }
        currentnode = n.next();
    }
}

void ContainingTree::InitializeMissingBranchLengths()
{
	FindAndSetRoot();
	Update();
	GetNodeDepths();
	//cout<<"InitializeMissingBranchLengths"<<endl;
	//ReportTreeHealth();
	SetEdgeLengths(true);
	MarkNodes(false); //We'll mark the nodes below which the branch lengths have been changed
	getPathLengths(Root); //make sure we store the path lengths
	GetNodeDepths(); //heights above root
	float maxpathlength = GetMaxPathLength();
	if (maxpathlength==0) { //there are no path lengths stored; initialize tree with all brlen = 1
		NodeIterator <Node> n (GetRoot());
		NodePtr currentnode = n.begin();
		while (currentnode)
		{
			if (currentnode!=Root) {
				currentnode->SetEdgeLength(1.0*currentnode->GetDepth() - 1.0*(currentnode->GetAnc())->GetDepth());
				//cout<<endl<<currentnode<<" "<<currentnode->GetEdgeLength();
			}
			/*			if (((currentnode->GetEdgeLength())==0) || ((currentnode->GetEdgeLength())!=(currentnode->GetEdgeLength())) || gsl_isnan(currentnode->GetEdgeLength()) || gsl_isinf(currentnode->GetEdgeLength())) { //test for ==0 or ==nan
			 double totalbrlen=0;
			 double countofbranches=0;
			 if (currentnode->GetAnc()!=NULL) {
			 if ((((currentnode->GetAnc())->GetEdgeLength())==0) || (((currentnode->GetAnc())->GetEdgeLength())!=((currentnode->GetAnc())->GetEdgeLength())) || gsl_isnan((currentnode->GetAnc())->GetEdgeLength()) || gsl_isinf((currentnode->GetAnc())->GetEdgeLength())) {
			 }
			 else {
			 totalbrlen+=(currentnode->GetAnc())->GetEdgeLength();
			 countofbranches++;
			 }
			 }
			 if (currentnode->GetSibling()!=NULL) {
			 if ((((currentnode->GetSibling())->GetEdgeLength())==0) || (((currentnode->GetSibling())->GetEdgeLength())!=((currentnode->GetSibling())->GetEdgeLength())) || gsl_isnan((currentnode->GetSibling())->GetEdgeLength()) || gsl_isinf((currentnode->GetSibling())->GetEdgeLength())) {
			 }
			 else {
			 totalbrlen+=(currentnode->GetSibling())->GetEdgeLength();
			 countofbranches++;
			 }
			 }
			 if (currentnode->GetChild()!=NULL) {
			 if ((((currentnode->GetChild())->GetEdgeLength())==0) || (((currentnode->GetChild())->GetEdgeLength())!=((currentnode->GetChild())->GetEdgeLength())) || gsl_isnan((currentnode->GetChild())->GetEdgeLength()) || gsl_isinf((currentnode->GetChild())->GetEdgeLength())) {
			 }
			 else {
			 totalbrlen+=(currentnode->GetChild())->GetEdgeLength();
			 countofbranches++;
			 }
			 }
			 double newlength=1.0;
			 if (countofbranches>0) {
			 newlength=totalbrlen/countofbranches; //assume missing brlen similar to that of neighbors
			 }
			 currentnode->SetEdgeLength(newlength); */
			currentnode = n.next();
		}
		
	}
	
	else { //we have to be more clever
		PreorderIterator <Node> n (GetRoot());
		NodePtr currentnode = n.begin();
		while (currentnode)
		{
			if (currentnode==GetRoot()) {
				currentnode->SetEdgeLength(0);
				getPathLengths(currentnode);
			}
			else if (((currentnode->GetEdgeLength())==0)) {
				float pathlength = currentnode->GetPathLength();
				if (currentnode->IsLeaf()) { //is a leaf; we've done its ancestors
					currentnode->SetEdgeLength(maxpathlength-pathlength);
					currentnode->SetMarked(true);
					getPathLengths(currentnode);
				}
				else {
					currentnode->SetEdgeLength(0.5*(maxpathlength-pathlength));
					currentnode->SetMarked(true);
					getPathLengths(currentnode);
				}
			}
			currentnode = n.next();
		}
		
		
	}
}
	
void ContainingTree::RandomlyModifySingleBranchLength(double markedmultiplier, double brlensigma)
{
	FindAndSetRoot();
	Update();
	GetNodeDepths();	
	SetEdgeLengths(true);
	//cout<<"RandomlyModifySingleBranchLength"<<endl;
	//ReportTreeHealth();
	int numberofunmarkednodes=-2+(2*GetNumLeaves());
	int numberofmarkednodes=0;
	NodeIterator <Node> n (GetRoot());
	NodePtr currentnode = n.begin();
	//cout<<"numberofunmarkednodes = "<<numberofunmarkednodes<<" numberofmarkednodes = "<<numberofmarkednodes<<endl;
	while (currentnode)
	{
		if (currentnode->IsMarked()) {
			numberofunmarkednodes--;
			numberofmarkednodes++;
			//cout<<"numberofunmarkednodes = "<<numberofunmarkednodes<<" numberofmarkednodes = "<<numberofmarkednodes<<endl;
		}
		currentnode = n.next();
	}
	double probabilityperunmarked=1.0/(numberofunmarkednodes+markedmultiplier*numberofmarkednodes);
	bool changedbrlen=false;
	while (!changedbrlen) {
		currentnode = n.begin();
		while (currentnode)
		{
			double thresholdprobability=probabilityperunmarked;
			if (currentnode->IsMarked()) {
				thresholdprobability*=markedmultiplier;
			}
			double testvalue=gsl_ran_flat (r,0,1);
			//cout<<"Testvalue = "<<testvalue<<" thresholdprobability = "<<thresholdprobability<<endl;
			if (testvalue<thresholdprobability) { //adjust brlen
				//currentnode->SetEdgeLength(gsl_ran_lognormal(r,currentnode->GetEdgeLength(),brlensigma)); //Original brlen moving, but resulted in too high values
				//use gamma distribution
				//with gamma, mean = a * b and var = a * b * b
				//we want mean = currentnode->GetEdgeLength() and var = brlensigma^2
				//so solve for best values of a and b
				double a=(currentnode->GetEdgeLength())*(currentnode->GetEdgeLength())/(brlensigma*brlensigma); //a = x^2/ var
				double b=(brlensigma*brlensigma)/(currentnode->GetEdgeLength()); //b = var / x
				currentnode->SetEdgeLength(gsl_ran_gamma(r,a,b));
				changedbrlen=true;
				//cout<<"Should break here"<<endl;
				break;
			}
			currentnode = n.next();
		}
		
	}
	//cout<<"Done adjusting"<<endl;
}

void ContainingTree::NodeSlideBranchLength(double markedmultiplier)
{
	//Find an internal, nonroot node, and move it up or down, adjusting length of its two child edges, too
	FindAndSetRoot();
	Update();
	GetNodeDepths();	
	SetEdgeLengths(true);
	if (GetNumLeaves()>2) { //so there are internal nodes to move
		//cout<<"RandomlyModifySingleBranchLength"<<endl;
		//ReportTreeHealth();
		int numberofunmarkednodes=GetNumLeaves()-2; //number of internal nodes
		int numberofmarkednodes=0;
		NodeIterator <Node> n (GetRoot());
		NodePtr currentnode = n.begin();
		//cout<<"numberofunmarkednodes = "<<numberofunmarkednodes<<" numberofmarkednodes = "<<numberofmarkednodes<<endl;
		while (currentnode)
		{
			if (currentnode->IsMarked()) {
				if (!(currentnode->IsLeaf()) && currentnode!=GetRoot()) {
					numberofunmarkednodes--;
					numberofmarkednodes++;
				}
			}
			currentnode = n.next();
		}
		double probabilityperunmarked=1.0/(numberofunmarkednodes+markedmultiplier*numberofmarkednodes);
		bool changedbrlen=false;
		while (!changedbrlen) {
			currentnode = n.begin();
			while (currentnode)
			{
				if (!(currentnode->IsLeaf()) && currentnode!=GetRoot()) {
					double thresholdprobability=probabilityperunmarked;
					if (currentnode->IsMarked()) {
						thresholdprobability*=markedmultiplier;
					}
					double testvalue=gsl_ran_flat (r,0,1);
					//cout<<"Testvalue = "<<testvalue<<" thresholdprobability = "<<thresholdprobability<<endl;
					if (testvalue<thresholdprobability) { //adjust brlen
						double belowedgelength=currentnode->GetEdgeLength();
						
						double aboveedgelength=0.9*belowedgelength;
						if ((currentnode->GetSibling())!=NULL) { //has sibling
							aboveedgelength=(currentnode->GetSibling())->GetEdgeLength();
						}
						else { //is the rightmost sibling, so get its left sib
							aboveedgelength=((currentnode->GetAnc())->GetChild())->GetEdgeLength();
						}
						double movedistance=gsl_ran_flat(r,-0.5*(GSL_MIN(belowedgelength,aboveedgelength)),0.5*(GSL_MIN(belowedgelength,aboveedgelength)));
						currentnode->SetEdgeLength(belowedgelength-movedistance);
						NodePtr nextnode=currentnode->GetChild();
						while (nextnode!=NULL) {
							nextnode->SetEdgeLength(aboveedgelength+movedistance);
							nextnode=nextnode->GetSibling(); //so, doesn't assume binary tree
						}
						changedbrlen=true;
						//cout<<"Should break here"<<endl;
						break;
					}
				}
				currentnode = n.next();
			}
			
		}
	}
	//cout<<"Done adjusting"<<endl;
	else {
		cout<<"Warning: attempting a move where there are no internal nodes (besides the root)"<<endl;
	}
}

void ContainingTree::ModifyTotalBranchLength(double brlensigma)
{
	FindAndSetRoot();
	Update();
	GetNodeDepths();	
	SetEdgeLengths(true);
	double roottotipheight=0;
	NodePtr currentnode=GetRoot();
	while (currentnode!=NULL) {
		roottotipheight+=currentnode->GetEdgeLength();
		currentnode=currentnode->GetChild();
	}
	//use gamma distribution
	//with gamma, mean = a * b and var = a * b * b
	double a=(roottotipheight)*(roottotipheight)/(brlensigma*brlensigma); //a = x^2/ var
	double b=(brlensigma*brlensigma)/(roottotipheight); //b = var / x
	double multiplier=(gsl_ran_gamma(r,a,b))/roottotipheight; //we'll scale whole tree
	
	NodeIterator <Node> n (GetRoot());
	currentnode = n.begin();
	currentnode = n.begin();
	while (currentnode)
	{
		currentnode->SetEdgeLength((currentnode->GetEdgeLength())*multiplier);
		currentnode = n.next();
	}
	//cout<<"Done adjusting"<<endl;
}
