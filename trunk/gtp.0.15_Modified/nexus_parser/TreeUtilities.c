NODETYPE * sister(NODETYPE * n)
{
if (n->anc == NULL)
	return NULL;
if (n->anc->firstdesc==n)
	return n->sib;
else
	return prevSib(n);
}


NODETYPE *createSubtree(NODETYPE *srcNode, int SubtreeSize)

/* Returns a pointer to a newly allocated tree, which is created by copying a 
subtree from tree srcRoot. Copies pertinent time information from source nodes.
Each node also stores a pointer to the node on the source tree from whence it came.
This permits rapid updating of information about time.

At the moment this routine ignores SubtreeSize, and makes a subtree from three branches
surrounding srcNode.

*/

{
NODETYPE *root, *cnode,*node,*child;
if (!isTip(srcNode) && !isRoot(srcNode)) /* only allow subtrees from internal nodes! */
	{
	root=newnode();
	copyNodeInfo(srcNode->anc,root);
	root->nodePtr=srcNode->anc;
	cnode=newnode();
	AddChild(root,cnode);
	copyNodeInfo(srcNode,cnode);
	cnode->nodePtr=srcNode;
	child=srcNode->firstdesc;
	SIBLOOP(child)  /* for each child of the source node, create a child on the copied tree */
		{
		node=newnode();
		AddChild(cnode,node);
		copyNodeInfo(child,node);
		node->nodePtr=child;		
		}
	return root;
	}
else
	return NULL;
}

/**********************************/



void AddChild(NODETYPE * parent, NODETYPE * theChild)
        {
	NODETYPE *aChild;
	if (parent->firstdesc)
	    {
	    aChild=parent->firstdesc;
	    if (aChild)
		    {
		    while(aChild->sib)
			    aChild=aChild->sib;
		    aChild->sib=theChild;
		    }
	    }
	else
	    parent->firstdesc=theChild;
	theChild->anc=parent;
	theChild->sib=NULL;
        return;
        }
void RemoveTaxonLvAnc(NODETYPE * rmTaxon)

/* remove a tip or clade, but leave its ancestor node in place */

{
NODETYPE * prev;
prev=prevSib(rmTaxon);
if(prev)	/* either rmTaxon is the firstdesc or its got a prev sib */
	prev->sib=rmTaxon->sib;
else
	rmTaxon->anc->firstdesc=rmTaxon->sib;
}


NODETYPE * RemoveTaxon(TREE T,NODETYPE * rmTaxon)

/* Removes a taxon, or clade, including the stem lineage
 * Does not remove the node below the stem lineage if that becomes of degree two
 * Does not deallocate memory for the subtree that is deleted, or change links on that subtree.
 * Won't allow removal of the root node
 * If the node is one of only two children of the root, the root is removed as well.
 * RETURNS A POINTER TO THE PRUNED TREE!
 */
     {
     NODETYPE *n, *prev, *parent,*grandparent,*sis,*root;
     if (T)
     	root=T->root;
     else
	root=NULL;	/* used for cases in which don't know the tree and don't care */
     if (rmTaxon==NULL)
	return root;
     if (!isRoot(rmTaxon))
        {
	parent=rmTaxon->anc;
	grandparent=parent->anc; /* might be NULL if parent is the root */
	if (node_tomy(parent)==2)
		{
		sis=sister(rmTaxon);
		if (isRoot(parent)) 
			{
			sis->anc=NULL;	/* make sure this node acquires 'root' status */
			sis->sib=NULL;
			return sis; /* new root of tree is this sister node */
			}
		else
			{
			sis->anc=grandparent;
			prev=prevSib(parent);
			if (prev)
				{
				prev->sib=sis;
				}
			else
				{
				grandparent->firstdesc=sis;
				}
			sis->sib=parent->sib;
			sis->length+=parent->length;
			}
		}
	if (node_tomy(parent)>2)
		{
		prev=prevSib(rmTaxon);
		if(prev)	/* either rmTaxon is the firstdesc or its got a prev sib */
			prev->sib=rmTaxon->sib;
		else
			parent->firstdesc=rmTaxon->sib;
		}

        return root;
        }
     }

NODETYPE * prevSib(NODETYPE* node)

/* returns the sib that points to this sib, or null if this sib is the first desc or if this sib is root */

	{
	NODETYPE *prev, *n;
	prev=NULL;
	if(!isRoot(node))
	    {
	    n=node->anc->firstdesc;
	    while(n != node)  
			    {
			    if (n->sib == NULL)
				     return NULL;   
			    prev=n;
			    n=n->sib;
			    }
	    return prev;
	    }
	else
	    return NULL;
	}

NODETYPE * ReRoot(NODETYPE * atNode)
    {
	NODETYPE *n, *r;
	if(isRoot(atNode))
	    return atNode; /* don't change the root */
	n=atNode->anc;
	if (!isRoot(n))
	    {
	    r=newnode();
	    RemoveTaxon(NULL,atNode);
	    AddChild(r, atNode);
	    Flip(n);
	    AddChild(r, n);
	    n->length=0; /* leave all the length on the left root's branch */
	    init_node_ids(r, 0);
	    return r;
	    }
	else
	    return n; /* don't change the root here either */
    }
void Flip(NODETYPE *a)
    {
	NODETYPE * b,  *saveAnc;
	float saveLength;
	b=a->anc;	
	if (!isRoot(b))
	    {
	    Flip(b);  /* recurse until the root, then back up */		
	    }
	RemoveTaxon(NULL,a);
	AddChild(a, b);
	b->length=a->length; /* flip the branch lengths too */
	if (node_tomy(b)==1)  /* then delete this node */
	    {
	    saveLength=b->length;
	    RemoveTaxon(NULL,b);
	    AddChild(a, b->firstdesc);
	    b->firstdesc->length+=saveLength;
	    Node_Destructor(b);
	    /* deallocate node b HERE */
	    };

	return;
    }

/***********************************************************************************/
void DisposeTree(NODETYPE *node)

	/* Frees up the tree memory and its taxon names */
{
	NODETYPE *child;
	if (!node) return;
	child=node->firstdesc;
	SIBLOOP(child) {
			DisposeTree(child);
			}
	Node_Destructor(node);
	return;
}

/***********************************************************************************/
void copyNodeInfo(NODETYPE *source,NODETYPE *dest)

/* Copies SOME information about one node to another node */

{

	dest->taxon_name=DupStr(source->taxon_name);
	dest->length=source->length;
	dest->time=source->time;
	dest->minAge=source->minAge;
	dest->maxAge=source->maxAge;
	dest->id=source->id;
	dest->free=source->free;
	dest->numdesc=source->numdesc;
	dest->estRate=source->estRate;
	return;
}


/***********************************************************************************/
		root=makegroup();		
		}
	else 
		return(NULL);
	return(root);
}

/***********************************************************************************/
NODETYPE * find_id(NODETYPE *node,int id)
/* returns the node that has an id that matches 'id' or NULL
if it fails to find */


{
	NODETYPE *child, *found_node=NULL;

	if (node->id == id)
			return node;
	child=node->firstdesc;
	SIBLOOP(child) 
		if (found_node=find_id(child,id) )
			return found_node;
	return NULL;
}

/***********************************************************************************/
void TreeToTaxaPtrList(NODETYPE *node,  PtrList NodeList)	
{

/* Moves through clade from node, compiling list of terminals NODES (!);
on entry taxaList must be a valid pointer to a possibly empty list */

	NODETYPE *child;
	if (isTip(node)) 
		{
		pListAddItem(NodeList, node);
		return;
		}

	else	/* interior node */
		{
		child=node->firstdesc;
		SIBLOOP(child)
			{
			TreeToTaxaPtrList(child, NodeList);
			}
		return;
		}
}
