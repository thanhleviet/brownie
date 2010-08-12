#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mySmallTreeLib.h"
#include "my_vector.h"
#include "my_slist.h"
#include "my_hash.h"
#include "my_structures.h"
#include <math.h>

Vector translVec;
Vector treesVec;
int translationTable=0;
node gRootNode;
node gn1;
long gId;
static void initNodeIdsHelper(node n);
static void initNodeIdsHelper2(node n);

void initNodeIds2(node n, int id)
{
gId=id;
initNodeIdsHelper2(n);
}
static void initNodeIdsHelper2(node n)
{
	node child;
	if (!n) 
		return;
	n->id2=gId++;
	child=n->firstdesc;
	SIBLOOP(child) 
		initNodeIdsHelper2(child);
	return;
}
void initNodeIds(node n, int id)
{
gId=id;
initNodeIdsHelper(n);
}
static void initNodeIdsHelper(node n)
{
	node child;
	if (!n) 
		return;
	n->id=gId++;
	child=n->firstdesc;
	SIBLOOP(child) 
		initNodeIdsHelper(child);
	return;
}


int isBinaryTree(node n)

// it's binary if each interior node has exactly two children

{
	int retval=1;
	node child,sibling;
	child=n->firstdesc;
	if (!isTip(n)) 
		{
		sibling=child->sib;
		if (sibling == NULL) return 0; // node has only one kid --> not binary
		else
			if (sibling->sib) return 0; // node has three or more kids --> not binary	
		SIBLOOP(child) 
			if (!isBinaryTree(child)) return 0;
		}
	return retval;


}

nodeArray newTipNodeArray(node root)
// OK as long as this function gets setup only once; else static init in helper will screw up on subseq calls
{
long N;
nodeArray A;
N=numdesc(root);
A=(nodeArray)malloc(N*sizeof(node));
newTipNodeArrayHelper(A,root);
return A;
}
static void newTipNodeArrayHelper(nodeArray A, node n)
{
node child;
static long ix=0;
if (isTip(n))
	A[ix++]=n;
child = n->firstdesc;
SIBLOOP(child)
	newTipNodeArrayHelper(A,child);
}

Vector nexus2TreesVec(char *buffer)

// Function for use with a parser that expects possibly multiple trees in a block with possibly different label sets

{
extern char* myinputptr;
treesVec=vectorNew(100);
myinputptr=buffer;
yyparse();
return treesVec;
}

node nexus2rootNode(char *buffer)

// As of now, we configure to read from a file throughout...[but see the following note if you want to change]

// Depending on how nexusLexer.l is configured, this routine either invokes a parser on STDIN or invokes a parser
// based on the string pointed to by buffer. STDIN or string must contain a valid nexus file with a tree block and one
// tree. It is parsed and a pointer to the root node of a tree structure is returned.

{
extern char* myinputptr;
translVec=vectorNew(100);
myinputptr=buffer;
yyparse();
return gRootNode;
}

#define	INIT_BUFFER_SIZE	10000 

/* reads from a file stream a nexus file and returns the string buffer */

char * slurpNexus (FILE * inFileStream)

{
	char	*buffer;
	int		c;
	long	count=0,i=0,bufSize=INIT_BUFFER_SIZE;

	
	buffer=(char*)malloc(bufSize*sizeof(char));

	while ((c=fgetc(inFileStream)) != EOF)	/* have to define c as int so that EOF can be detected*/
		{
		if (i >= bufSize-1) /* have to save room for terminating null */
			{
			bufSize*=1.6;
			buffer=(char*)realloc(buffer,bufSize*sizeof(char));
			if (!buffer)
				printf("Failure to grow buffer in slurpNexus\n");
			}
		buffer[i]=c;
		++i;
		
		}
		buffer[i]='\0';
return buffer;	
}

unsigned long numNodes(node n)
{

	unsigned long sum=0;
	node child;
	child=n->firstdesc;
	SIBLOOP(child) 
		sum += numNodes(child); /* add one for each child and all that children's*/
	return (1+sum);	/* the 1 is to count this node, which must be internal */
}

/***********************************************************************************/
void preOrderVoid(node n,void (*f)(node))
{
	node child;
	(*f)(n);
	if (!isTip(n))
		{
		child=n->firstdesc;
		SIBLOOP(child) 
			preOrderVoid(child,f);
		}
	return ;	
}
double preOrder(node n,double (*f)(node))
{
	double sum=0;
	node child;
	sum+=(*f)(n);
	if (!isTip(n))
		{
		child=n->firstdesc;
		SIBLOOP(child) 
			sum += preOrder(child,f);
		}
	return (sum);	
}
/***********************************************************************************/
unsigned long  maxorderEffective(node n)  
// the maximum number of nodes between this node and a descendant tip
// corrected in case of internal named node
{
        unsigned long  max,temp;
        node child;
        if (!n) return(-1);
        if (isTip(n)  ) 
		{n->order=0; return (0);}
	else
		{
		if (n->label)
			{
			n->order=0;
printf("%s %li\n",n->label,n->order);
			return 0;
			}
        	max=0;
        	child=n->firstdesc;
        	SIBLOOP(child) {
                	        temp=maxorder(child);
                        	if (temp > max) max = temp;
                        	}
        	n->order=max+1;
        	return (max+1);
		}
}
unsigned long  maxorder(node n)  // the maximum number of nodes between this node and a descendant tip
{
        unsigned long  max,temp;
        node child;
        if (!n) return(-1);
        if (isTip(n)  ) {n->order=0; return (0);}
        max=0;
        child=n->firstdesc;
        SIBLOOP(child) {
                        temp=maxorder(child);
                        if (temp > max) max = temp;
                        }
        n->order=max+1;
        return (max+1);
}

unsigned long  numdescEffective(node n)

// For every node: 
//	1. If it is named, set the numdesc = some function of the real numdesc, say sqrt(numdesc)
//	....or, lately, just a fixed value to leave space around it, say 10;
//	2. If it is not named, calculate the number of offspring in the usual way (though noting 
//         that some descendant named nodes might be counted as the fake value.

// Entire tree is set up this way. --Useful for collapsing clades

{
	unsigned long sum=0;
	node child;
	if (!n) return(-1);
	if (isTip(n)) 
		{
		n->numdescEffective=1; 
		return (1);
		}
	child=n->firstdesc;
	SIBLOOP(child) 
		sum+=numdescEffective(child);
	if (n->label)
//		n->numdescEffective=2*sqrt(sum); // just an arbitrary function for looks
		n->numdescEffective=10;
	else
		n->numdescEffective=sum;
//printf("%s %li\n",n->label,n->numdesc);
	return (n->numdescEffective);

}
unsigned long  numdesc(node n)

/* determines the number of leaves descended from every node and stores them at node */

{
	unsigned long sum=0;
	node child;
	if (!n) return(-1);
	if (isTip(n)) 
		{
		n->numdesc=1; 
		return (1);
		}
	child=n->firstdesc;
	SIBLOOP(child) 
		sum+=numdesc(child);
	n->numdesc=sum;
	return (sum);
}

node newnode(char *label, double number)

// Pass label=NULL if don't want to assign a label to this node; otherwise pass the string

	{
	node n;
	n=(node)malloc(sizeof(struct nodetype));
	if (n)
		{
            n->label=NULL; // try this...
		n->anc=NULL;
		n->firstdesc=NULL;
		n->sib=NULL;
		if (label != NULL)
			n->label=DupStr(label);
		n->number=number;
		n->nodeFlag=0;
		n->data=NULL;
		return n;
		}
	else
		return NULL;
	}

void treeFree(node n) // code not checked
/* Frees up the tree memory and its taxon names */
{
    node child;
    node sibling;
    if (!n) return;
   // printf("\ncurrent node has anc %i, firstdesc %i, sib %i, and number %i\n",n->anc,n->firstdesc,n->sib,n->number);
    child=n->firstdesc;
  //  SIBLOOP(child) {
    //    treeFree(child);
  //  }

    //same as deletetraverse from Rod Page's TreeLib
    sibling=n->sib;
    treeFree(child);
    treeFree(sibling);
    nodeFree(n);
    return;
}

void nodeFree(node n)
{
    if (n->label) free(n->label);
    free(n);
// should free the data structure too if its used...
}

/* find the last sib of given node */

node lastSib(node n)
	{
	if (n)
		while(n->sib)
			n=n->sib;
	return n;
	}




/* add a node as the sib of a list of sibs */

void appendSib(node L, node n)
	{
	if (L)
		{
		while(L->sib)
			L=L->sib;
		L->sib=n;
		}
	return;
	}

// insert a node as sister to an existing node
// Return the root of the new tree, which might be different if we insert node as sister to old root

node insertSister(node root, node treenode, node addnode )
	{
	node new,next,prev;
	new = newnode(NULL,0);
	new->firstdesc=addnode;
	addnode->anc=new;
	addnode->sib=treenode; // root
	treenode->anc=new;
	if (treenode == root)
		return new;	
	else
		{
		next=treenode->sib;
		prev=prevSib(treenode);
		treenode->sib=NULL;
		new->anc=treenode->anc;
		new->sib=next;
		if (!prev) // the treenode is in the firstdesc position
			new->anc->firstdesc=new;
		else
			prev->sib=new;
		return root;
		}
	}

node removeNodeAndAnc(node root, node delnode)

// disconnect a node from its parent and delete its parent node if that node now only has one descendant
// (taking care to fix the links in the original tree caused by this deletion)
// Ignore the request if the node is the root.
// Return the root of the modified tree (it might change)

	{
	node parent,sis;
	int first,polytomy;
	if (isRoot(delnode)) return root;
	parent = delnode->anc;
	if (parent->firstdesc==delnode) first=1; else first=0; //is delnode in firstdesc position?
	if (node_tomy(parent)>2) polytomy=1; else polytomy=0;
	RemoveTaxonLvAnc(delnode);
	if (isRoot(parent) && !polytomy)
			root=root->firstdesc; // this will have been set to correct only desc node by removeTaxonLvAnc
	deleteMonoNode(parent);  // this function only deletes the node if it has one descendant	
	return root;
	}

void deleteMonoNode(node n) 

// if a node has only one descendant delete it and connect parent and child; if it is a tip, ignore...
// if it is the root, just delete it and set its child's anc field to NULL

	{
	node parent, child;
	if (isTip(n)) return;
	if (node_tomy(n)==1 )
	    {
	    child = n->firstdesc;
	    if (!isRoot(n))
		{
		parent = n->anc;
		RemoveTaxonLvAnc(n);
		AddChild(parent,child);
		}
	    else
		child->anc=NULL;
	    nodeFree(n);
	    }
	}
/* make a new node that is ancestral to a list of sibs, L */

node makeAnc(node L, char * label, double number)
	{
	node n;
	n=newnode(label,number);
	if (n)
		{
		n->firstdesc=L;
		L->anc=n;
		while (L->sib) /* set all the sibs ancestor */
			{
			L=L->sib;
			L->anc=n;
			}
		return n;
		}
	else
		return NULL;
	}

void printtree(node n)
	{
	node child;
	int sb;
	if (n->sib) sb=n->sib->id2;
	else sb=-1;
	if (isRoot(n))
		printf("Node %s: number=%f id=%li id2=%li sib=%\n",n->label,n->number,n->id,n->id2,sb);
	else
		printf("Node %s: number=%f id=%li id2=%li anc id2=%li sib=%i\n",n->label,n->number,n->id,n->id2, n->anc->id2,sb);
	if (child=n->firstdesc)
	   for (;child;child=child->sib)
		printtree(child);
	return;
	}

node find_taxon_name(node n,char *name)
/* returns the node that has a taxon name that matches 'name' or NULL
if it fails to find */


{
        node child, found_node;

        if (n->label)
                if (isEqual(name,n->label))
                        return n;
        child=n->firstdesc;
        SIBLOOP(child)
                if (found_node=find_taxon_name(child,name) )
                        return found_node;
        return NULL;
}
void subtreeLight(node root, slist taxonList)

// For the subtree induced by the taxonList (extending down to its MRCA), set the nodeFlag for all its nodes

{
node MRCA,p,p1;
slistEntry e;
char * tax;
if (taxonList->numElements < 2)
	fatal("Must have two or more taxa in subtreeLight command\n");
if (!(MRCA = mrca(root,taxonList)))
	fatal("MRCA not found in subtreeLight()\n");
e=taxonList->head;
while (e)
	{
	tax = e->elem;
printf("%s\n",tax);
	p1 = find_taxon_name(root,tax);
	for (p=p1;p!=MRCA;p=p->anc)
		p->nodeFlag=1;
	MRCA->nodeFlag=1;
	e=e->next;
	}
}
node mrca(node root, slist taxonList)

// NOTE. This routine resets each node's nodeFlag to 0 upon exit.
// Also note that in older versions, I destructively changed the input taxonList; now it just traverses it.
// Something to beware of elsewhere!
{

char * tax1, *tax2;
node 	p1, p2, p, theMRCA;
slistEntry e,e1,e2;
if (taxonList->numElements < 2)
	fatal("Must have two or more taxa in mrca command\n");
//tax1 = slistPopLastElem(taxonList);
e1=taxonList->head;
tax1 = e1->elem;
p1 = find_taxon_name(root,tax1);
if (!p1)
	fatal("Taxon name not found on tree in mrca command\n");
theMRCA=p1;
resetFlag(root); // initialize the tree to all 0's
for (p=p1;p;p=p->anc)
        p->nodeFlag=1;

//while (taxonList->numElements)
e2=e1;
while (e2=e2->next)
	{
	//tax2 = slistPopLastElem(taxonList);
	tax2=e2->elem;
	p2 = find_taxon_name(root,tax2);
	if (!p2)
		fatal("Taxon name not found on tree in mrca command\n");
	for (p=p2;p;p=p->anc)
        	if (p->nodeFlag==1)
               		if (isNodeDescendant(theMRCA,p)) 
				{
				theMRCA=p; // if the current MRCA is actually younger than p, make p the MRCA. 
				break;	// go on to the next taxon in list
				}
	}

resetFlag(root); 
return theMRCA;
}

/**********************************************************************/

int isNodeDescendant(node nodeA, node nodeB)
/*
 * Returns 1 if nodeA is the descendant of nodeB or is nodeB; 0 if not (note that A and B might not be either) 
 */

{
node n;
for(n=nodeA;n;n=n->anc)
        /* worst case, terminates when node = NULL at ancestor of root */
        if (n==nodeB) return 1;
return 0;
}

void resetFlag(node n)
{

	node child;
	n->nodeFlag=0;
	child=n->firstdesc;
	SIBLOOP(child) 
		resetFlag(child); 
}

void fatal(char * msg)
{
printf("%s\n",msg);
exit(1);
}


node  copyTree(node a)
// returns a node that is either a tip, or the root of a properly formatted tree, but its ancestor and sibs are undefined
{
node child,first,newfirst,newn,n,prev;
//printf("Numdesc:%i\n",numdesc(a));
newn = newnode(a->label,a->number);
//printf("New Label created:%s\n",newn->label);
//copy....
if(!isTip(a))
	{
	first=a->firstdesc;
	newfirst=copyTree(first);
//printf("Inserted as first descendant Label:%s\n",newfirst->label);
	newn->firstdesc=newfirst;
	newfirst->anc = newn;
	prev=newfirst;
	child=first->sib; // start loop with the second sib in the sib list...
	SIBLOOP(child)
		{
		n = copyTree(child);
//printf("Inserted as sib Label:%s\n",n->label);
		prev->sib=n;
		prev=prev->sib;
		n->anc = newn;
		}
	}
//printf("NumdescNew:%i\n",numdesc(newn));
return  newn;
}

double calcMaxToTip(node n)

/* Calculates maximum distance from root to tip when lengths are available */

{
        double max=0.0,temp,thisLength;
        node  child;
        if (!n) return(0.0);

        if (isRoot(n))
                thisLength=0.0;
        else
                thisLength=n->number;        /* don't add length under the root */
        if (isTip(n))
                return (thisLength);  /* treat a tip and a compact node the same way */
        child=n->firstdesc;
        SIBLOOP(child) {
                        temp=calcMaxToTip(child);
                        if (temp > max) max = temp;
                        }
        return thisLength+max;
}
double treeLength(node n)

/* Sums the branch lengths over tree assuming they are stored in the 'number' field.
   Does not include the length of the branch subtending the root */

{
	node child;
	double dur;
	if (isRoot(n))
		dur=0.0;
	else
		dur=n->number;
	child=n->firstdesc;
	SIBLOOP(child)
		dur+=treeLength(child);
	return dur;
}
void suppressBinaryNodes(node n)
// traverse a tree and remove any internal nodes with only one descendant
// ignore if node is the root
{
	node child,parent;
	child = n->firstdesc;
	if (node_tomy(n)==1 && !isRoot(n))
		{
		printf("Found a binary node\n");
		parent = n->anc;
		RemoveTaxonLvAnc(n);
		AddChild(parent,child);
		nodeFree(n);
		printf("Freeing binary node\n");
		}
	SIBLOOP(child)
		suppressBinaryNodes(child);		

}
node sister(node n)
// if the parent has just two children, this unambiguously returns the sister group;
// if the parent has more than two children, it returns one of the other children, but which one is anyone's guess...

{
if (n->anc == NULL)
	return NULL;
if (n->anc->firstdesc==n)
	return n->sib;
else
	return prevSib(n);
}

/**********************************/

void AddChild(node parent, node theChild)
        {
	node aChild;
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

void RemoveTaxonLvAnc(node rmTaxon)

/* remove a tip or clade, but leave its ancestor node in place */

{
node prev;
prev=prevSib(rmTaxon);
if(prev)	/* either rmTaxon is the firstdesc or its got a prev sib */
	prev->sib=rmTaxon->sib;
else
	rmTaxon->anc->firstdesc=rmTaxon->sib;
rmTaxon->anc=NULL;
rmTaxon->sib=NULL;
}



node prevSib(node nd)

/* returns the sib that points to this sib, or null if this sib is the first desc or if this sib is root */

	{
	node prev, n;
	prev=NULL;
	if(!isRoot(nd))
	    {
	    n=nd->anc->firstdesc;
	    while(n != nd)  
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

node ReRoot(node root, node atNode)

// at the moment does NOT reset node ids...

    {
	node n, r;
	if(isRoot(atNode))
	    return atNode; /* don't change the root */
	n=atNode->anc;
	if (!isRoot(n))
	    {
	    r=newnode(NULL,0.0);
//printf("ReRoot:RemoveTaxonLvAnc:%i\n",atNode->id2);
	    RemoveTaxonLvAnc(atNode);
//printf("ReRoot:AddChild:parent=%i child=%i\n",r->id2,atNode->id2);
	    AddChild(r, atNode);
//printf("ReRoot:Flip:node=%i\n",n->id2);
	    Flip(n);
//printf("ReRoot:AddChild:parent=%i child=%i\n\n",r->id2,n->id2);
	    AddChild(r, n);
	    return r;
	    }
	else
	    return n; /* don't change the root here either */
    }
void Flip(node x)
    {
	node xanc,parent;
//	printf("Flip:entering node %i\n",x->id2);
	xanc=x->anc;
	if (xanc->anc) 
		{
//		printf("Flip:recursing in ...\n");
		Flip(xanc);
		}
//	else ...this else is a BUG
		{
//printf("Flip:RemoveTaxonLvAnc:%i\n",x->id2);
		RemoveTaxonLvAnc(x);
//printf("Flip:AddChild:parent=%i child=%i\n\n",x->id2,xanc->id2);
		AddChild(x,xanc); // notice this appears to be reversed, as it should be...
		if (node_tomy(xanc) == 1) // A binary node has been added during this procedure...deleting
			{
//printf("Flip:deleting binary node\n");
			parent=xanc->anc; // save this because next call removes link
			RemoveTaxonLvAnc(xanc);
			AddChild(parent,xanc->firstdesc);
			nodeFree(xanc);
			}
		}
    }

/***********************************************************************************/
void copyNodeInfo(node source,node dest)

/* Copies SOME information about one node to another node */

{

	dest->label=DupStr(source->label);
	dest->number=source->number;
	dest->id=source->id;
	return;
}

/***********************************************************************************/
node find_id2(node n,int id)
/* returns the node that has an id that matches 'id' or NULL
if it fails to find */

{
	node child, found_node=NULL;
	if (n->id2 == id)
			return n;
	child=n->firstdesc;
	SIBLOOP(child) 
		if (found_node=find_id2(child,id) )
			return found_node;
	return NULL;
}
/***********************************************************************************/
node find_id(node n,int id)
/* returns the node that has an id that matches 'id' or NULL
if it fails to find */

{
	node child, found_node=NULL;
	if (n->id == id)
			return n;
	child=n->firstdesc;
	SIBLOOP(child) 
		if (found_node=find_id(child,id) )
			return found_node;
	return NULL;
}

#if 0
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
#endif
/***********************************************************************************/
int node_tomy(node n)

/* number of immediate descendants of this node (including this one!) */

{
	node child;
	int tomy=0;
	child=n->firstdesc;
	SIBLOOP(child) 
		++tomy;
	return tomy;
}

void make_parens(node n)

/* writes a parens formatted tree description with labels and durations or
lengths.  flag=0: print lengths; flag =1: print durations as lengths,  
flag=2: print rates as lengths, flag=3: print node id's as lengths */

{
  int width;

  if (isTip(n)) 
    {
    if (n->label)
      printf("%s",n->label);
    }
  else printf("(");

  if (n->firstdesc) make_parens(n->firstdesc);

  if (!isTip(n))
    {
      printf(")");
      if (n->label) 
	    printf("%s",n->label);
    }

  if (n->sib) printf(","),make_parens(n->sib);

}
