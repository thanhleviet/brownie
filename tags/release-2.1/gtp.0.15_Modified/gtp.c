/*

 This program uses gene tree reconciliation (Goodman et al., 1979; etc.) to calculate duplication scores
 for a given species tree and one or more gene trees. It implements the algorithm described in Zmasek & Eddy (2001),
 which is simple and apparently quite fast, despite its worse-than-linear running time.

 It can be used as the inner code to a full 'gene tree parsimony' search for the most parsimonious species
 tree, as described in Sanderson and McMahon (in review). Hence the name of the program, 'gtp'.

 Input: A nexus file consisting of one or more tree blocks. The entire collection of trees in the file is treated as if it
	were in a single block. The first tree must be the rooted binary species tree, the remaining trees must be rooted
	binary gene trees labeled with the same taxon names as in the species tree. Each tree may be annotated with a tree
	weight if it is part of block of equally parsimonious trees. This weight is added in a non-standard way: by
	adding a nexus-style branch length notation to the root branch of the tree. Note that all EPTs must be listed in
	contiguous blocks in the file. IMPORTANT: the weights, if used, should be entered with 4+ decimals precision to
	ensure that when summed together the entire set is at least ~0.99, even allowing for roundoff error. The program
	relies on this to spot the boundaries between EPT sets. Sorry for the hack.

 Output:

	Detail mode: a table in which each row has the following fields:

 1. The gene tree block number (or gene tree number if there is only one tree in the block)
 2. The number of gene trees in the block, when multiple optimal gene trees present
 3. The total number of duplications assuming the gene tree was rooted as given in the input file
 4. The number of strong duplications assuming the gene tree was rooted...
 5. The minimum number of duplications across all re-rootings of the gene tree
 6. The minimum number of strong duplications across all re-rootings of the gene tree

	In addition, a summary row is printed after this table with the following fields:

 1. Summed value of column 3 in the Detail table
 2. Summed value of column 4 in the Detail table
 3. Summed value of column 5 in the Detail table
 4. Summed value of column 6 in the Detail table
 5. The number of gene trees in which column 3 in the Detail table is exactly 0
 6. The number of gene trees in which column 4 in the Detail table is exactly 0
 7. The number of gene trees in which column 5 in the Detail table is exactly 0
 8. The number of gene trees in which column 6 in the Detail table is exactly 0

	Notes:
	A strong duplication is one that is followed by a speciation event. In addition, in 'unrooted' mode, the score is a
	minimum value across all possible rerootings of the gene tree.

	When multiple equally parsimonious (or likely) trees are present, all scores are weighted
	averages of the scores for individual trees in the EPT set. The collection of such trees
	is called a 'gene tree block'.

 */

/* Program modified by Brian O'Meara July 2006
*
*/


#define VERSION 0.15


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_structures.h"
#include "my_hash.h"
#include "mySmallTreeLib.h"
#include "gtp.h"
#include <math.h>
#include <string.h>

static long numStrongDup(node n);
static long gtp(node n);
static void initAncArray(node n);
static void initgTree(node n, Hash h);
static void makeTipNameHashHelper(node n);
void printUsage(void);
Hash gTipNameHash;
long *gAncArray;

static void initAncArray(node n)
{
    node child;
    if (!isRoot(n))
        gAncArray[n->id]=(n->anc)->id;
    child=n->firstdesc;
    SIBLOOP(child)
        initAncArray(child);
}

Hash makeTipNameHash(node root)
{
    gTipNameHash = hashNew(100);
    makeTipNameHashHelper(root);
    return gTipNameHash;
}
static void makeTipNameHashHelper(node n)
{
    node child;
    Entry e;
    long *id;
    if (isTip(n))
    {
        id = (long*)malloc(sizeof (long));
        *id = n->id;
        hashInsert(gTipNameHash, n->label,id,&e);
    }
    child=n->firstdesc;
    SIBLOOP(child)
        makeTipNameHashHelper(child);
}
static void initgTree(node n, Hash h)
{
    node child;
    Entry e;
    char *key;
    long id;
    if (isTip(n))
    {
        key = n->label;
        e=hashKeyExists(h,key);
        if (e)
        {
            id=*(long*)(e->val);
            n->id=id;
        }
    }
    child=n->firstdesc;
    SIBLOOP(child)
        initgTree(child,h);
}

static long numStrongDup(node n)
{
    node g1,g2;
    int f1,f2;
    long numDup=0;

    if (isTip(n))
    {
        n->nodeFlag=1; // have to set this up for duplication at tip node...careful it doesn't affect other funcs
        n->nodeFlag2=0;
        return 0;
    }
    else
    {
        g1=n->firstdesc;
        g2=g1->sib;
        numDup+=numStrongDup(g1);
        numDup+=numStrongDup(g2);
        if (g1->nodeFlag==0 || g2->nodeFlag==0)
            n->nodeFlag2=1; // obviously if either descendant is a S, then the spec seen flag should be set
        else // both are D, however...
        {
            f1 = g1->nodeFlag2;
            f2 = g2->nodeFlag2;
            if (f1 || f2)
                n->nodeFlag2=1; // one of the desc's nodes has seen a spec event, so set it here too
            else
                n->nodeFlag2=0; // else it's still never been seen
        }
        if (n->nodeFlag2 && n->nodeFlag) ++numDup; // only if a D at this node AND a spec's been seen!
        return numDup;
    }
}

static long strongDup(node n)

// return the number of 'strong' duplications. A strong dup is one that is followed by a speciation event,
// as opposed to weak dups that only occur after last speciation, and might reflect population diff or alleles...

// NB! Using the total number of dups or the number of strong dups as a sp tree search criteria should give the same results!
{
    node child;
    long numDup=0;
    if (isTip(n))
    {
    }
    child=n->firstdesc;
    SIBLOOP(child)
        numDup+=strongDup(child);
    return numDup;
}

static long gtp(node n)
{
    node child,g1,g2; // see Zmasek and Eddy
    long a,b,asave,bsave,numdup=0;
    Entry e;
    char *key;
    long id;
    child=n->firstdesc;
    SIBLOOP(child)
        numdup+=gtp(child);
    if (!isTip(n))
    {
        child=n->firstdesc;
        g1=child;
        g2=child->sib;
        a=g1->id;
        b=g2->id;
        asave=a;
        bsave=b;
        while (a!=b)
        {
            if (a>b)
                a=gAncArray[a];
            else
                b=gAncArray[b];
        }
        n->id=a;
        if (a==asave || a==bsave)
        {
            ++numdup; // Duplication
            n->nodeFlag=1;
        }
        else
            n->nodeFlag=0; // speciation;
    }
    return numdup;
}

//function to return score of number of strong dupes, added by BCO to bypass main and loading files
double ReturnScore(char *buffer,int unrooted) //unrooted=1 if unrooted, 0 otherwise
{
    //printf("\n\ninput is %s\n\n",buffer);
    char * bufferbkup=buffer;
    Vector treesVec;
    Entry e;
    node root,gTreeRoot,spTreeRoot,rroot,rr;
 // char *buffer,*p, buf1[64], buf2[64],*fnInput,*layoutStr,**sargv,*dummy;
    char *p, buf1[64], buf2[64],*fnInput,*layoutStr,*dummy;
	// FILE *f;
    int err=0,c, i,details=0,curEPT=0;
    long id,numRootings,numTrees,lastId;
    double numDup,numSDup,minD,maxD,minSD,maxSD,numDupSave,numSDupSave;
    double grandMinDup=10000,grandMinSDup=10000,grandMinUnrootedDup=10000,grandMinUnrootedSDup=10000;
    double wt, totalD=0,totalSD=0,totalMinD=0,totalMinSD=0,totalMaxD=0,totalMaxSD=0;
    double trackWt;
    double min0,minS0,minU0,minUS0;
    long totalMin0,totalMinS0,numWtTrees=0,totalWtTrees=0;
    long totalUMin0,totalUMinS0;
    Hash spHash;
    //printf("Line 240\n\n");
    node tipNodes[20];
    
 // fnInput=layoutStr=NULL;

#if 0 // temp code for tree building...
    sprintf(buf1,"tax1");
    rr=newnode(buf1,0);
    tipNodes[1]=rr;
    rroot=newnode(NULL,0);
    root=rroot;
    //printf("Line 251\n\n");
    for (i=2;i<=10;i++)
    {
        root=insertSister(root,root,rr);
        sprintf(buf1,"tax%i",i);
        rr=newnode(buf1,0);
        tipNodes[i]=rr;
    }
    printtree(root);
    for (i=1;i<=9;i++)
    {
        root=removeNodeAndAnc(root,tipNodes[i]);
        printtree(root);
    }
    exit(1);

#endif
    //printf("Line 268\n\n");
    // buffer = slurpNexus(f);
    	//printf("Now Printing Buffer:\n%s\nEnd of Buffer",buffer); //Added by BCO
    treesVec=nexus2TreesVec(buffer); // Actually sets up a tree data structure from input file and calls to mySmallTreeLib
    numTrees=vectorLastIndex(treesVec)+1;

  //  int faultloop=0;

   // while ((numTrees < 2) && (faultloop<2)) { //Put this in because there sometimes seems to be an error in the parsing code (it gives a syntax error when first reading a file, but when reading it a second time, works fine
     //   buffer=bufferbkup;
      //  treesVec=nexus2TreesVec(buffer);
      //  numTrees=vectorLastIndex(treesVec)+1;
      //  faultloop++;
    //}
  //  if (faultloop>0) {
  //      printf("               Input to GTP faulted %i times\n",faultloop);
  //  }
    if (numTrees<2)  { //Something went wrong on input
        return -9999999; //a very bad number
                               //printf("Must have at least two input trees");
    }
    else
    {
        //printf("Line 279\n\n");
        for (i=0;i<=numTrees-1;i++)
        {
            if (!isBinaryTree((node)vectorGet(treesVec,i)))
            {
                printf ("Tree %i is not binary (this only works for binary trees) :exiting...\n",i);
                return -9999999; //a very bad number
            }
        }
        if (err) exit (1);
        spTreeRoot=(node)vectorGet(treesVec,0);
        initNodeIds(spTreeRoot, 0);
		//printtree(spTreeRoot);
		//printf("\n\n");
        
        //printf("Line 294\n\n");
        gAncArray = (long *)malloc(2*numdesc(spTreeRoot)*sizeof(long));
        //printf("Line 296\n\n");
        initAncArray(spTreeRoot);
        spHash=makeTipNameHash(spTreeRoot);
        
		// Note hack to weight trees: I use the 'number' field of the root node (obtained in the nexus tree command) to keep the wt
		// We'll assume that the number field of the root is set to 0.0 by default when there's no weight provided!

        trackWt=0.0;
        totalMin0=0;
        totalMinS0=0;
        totalUMin0=0;
        totalUMinS0=0;
        min0=0;
        minS0=0;
        minU0=0;
        minUS0=0;

        for (i=1;i<=numTrees-1;i++) // loop over the provided gene trees
        {
            //printf("Line 315\n\n");
            gTreeRoot=(node)vectorGet(treesVec,i);
            initgTree(gTreeRoot,spHash);
            wt=gTreeRoot->number; // weight the dups by a tree weight if it's stored
            if (wt == 0.0) // hack assumes we never deliberately set a tree weight to zero--otherwise delete from input file stupid
                wt = 1.0;
            numDup=wt*gtp(gTreeRoot);
            numSDup=wt*numStrongDup(gTreeRoot);
            totalD+=numDup;
            totalSD+=numSDup;
            numDupSave=numDup;
            numSDupSave=numSDup;
            min0+=numDup; // this is the total number of dups across the MP trees for this rep, gets reset below
            minS0+=numSDup; // this is the total number of dups across the MP trees for this rep

            rroot=gTreeRoot; // simplifies final destructor below
            if (unrooted)
            {
                numRootings = 2*numdesc(gTreeRoot)-1; 
				//=number of internal nodes on binary rooted tree minus 1 (don't reroot on the root!).
				// (note that some of these rootings are identical; in fact there are only 2n-3 distinct rootings on a binary tree...
                lastId = numRootings -1 ; // and the ids go from 0..numRootings-1
                initNodeIds2(gTreeRoot,0); // these ids will remain in place in field ->id2 as trees are rerooted
                                           // we have to leave ->id field clear because it's used by gtp algorithm
                rroot=gTreeRoot;
                minD=10000;minSD=10000;maxD=-10000;maxSD=-10000;
                for (id=1;id<=lastId;id++) // can skip id=0; that's rerooting at the root!
                {
                    rr=find_id2(rroot,id); // notice when id=0 (the root), this tree will not be rerooted
                    rroot=ReRoot(rroot,rr);
					//make_parens(rroot);printf("\n");
                    //printf("\n\nNow printing on rerooting %i\n\n",id);
                    //printtree(rroot);
                    //printf("\n\n");
                    initgTree(rroot,spHash);
                    numDup=wt*gtp(rroot);
                    numSDup=wt*numStrongDup(rroot);
                    if (numDup < minD) minD=numDup;
                    if (numSDup < minSD) minSD=numSDup;
                } 
            }
            totalMinD+=minD;
            totalMinSD+=minSD;

            
			// some moderately nasty code to keep track of multiple MP trees by using their weights! Note their weights better be there...
			// assumes the weights within an MP set add to 1!
			// min0,minS0,minU0,minUS0 are the weighted summed dup scores in THIS set of weighted trees. For the unrooted versions, they correspond
			// to the minimum scores across the rerootings. (note the rooted versions are stored prior to the unrooting code...)
            minU0+=minD;
            minUS0+=minSD;

            trackWt+=wt;
            ++numWtTrees;
            if (wt == 1.00 || trackWt >= 0.95) // reset
            {
                ++curEPT;
                trackWt=0.0;
                if (min0 == 0.0) totalMin0++; // all the mp trees had min dups = 0
                if (minS0 == 0.0) totalMinS0++; // all the mp trees had smin dups = 0
                if (minU0 == 0.0) totalUMin0++; // all the mp trees had min dups = 0
                if (minUS0 == 0.0) totalUMinS0++; // all the mp trees had smin dups = 0
                totalWtTrees+=numWtTrees; // for a later consistency check...
                numWtTrees=0;
                min0=0;
                minS0=0;
                minU0=0;
                minUS0=0;
            }
            treeFree(rroot);
}
//printf("\nfree treesVec\n");
vectorFree(treesVec);
//printf("free spHash\n");
hashFree(spHash);
//printf("free gAncArray\n");
free (gAncArray);
//printf("free spTreeRoot\n");
treeFree(spTreeRoot);
if (unrooted) {
    return totalMinSD;
}
else {
    return totalSD;
}

	}
}


/***********************************************************/

void printUsage(void)
{

    printf("Usage: gtp [-r] [-u] [-d] -f treefile\n");
    printf("\t-f\tInput file name\n");
    printf("\t-u\tFind min and max num of dupls across all rootings of gene tree\n");
    printf("\t-r\tFind number of duplications assuming gene tree is rooted as is (default)\n");
    printf("\t-d\tPrint details for each gene tree\n");
    printf("\t-v\tPrint version number\n");

    printf("\nNotes:\n");
    printf("1) input file must consist of a nexus tree block with a species tree followed by\n");
    printf(" one or more gene trees, labeled with the same taxon names as are used in the species tree.\n");
    printf("2) to assign optional weights to gene trees, annotate the tree description in the input file\n");
    printf(" by appending the weight as a 'branch length' to the right of the last right parens of the tree.\n");
    printf(" Note that the weights must be given to sufficient decimals precision that their sum in the EPT set\n");
    printf(" is at least equal to 0.99. This often requires 4 or more decimals!\n");
}

//int main(void)
//#{
//    printf("It works!\n");
//    return 0;
//}
