
#include <stdio.h>
#include "my_slist.h"
#include "my_vector.h"
#ifndef _MYSMALLTREE
#define _MYSMALLTREE
#define SIBLOOP(c)			for (; (c); (c)=(c)->sib)
#define isTip(c)			 ( (c)->firstdesc == NULL )
#define isRoot(c)			( (c)->anc == NULL )
#define isEqual(a,b) (!(strcmp((a),(b))))  // currently this is case sensitive!

struct nodetype {
                struct nodetype                 *anc;
                struct nodetype                 *firstdesc;
                struct nodetype                 *sib;
                char                            *label;
                double                          number;
		unsigned long			numdesc;
		unsigned long			numdescEffective; // a different number that is used to guess sizes for layouts when clade collapsed
		unsigned long			order;
		unsigned long			id;
		unsigned long			id2;
		int				nodeFlag;
		int				nodeFlag2;
		void 				*data;
		/*GL*/float*			labelColor;
		/*GL*/float*			treeColor;
		};

typedef         struct nodetype * node;
typedef		node * nodeArray;

void treeFree(node n);
unsigned long  maxorderEffective(node n);  
unsigned long  numdescEffective(node n);
node lastSib(node n);
node  copyTree(node a);
void resetFlag(node n);
node mrca(node root, slist taxonList);
int isNodeDescendant(node nodeA, node nodeB);
void fatal(char *s);
node find_taxon_name(node n,char *name);
char * slurpNexus(FILE *f);
unsigned long  numNodes(node n);
unsigned long  maxorder(node n);
unsigned long  numdesc(node n);
void preOrderVoid(node n,void (*f)(node));
double preOrder(node n,double (*f)(node));
Vector nexus2TreesVec(char *buffer);
node nexus2rootNode(char *buffer);
node newnode(char *label, double number);
void appendSib(node L, node n);
node makeAnc(node L,char * label, double number);
void printtree(node n);
char * DupStr(char *);
double calcMaxToTip(node n);
void subtreeLight(node root, slist taxonList);
double treeLength(node n);
nodeArray newTipNodeArray(node root);
static void newTipNodeArrayHelper(nodeArray A,node child);
void initNodeIds(node n, int id);
void initNodeIds2(node n, int id);
int isBinaryTree(node n);
node prevSib(node n);
node sister(node n);
void AddChild(node parent, node child);
void RemoveTaxonLvAnc(node rmTaxon);
void RemoveTaxon(node rmTaxon);
node ReRoot(node root,node atNode);
void Flip(node a);
void copyNodeInfo(node source,node dest);
int node_tomy(node n);
void nodeFree(node n);
node find_id(node n,int id);
node find_id2(node n,int id);
void suppressBinaryNodes(node n);
void make_parens(node n);
node insertSister(node root, node treenode, node addnode );
node removeNodeAndAnc(node root, node delnode);
void deleteMonoNode(node n);
#endif
