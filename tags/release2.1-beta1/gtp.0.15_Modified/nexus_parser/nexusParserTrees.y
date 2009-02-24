%{
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mySmallTreeLib.h"
#include "my_vector.h"
#include "my_slist.h"
#include "my_hash.h"
#include "my_structures.h" // for DupStr

char cbuf[64];
extern node gRootNode;
extern Vector translVec; 
extern Vector treesVec;
extern int translationTable;

slistEntry se;
slist taxList=NULL;
Hash  mrcaHash;
Entry e;
%}


%token PUNCT NEXUS NXBEGIN END TREE QUOTED
%token  INT  REAL ALPHANUM TRANSLATE MRCA 

%union {
	char sval[1000];
	struct nodetype *nodeval;
	struct slistType *sList;	
	long intval;
	double flval;
	}
%type <nodeval> treeroot leaf tree siblist
%type <sval> identifier taxon_name ALPHANUM number INT REAL

%%
nexus: NEXUS blocklist
	|nexus NEXUS blocklist;
	;
blocklist: block
	|blocklist block
	;
block:	NXBEGIN identifier ';' treelist END ';'
	;
treelist: treeroot   {vectorPushBack(treesVec,$1);}
	| treelist treeroot {vectorPushBack(treesVec,$2);}
	;
treeroot: TREE identifier '=' tree ';' {$$=$4; /*printf("Processed tree command\n"); printtree($4);*/}
	;



siblist: leaf 					{$$=$1;}
	| tree 					{$$=$1;}
	| siblist ',' leaf 			{appendSib($1,$3);$$=$1;}
	| siblist ',' tree			{appendSib($1,$3);$$=$1;}
	;

/* Note that the following is context sensitive, depending on whether the TRANSLATE command has been encountered...
   Does this violate something key about the grammar rules? */

leaf: 	taxon_name		{$$=newnode($1,0.0);}
	| taxon_name ':' number		{$$=newnode($1,strtod($3,NULL));}
	;

tree: '(' siblist ')'  				{$$=makeAnc($2,NULL,0.0);}
	| '(' siblist ')' ':' number  		{$$=makeAnc($2,NULL,strtod($5,NULL));}
	| '(' siblist ')' taxon_name 			{$$=makeAnc($2,$4,0.0);}
	| '(' siblist ')' taxon_name ':' number 	{$$=makeAnc($2,$4,strtod($6,NULL));}
	;


number: INT	{strcpy($$,$1);}	
	| REAL	{strcpy($$,$1);}	
	;

taxon_name: identifier		{strcpy($$,$1);}

identifier: ALPHANUM		{strcpy($$,$1);}
	| INT			{;}
	| QUOTED		{;}
	;

%%

yyerror(char *s)
{
	//fprintf(stderr,"%s\n",s);
}

