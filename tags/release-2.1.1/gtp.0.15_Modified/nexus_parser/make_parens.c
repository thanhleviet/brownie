void make_parens(NODETYPE *node, int flag)

/* writes a parens formatted tree description with labels and durations or
lengths.  flag=0: print lengths; flag =1: print durations as lengths,  
flag=2: print rates as lengths, flag=3: print node id's as lengths */

{
  extern long gNumSites;
  double value, duration;
  int width;
  if (!isRoot(node))
	{
	if (flag==0)
		value = node->length;
	else if (flag == 1)
		value = node->anc->time - node->time; /* duration */
	else if (flag == 2)
		value = node->estRate/gNumSites;		/*rate*/
	}

  if (isTip(node)) 
    {
    if (*(node->taxon_name)=='\0')
	{
	width = log10(node->id)+1; 
	printf("tx%-*i", width, node->id);
	}
    else
      printf("%s",node->taxon_name);
    if (flag == 3)
	printf(":%i",node->id);
    else
    	printf(":%-8.6f",value);
    }
  else printf("(");

  if (node->firstdesc) make_parens(node->firstdesc,flag);

  if (!isTip(node))
    {
      printf(")");
      if (*(node->taxon_name)!='\0') 
	    printf("%s",node->taxon_name);
      if (!isRoot(node)) 
	 {
         if (flag == 3)
		printf(":%i",node->id);
    	 else
	    printf(":%-8.6f",value);
	}
    }

  if (node->sib) printf(","),make_parens(node->sib,flag);

}
