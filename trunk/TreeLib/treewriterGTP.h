/*
 * TreeLib
 * A library for manipulating phylogenetic trees.
 * Copyright (C) 2001 Roderic D. M. Page <r.page@bio.gla.ac.uk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307, USA.
 */
 
// $Id: treewriter.h,v 1.1.1.1 2006/05/21 15:38:03 bcomeara Exp $

//MODIFIED BY BCO to output to Sanderson's GTP

#ifndef TREEWRITERGTP_H
#define TREEWRITERGTP_H

#include "TreeLib.h"
#include "nodeiterator.h"

#include <iostream>

/**
 * @class NewickTreeWriterGTP
 * Base class for writing Newick format tree descriptions.
 *
 */
class NewickTreeWriterGTP
{
public:
	NewickTreeWriterGTP (Tree *tree) { t = tree; writeEdgeLengths = true; };
        virtual ~NewickTreeWriterGTP () {};
	
	/**
	 * Toggle on or off the writing of edge lengths in the tree description.
	 */	
	virtual void SetWriteEdgeLengths (bool on) { writeEdgeLengths = on; };
	/**
	 * Set output stream to which the tree description is written
	 * @param s pointer to the stream
	 */	
	virtual void SetStream (std::ostream *s) { f = s; };
	/**
	 * Write the tree description.
	 */	    
    virtual void WriteGTP ();
protected:
	std::ostream *f;
	Node *cur;
	std::stack < Node *, std::vector<Node *> > stk;
	Tree * t;
	
	/**
	 * Write the symbol signalling the end of the tree description. By default this is ';'
	 */	
	virtual void WriteEndOfTree ();
	/**
	 * Write '('
	 */	
    virtual void WriteLeftParenthesis ();
	/**
	 * Write ')'
	 */	
    virtual void WriteRightParenthesis ();
	/**
	 * Write the sibling symbol. By default this is ','
	 */	
    virtual void WriteSiblingSymbol ();
	/**
	 * Write leaf node.
	 */	
    virtual void WriteLeaf ();
	/**
	 * Write internal node.
	 */	
    virtual void WriteInternal ();
    
    bool writeEdgeLengths;

};

#endif


