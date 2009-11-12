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
 
 // $Id: treewriter.cpp,v 1.1.1.1 2006/05/21 15:38:03 bcomeara Exp $

//MODIFIED BY BCO to output to Sanderson's GTP

#include "treewriterGTP.h"

//------------------------------------------------------------------------------
void NewickTreeWriterGTP::WriteLeftParenthesis ()
{
	*f << '(';
}

//------------------------------------------------------------------------------
void NewickTreeWriterGTP::WriteRightParenthesis ()
{
	*f << ')';
}

//------------------------------------------------------------------------------
void NewickTreeWriterGTP::WriteSiblingSymbol ()
{
	*f << ',';
}

//------------------------------------------------------------------------------
void NewickTreeWriterGTP::WriteLeaf ()
{
	*f << NEXUSString (cur->GetLabel());
	//if (t->GetHasEdgeLengths () && writeEdgeLengths)
	//	*f << ":" << cur->GetEdgeLength();
}


//------------------------------------------------------------------------------
void NewickTreeWriterGTP::WriteInternal ()
{
	if (cur->GetLabel() != "")
		*f << NEXUSString (cur->GetLabel());
	//if (t->GetHasEdgeLengths () && writeEdgeLengths)
	//	*f << ":" << cur->GetEdgeLength();
}

//------------------------------------------------------------------------------
void NewickTreeWriterGTP::WriteEndOfTree ()
{
	*f << ';';
}


//------------------------------------------------------------------------------
void NewickTreeWriterGTP::WriteGTP ()
{
	cur = t->GetRoot();

    while (cur)
    {
        if (cur->GetChild())
        {
            WriteLeftParenthesis ();
            stk.push (cur);
            cur = cur->GetChild();
        }
        else
        {
            WriteLeaf ();
            while (!stk.empty() && (cur->GetSibling() == NULL))
            {
                WriteRightParenthesis ();
                cur = stk.top();
                WriteInternal ();
                stk.pop();
            }
            if (stk.empty())
                cur = NULL;
            else
            {
                WriteSiblingSymbol ();
                cur = cur->GetSibling();
            }
        }
    }
    WriteEndOfTree ();
}
