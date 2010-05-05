#include "nexusdefs.h"
#include "discretedatum.h"
#include "discretematrix.h"
#include "nexustoken.h"
#include "nexus.h"
#include "taxablock.h"
#include "charactersblock.h"
#include "charactersblock2.h"

/*This class is identical to the DataBlock class
  It was created so that a single nexus file could have two
  character blocks loaded
Brian O'Meara
April 2007
*/

/**
 * @class      CharactersBlock2
 * @file       charactersblock2.h
 * @file       charactersblock2.cpp
 * @author     Paul O. Lewis
 * @copyright  Copyright © 1999. All Rights Reserved.
 * @see        Association
 * @see        AssocList
 * @see        CharactersBlock
 * @see        DiscreteDatum
 * @see        DiscreteMatrix
 * @see        LabelList
 * @see        LabelListAssoc
 * @see        LabelListBag
 * @see        Nexus
 * @see        NexusBlock
 * @see        NexusReader
 * @see        NexusToken
 * @see        SetReader
 * @see        TaxaBlock
 *
 * This class handles reading and storage for the Nexus block DATA.
 * It is derived from the CharactersBlock class, and differs from
 * CharactersBlock only in name and the fact that newtaxa is initially
 * true rather than false.
 */

/**
 * @constructor
 *
 * Performs the following initializations:
 * <table>
 * <tr><th align="left">Variable <th> <th align="left"> Initial Value
 * <tr><td> id             <td>= <td> "DATA"
 * <tr><td> newtaxa        <td>= <td> true
 * </table>
 */
CharactersBlock2::CharactersBlock2( TaxaBlock& tb, AssumptionsBlock& ab )
	: CharactersBlock( tb, ab )
{
	id = "CHARACTERS2";
	newtaxa = false;
}

/**
 * @method Reset [void:protected]
 *
 * Calls Reset function of the parent class (CharactersBlock) and
 * resets newtaxa to true in preparation for reading another DATA block.
 */
void CharactersBlock2::Reset()
{
	CharactersBlock::Reset();
}


