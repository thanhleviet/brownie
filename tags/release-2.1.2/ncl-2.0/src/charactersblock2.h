#ifndef __CHARACTERSBLOCK2_H
#define __CHARACTERSBLOCK2_H

//
// CharactersBlock2 class 
//Just allows us to have a  second characters block in the file
//This added by BCO
//
class CharactersBlock2 : public CharactersBlock
{
	protected:
		void Reset();

	public:
		CharactersBlock2( TaxaBlock& tb, AssumptionsBlock& ab );
};

#endif
