#include <strstream>
#include <fstream.h>
#include <iomanip.h>
#include <unistd.h>
#include <stdio.h>
#include <set>

#include "nexusdefs.h"
#include "xnexus.h"
#include "nexustoken.h"
#include "nexus.h"
#include "taxablock.h"
#include "assumptionsblock.h"
#include "treesblock.h"
#include "discretedatum.h"
#include "discretematrix.h"
#include "charactersblock.h"
#include "gport.h"
#include "profile.h"
#include "nodeiterator.h"
#include "setreader.h"
#include "treeorder.h"
#include "treedrawer.h"
#include "ntree.h"
#include "stree.h"
#include "containingtree.h"
#include "quartet.h"
#include <gsl/gsl_sf_gamma.h>
#include "TreeLib.h"
#include "gtree.h"
#include "treereader.h"
#include "treewriter.h"
#include <time.h>
#include <map>
#include "brownie.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_version.h>
#include "optimizationfn.h"
#include "cdfvectorholder.h"
#include <sstream>
#include <iostream>

//took out this section since GTP is built in
//extern "C" {
//#include "gtp.c"
//    static long gtp(node n);
//    static long strongDup(node n);
//    static long numStrongDup(node n);
//    static void initgTree(node n, Hash h);
//    static void makeTipNameHashHelper(node n);
//    Hash makeTipNameHash(node root);
//    static void initAncArray(node n);
//    double ReturnScore(char *buffer, int unrooted);
//    void printUsage(void);
//}

using std::string;
using std::ostringstream;
gsl_rng * r;  /* global generator */


/**
*  Brownie
 *  Copyright Brian O'Meara, 2006
 *  http://wwww.brianomeara.info/brownie
 *  Released under the GNU Public License
 *  There is no warranty. However, please let
 *  me know about any bugs, and I'll try to
 *  squash them.
 */

/**
* This file is a heavily modified version of the
 * basicccmdline.cpp file that comes with the
 * Nexus class library.
 * Brian O'Meara. February 2006.
 */

/**
* This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */



/**
* @class      BROWNIE
 * @file       brownie.h
 * @file       brownie.cpp
 * @variable   command_maxlen [static int] maximum length of a command entered
 */

/**
* @constructor
 *
 * Initializes the id data member to "BROWNIE" and calls the FactoryDefaults
 * member function to perform the remaining initializations.
 */
BROWNIE::BROWNIE()
{
    id = "BROWNIE";
    FactoryDefaults();
}

/**
* @destructor
 *
 * Closes logf or echof if it is open.
 */
BROWNIE::~BROWNIE()
{
    if( logf_open )
        logf.close();
    if (echof_open) {
        echof<<"end;";
        echof.close();
    }
}

/**
* @method CharLabelToNumber [int:protected]
 * @param s [nxsstring] the character label to be translated to character number
 *
 * The code here is identical to the base class version (simply returns 0),
 * so the code here should either be modified or this derived version
 * eliminated altogether.  Under what circumstances would you need to
 * modify the default code, you ask?  This function should be modified
 * to something meaningful if this derived class needs to construct and
 * run a SetReader object to read a set involving characters.  The SetReader
 * object may need to use this function to look up a character label
 * encountered in the set.  A class that overrides this method should
 * return the character index in the range [1..nchar]; i.e., add one to the
 * 0-offset index.
 */
int BROWNIE::CharLabelToNumber( nxsstring /*s*/ )
{
    return 0;
}

/**
* @method EnteringBlock [virtual void:public]
 * @param blockName [nxsstring] the name of the block just entered
 *
 * Called by the Nexus object when a block named blockName is entered.
 * Allows program to notify user of progress in parsing the NEXUS file.
 * Virtual function that overrides the pure virtual function in the
 * base class Nexus.
 */
void BROWNIE::EnteringBlock( nxsstring blockName )
{
    message = "Reading ";
    message += blockName;
    message += " block...";
    PrintMessage();
}

/**
* @method ExitingBlock [virtual void:public]
 * @param blockName [nxsstring] the name of the block just exited
 *
 * Called by the Nexus object when exiting a block named blockName.
 * Allows program to notify user of progress in parsing the NEXUS file.
 * Virtual function that overrides the pure virtual function in the
 * base class Nexus.
 */
void BROWNIE::ExitingBlock( nxsstring /*blockName*/ )
{
}

/**
* @method FactoryDefaults [void:protected]
 *
 * Sets all data members to their factory default settings:
 * <table>
 * <tr><th align="left">Variable <th> <th align="left"> Initial Value
 * <tr><td> inf_open          <td>= <td> false
 * <tr><td> logf_open         <td>= <td> false
 * <tr><td> quit_now          <td>= <td> false
 * <tr><td> message           <td>= <td> ""
 * <tr><td> next_command[0]   <td>= <td> '\0'
 * <tr><td> trees             <td>= <td> NULL
 * <tr><td> taxa              <td>= <td> NULL
 * <tr><td> assumptions       <td>= <td> NULL
 * <tr><td> characters        <td>= <td> NULL
 * </table>
 */
void BROWNIE::FactoryDefaults()
{
    inf_open = false;
    logf_open = false;
    echof_open=false;
    quit_now = false;
    message = "";
    next_command[0] = '\0';
    gslseedtoprint=time(NULL);
    gslseed=gslseedtoprint;
    gsl_rng_set(r,gslseed);
    trees = NULL;
    taxa = NULL;
    assumptions = NULL;
    characters = NULL;
    chosenchar=1;
    chosentree=1;
    tipvariancetype=0;
    progressbartotal=0;
    progressbarcount=0;
    progressbarprinted=0;
    debugmode=false;
    maxiterations=1000;
    stoppingprecision=1e-4;
    randomstarts=15;
    treefilename="besttrees.tre";
    badgtpcount=0;
    stepsize=.1;
	npercent=0.95;
    detailedoutput=false;
    outputallatonce=true;
    chosenmodel=1;
    citationarray[0]=false;
    maxnumspecies=100;
    minnumspecies=1;
	minsamplesperspecies=3;
	rearrlimit=-1;
    steepest=false;
	exhaustive=false;
    status=true;
	jackknifesearch=false;
    nreps=5;
	jreps=100;
	pctdelete=1.0/3.0;
	jackrep=0;
    showtries=false;
    structwt=0.5;
	triplettoohigh=false;
	gtptoohigh=false;
	tripletdistthreshold=0.2; //Sets how often to use NJ tree distances for starting assignments (higher number=more often) and how often to use triplet support
    pthreshold=1;
	chosensubsampling=2.0;
    movefreqvector.push_back(0.80); //initial values of the movefreqvector; set up so first try doing swaps and rerootings, then reassignments
    movefreqvector.push_back(0.01);
    movefreqvector.push_back(0.01);
    movefreqvector.push_back(0.01);
    movefreqvector.push_back(0.17);
    sppnumfixed=false;
    unrooted=0; //by default, use rooted gene trees
    bestscore=GSL_POSINF;
	bestscorelocal=GSL_POSINF;
    for (int arrayloc=1;arrayloc<20;arrayloc++) {
        citationarray[arrayloc]=false;
    }
    for (int state=0;state<=(maxModelCategoryStates-1);state++) {
        staterestrictionvector.push_back(state);
        timeslicetimes.push_back(GSL_POSINF);
        timeslicemodels.push_back(0);
    }
    //timeslicetimes[maxModelCategoryStates]=GSL_POSINF; //Use the BMC and BMS code; the size of this vector is the maximum number of time intervals minus 1, and the max num of time intervals is the max number of models
    //timeslicemodels[maxModelCategoryStates]=0; //Since we use BMC and BMS code, the size of this vector is the maximum number of model categories
	CDFvectorholder bob;
	CDFvector=bob.Initialize();
}
/**
* @method FileExists [bool:protected]
 * @param fn [const char*] the name of the file to check
 *
 * Returns true if file named fn already exists,
 * false otherwise.
 */
bool BROWNIE::FileExists( const char* fn )
{
    bool exists = false;

    //	if( access( fn, 0 ) == 0 )
    //		exists = true;

    FILE* fp = fopen( fn, "r" );
    if( fp != NULL ) {
        fclose(fp);
        exists = true;
    }

    return exists;
}

/**
* @method GetFileName [nxsstring:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called whenever a file name needs to be read from either
 * the command line or a file.  Expects next token to be "="
 * followed by the token representing the file name.  Call
 * this function after, say, the keyword "file" has been
 * read in the following LOG command:
 * <pre>
 * log file=doofus.txt start replace;
 * </pre>
 * Note that this function will read only "=doofus.txt "
 * leaving "start replace;" in the stream for reading
 * at a later time.
 */
nxsstring BROWNIE::GetFileName( NexusToken& token )
{
    // Eat the equals sign
    //
    token.GetNextToken();

    if( !token.Equals("=") ) {
        errormsg = "Expecting an equals sign, but found ";
        errormsg += token.GetToken();
        errormsg += " instead";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    // Now get the filename itself
    //
    token.GetNextToken();

    return token.GetToken();
}


gsl_rng * BROWNIE::ReturnR()
{
    return r;
}

void BROWNIE::ProgressBar(int total)
{
    if (total>0) { //so to start it, give a value of total >0 (# reps); after that, feed it progressbar(0);
        progressbartotal=total;

        cout<<"\nProgress:\n0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n|"<<flush;
        //cout<<"|....|....|....|....|....|....|....|....|....|....|"
    }
    else {
        progressbarcount++;
        double sampleratio=(1.0*progressbarcount)/(1.0*progressbartotal); // convert to floating point division
        double printratio=progressbarprinted/50.0;
        //cout<<sampleratio<<" "<<printratio<<endl;
        while (sampleratio>printratio) {
            cout<<"*"<<flush;
            progressbarprinted++;
            printratio=(1.0*progressbarprinted)/50.0;
        }
        if (progressbarcount==progressbartotal) { //stop it, reinitialize
            cout<<"|\n\n"<<flush;
            progressbarcount=0;
            progressbartotal=0;
            progressbarprinted=0;
        }
    }
}


/**
* @method GetNumber [nxsstring:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called whenever a number needs to be read from either
 * the command line or a file.  Expects next token to be "="
 * followed by the token representing the number.  Call
 * this function after, say, the keyword "tree" has been
 * read in the following LOG command:
 * <pre>
 * choose tree=5 start replace;
 * </pre>
 * Note that this function will read only "=5 "
 * leaving "start replace;" in the stream for reading
 * at a later time.
 */
nxsstring BROWNIE::GetNumber( NexusToken& token )
{
    // Eat the equals sign
    //
    token.GetNextToken();

    if( !token.Equals("=") ) {
        errormsg = "Expecting an equals sign, but found ";
        errormsg += token.GetToken();
        errormsg += " instead";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    // Now get the number itself
    //
    token.GetNextToken();

    return token.GetToken();
}

/**
* @method GetNumberOnly [nxsstring:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called whenever a number needs to be read from either
 * the command line or a file.  Does not expext an equals
 */
nxsstring BROWNIE::GetNumberOnly( NexusToken& token )
{

    token.GetNextToken();
    return token.GetToken();
}

/**
* @method HandleEndblock [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the END or ENDBLOCK command needs to be parsed
 * from within the BROWNIE block.  Basically just checks to make
 * sure the next token in  the data file is a semicolon.
 */
void BROWNIE::HandleEndblock( NexusToken& token )
{
    // get the semicolon following END or ENDBLOCK token
    //
    token.GetNextToken();

    if( !token.Equals(";") ) {
        errormsg = "Expecting ';' to terminate the END or ENDBLOCK command, but found ";
        errormsg += token.GetToken();
        errormsg += " instead";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }
}

/**
* @method HandleExecuteCmdLine
 * used to execute a datafile from the command line when starting the program
 */
void BROWNIE::HandleExecuteCmdLine(nxsstring fn)
{
    if( FileExists( fn.c_str() ) )
    {
        cerr << endl;
        cerr << "Opening " << fn << "..." << endl;

        ifstream inf( fn.c_str(), ios::binary | ios::in );

        inf_open = true;

        NexusToken ftoken(inf);

        inf_open = true;
        ifstream intreefile;
        intreefile.open(fn.c_str(),ios::in);
        if (!intrees.ReadTrees(intreefile))
        {
            message="No trees read from file\n";
            PrintMessage();
        }
        intreefile.close();

        try {
            Execute( ftoken );
        }
        catch( XNexus x )
        {
            NexusError( errormsg, x.pos, x.line, x.col );
            Reset();
        }



        //  if( !taxa->IsEmpty() ) {
        //      cerr << "  TAXA block found" << endl;
        //      if( logf_open )
        //          taxa->Report(logf);
        //  }

        if( inf_open )
            inf.close();
        inf_open = false;

        //     if( !trees->IsEmpty() ) {
        //         cerr << "  TREES block found" << endl;
        //         if( logf_open )
        //            trees->Report(logf);




        //     }



        //    if( !assumptions->IsEmpty() ) {
        //        cerr << "  ASSUMPTIONS block found" << endl;
        //        if( logf_open )
        //             assumptions->Report(logf);
        //    }

        //    if( !characters->IsEmpty() ) {
        //        cerr << "  CHARACTERS block found" << endl;
        //        if( logf_open )
        //            characters->Report(logf);
        //    }
    }
    else
    {
        cerr << endl;
        cerr << "Oops! Could not find specified file: " << fn << endl;
    }
}



/**
* @method HandleExecute [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Handles everything after the EXecute keyword and the terminating
 * semicolon.  Flushes all blocks before executing file specified,
 * and no warning is given of this
 */
void BROWNIE::HandleExecute( NexusToken& token )
{
    // Issuing the EXECUTE command from within a file is a no-no
    //
    if( inf_open ) {
        errormsg = "Cannot issue execute command from within a BROWNIE block";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    // Get the file name to execute
    //
    token.GetNextToken();
    if( token.Equals("?") ) {
        message="Usage: Exe <file-name>\n\n";
        message+="Executes nexus-formatted input file.\nThis file may contain continuous character data, trees, and Brownie blocks.\nThe file name, if put in single quotes, may contain the path to the file.\n\n";
        PrintMessage();
        token.GetNextToken(); //to remove the remaining token
    }
    else {
        nxsstring fn = token.GetToken();

        // get the semicolon terminating the EXECUTE command
        //
        token.GetNextToken();

        if( !token.Equals(";") ) {
            errormsg = "Expecting ';' to terminate the EXECUTE command, but found ";
            errormsg += token.GetToken();
            errormsg += " instead";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }

        // Before going through with this, make sure we're not going to overwrite
        // any stored blocks
        bool stuff_stored = !taxa->IsEmpty();
        stuff_stored = ( stuff_stored || !trees->IsEmpty() );
        stuff_stored = ( stuff_stored || !assumptions->IsEmpty() );
        stuff_stored = ( stuff_stored || !characters->IsEmpty() );
        if( stuff_stored && UserSaysOk( "Ok to delete?", "Data has already been read and stored" ) )
            PurgeBlocks();
        else if( stuff_stored ) {
            message = "\nExecute command aborted.";
            PrintMessage();
            return;
        }

        if( FileExists( fn.c_str() ) )
        {
            cerr << endl;
            cerr << "Opening " << fn << "..." << endl;

            ifstream inf( fn.c_str(), ios::binary | ios::in );

            inf_open = true;
            ifstream intreefile;
            intreefile.open(fn.c_str(),ios::in);
            if (!intrees.ReadTrees(intreefile))
            {
                message="No trees read from file\n";
                PrintMessage();
            }
            intreefile.close();

            NexusToken ftoken(inf);

            try {
                Execute( ftoken );
            }
            catch( XNexus x )
            {
                NexusError( errormsg, x.pos, x.line, x.col );
                Reset();
            }



            //    if( !taxa->IsEmpty() ) {
            //        cerr << "  TAXA block found" << endl;
            //        if( logf_open )
            //            taxa->Report(logf);
            //    }

            if( inf_open )
                inf.close();
            inf_open = false;

            //    if( !trees->IsEmpty() ) {
            //        cerr << "  TREES block found" << endl;
            //        if( logf_open )
            //             trees->Report(logf);




            //   }



            // if( !assumptions->IsEmpty() ) {
            //   cerr << "  ASSUMPTIONS block found" << endl;
            // if( logf_open )
            //   assumptions->Report(logf);
            //    }

            //  if( !characters->IsEmpty() ) {
            //       cerr << "  CHARACTERS block found" << endl;
            //       if( logf_open )
            //           characters->Report(logf);
            //   }
        }
        else
        {
            cerr << endl;
            cerr << "Oops! Could not find specified file: " << fn << endl;
        }
    }
}

/**
* @method HandleHelp [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the HELP command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleHelp( NexusToken& token )
{
    // Retrieve all tokens for this command, stopping only in the event
    // of a semicolon or an unrecognized keyword
    //
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading HELP command.\n\nTry typing \"";
            errormsg +=token.GetToken();
            errormsg +=" ?\" instead [without the quotes] if you want help for the command.";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }

	message = "\nGeneral commands:";
    message += "\n  help           -> shows this message";
    message += "\n  exe            -> executes nexus file";
    message += "\n  log            -> log output";
    message += "\n  echo           -> copies your commands into a batch file";
    //message += "\n  gettree        -> loads tree file";
    message += "\n  blocks         -> reports on blocks currently stored";
    message += "\n  showtree       -> displays currently loaded tree(s)";
	message += "\n  choose         -> chooses tree or char for analysis";
    message += "\n  taxset         -> stores a taxset";
	message += "\n  citation       -> outputs list of relevant papers for your analyses";
	message += "\n  quit           -> terminates application";
	message += "\n\nNumerical optimization settings:";
	message += "\n  [ALPHA] set    -> sets options";
    message += "\n  [ALPHA] numopt -> sets parameters for numerical optimization functions";
	message += "\n\nCharacter evolution:";
    message += "\n  ratetest       -> does censored rate test";
    message += "\n  vcv            -> outputs a variance-covariance matrix";
    message += "\n  [ALPHA] tipvariance    -> allows program to deal with variance in taxon means";
    message += "\n  [ALPHA] model  -> sets model of character evolution (OU, BM, etc)";
    message += "\n  [ALPHA] opt    -> gets score for chosen taxset for chosen model";
    message += "\n  [ALPHA] export -> exports a tree and data in Pagel format";
	//message += "\n\nSpecies delineation and  tree search:";
	//message += "\n  [BETA] hs      -> perform a heuristic search";
	//message += "\n  [BETA]jackknife-> perform a jackknife search";
	//message += "\n  [BETA] compare -> compare triplet overlap for coalescent trees";
    message += "\n\nType \"commandname ?\" [without the quotes]\nfor help on any command.";
    PrintMessage();
}

/**
* @method HandleBlocks [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the Blocks command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleBlocks( NexusToken& token )
{
    bool finishexecuting=true;
    // Retrieve all tokens for this command, stopping only in the event
    // of a semicolon or an unrecognized keyword
    //
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if ( token.Equals("?") ) {
            message="Usage: Blocks\n\n";
            message+="Reports status of loaded nexus blocks\n";
            PrintMessage();
            finishexecuting=false;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Blocks command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
    if (finishexecuting) {
        message = "\nNexus blocks currently stored:";
        PrintMessage();
        if( !taxa->IsEmpty() ) {
            cerr << "\n  TAXA block found" << endl;
            taxa->Report(cerr);
            if( logf_open )
                taxa->Report(logf);
        }

        if( !trees->IsEmpty() ) {
            cerr << "\n  TREES block found" << endl;
            trees->Report(cerr);
            if( logf_open )
                trees->Report(logf);
        }

        if( !assumptions->IsEmpty() ) {
            cerr << "\n  ASSUMPTIONS block found" << endl;
            assumptions->Report(cerr);
            if( logf_open )
                assumptions->Report(logf);
        }

        if( !characters->IsEmpty() ) {
            cerr << "\n  CHARACTERS block found" << endl;
            characters->Report(cerr);
            if( logf_open )
                characters->Report(logf);
        }
    }
}


/**
* @method HandleDebug [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 */
void BROWNIE::HandleDebug( NexusToken& token )
{
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            debugmode=true;
            break;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Debug command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
    message="Now entering debug mode\n";
    PrintMessage();
}






/**
* @method HandleGettrees [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the GETTREES command needs to be parsed
 * from within the BROWNIE block.
 *
 * If there are trees already in the executed Nexus block, no need
 *   to run this.
 */
void BROWNIE::HandleGettrees( NexusToken& token )
{
    if( !trees->IsEmpty() ) {
        if( UserSaysOk( "Ok to delete?", "Trees have already been read and stored" ) ) {
            Detach( trees );
            delete trees;
            trees = new TreesBlock(*taxa);
            Add( trees );
        }
        else {
            message = "\nGetTrees command aborted.";
            PrintMessage();
            return;
        }
    }

    nxsstring fn;
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("File") ) {
            fn=GetFileName(token);
            break;
        }
        else {
            fn=token.GetToken();
        }
    }
    //input stream
    ifstream intreefile;
    intreefile.open(fn.c_str(),ios::in);




    if (!intrees.ReadTrees(intreefile))
    {
        errormsg="Failed to read trees";
        throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }
    intreefile.close();
}

/**
* @method HandleShowtree [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the SHOWTREE command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleShowtree( NexusToken& token )
{
    nxsstring numbernexus;
    bool finishexecuting=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (finishexecuting) {
                intrees.ShowTrees(cout);
            }
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="Usage: ShowTree\nDisplays all loaded trees. No option yet for just displaying the chosen tree.";
            PrintMessage();
            finishexecuting=false;
        }
    }
}


/**
* @method HandleHeuristicSearch [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 */
void BROWNIE::HandleHeuristicSearch( NexusToken& token )
{
    nxsstring numbernexus;
    bool finishexecuting=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (finishexecuting) {
				if(exhaustive) {
					DoExhaustiveSearch();
				}
				else if (!jackknifesearch) {
					DoHeuristicSearch();
				}
				else {
					ofstream jacktreef;
					nxsstring jacktreename=treefilename;
					jacktreename+=".jack.tre";
					jacktreef.open(jacktreename.c_str());
					jacktreef<<"#nexus\nbegin trees;\n"; //change these couts to jacktreef
					jacktreef.close();
					for (jackrep=1;jackrep<=jreps;jackrep++) { //a jackknife search consists of several heuristic searches, with the weights of several input trees set to zero
																//Each gene tree's weight is used with the jackknife proportion to get a deletion vector; using the weight ensures that the genes have the same expectation of weighted representation in the final set (think of combining  bootstrap samples for one gene with regular samples for another).
//NOTE: Need to deal with tree weights properly in GetGTPScoreNew; have not added anything to GetTripletScore; need to deal with how tree weights are used to calculate when trees come from different genes [perhaps add a new weight vector?]

//Use this in place of the original for ..intrees to only deal with selected trees
//	for (int i=0;i<intrees.GetNumTrees();i++) {
//		bool usethistree=true;
//		if (jackknifesearch) {
//			if (jackknifevector[i]==0) {
//				usethistree=false;
//			}
//		}
//		if (usethistree) {
//          STUFF originally just under the intrees for loop
//      }
//  }
						DoHeuristicSearch();	
						jacktreef.open(jacktreename.c_str(), ios::out | ios::app);
						jacktreef<<jackknifetreestooutput;
						jacktreef.close();

					}
					jacktreef.open(jacktreename.c_str(), ios::out | ios::app);
					jackknifesearch=false;
					jacktreef<<"end;";
					jacktreef.close();
				}
            }
            break;
        }
        else if( token.Abbreviation("?") ) {
			if (!jackknifesearch) {
				message="Usage: HSearch [options]\nDoes a heuristic search ";
			}
			else {
				message="Usage: Jackknife [options]\nDoes a jackknife search, deleting input trees ";
			}
            message+="Available options:\n\n";
            message+="Keyword ---- Option type --------------------------- Current setting --";
			if (jackknifesearch) {
				message+="\nJReps        <integer-value>                         ";
				message+=jreps;
				message+="\nPctDelete    <real-value>                            ";
				message+=pctdelete;
			}
            message+="\nNReps        <integer-value>                         ";
            message+=nreps;
            //message+="\nTimeLimit    <integer-value>|None                    None";
            //message+="\nClock        Wall|CPU                                Wall";
            message+="\nRearrLimit   <integer-value>|None                    ";
			if (rearrlimit<0) {
					message+="None";
			}
			else {
				message+=rearrlimit;
			}
            message+="\nAssignFixed  No|Yes                                  Yes";
            message+="\nSppNumFixed  No|Yes                                  Yes";
            message+="\nMaxNumSpp    <integer-value>                         ";
            message+=maxnumspecies;
            message+="\nMinNumSpp    <integer-value>                         ";
            message+=minnumspecies;
			message+="\nMinSamp      <integer-value>                         ";
            message+=minsamplesperspecies;
            message+="\nStructWt     <double>                                ";
            message+=structwt;
            message+="\nPThreshold   <double>                                ";
            message+=pthreshold;
			message+="\nSubsample    <double>                                ";
			message+=chosensubsampling;
            message+="\nMoveFreq     (number number number number number)    (";
            for (int i=0;i<5;i++) {
                message+=movefreqvector[i];
                message+=" ";
            }
            message+=")";
            message+="\nFile         <file-name>                             ";
            message+=treefilename;
            //message+="\nRecordSearch No|Yes                                  No";
            //message+="\nSearchFile   <file-name>                             ";
            //message+="\nStatus       No|Yes                                  Yes";
           // message+="\nSteepest     No|Yes                                  No";
			//message+="\nExhaustive   No|Yes                                  No";
            message+="\nShowTries    No|Yes                                  ";
            message+="\n\nNReps: Number of random starting species trees to use";
            //message+="\nTimeLimit: Limit search to X seconds";
            //message+="\nClock: Count seconds for time limit using actual elapsed time ('Wall'->Clock on a wall) or CPU time"; //NOte to self: see discussion online of time() fn and clock() fn in C++
            message+="\nRearrLimit: Limit search to X rearrangements for each nrep";
            message+="\nAssignFixed: Assignment of gene samples to species is not optimized during a search";
            message+="\nSppNumFixed: The total number of species is not optimized during the search (if set to No, then AssignFixed is also set to No)";
            message+="\nMaxNumSpp: The maximum number of species to split the samples into (only relevant if SppNumFixed==No)";
			message+="\nMinSamp: The minimum  number of samples per species";
            message+="\nMoveFreq: Sets the relative proportion of times to try\n\t1) Species tree branch swaps\n\t2) Moving samples from one species to another\n\t3) Increasing the number of species\n\t4) Decreasing the number of species\n\t5) Attempt to reroot the species tree\n  If these don't sum to one, the program will automatically correct this.";
            //message+="\nRecordSearch: Output each step to a file: assignments go to a .txt file, and trees go to a .tre file";
            //message+="\nSearchFile: If RecordSearch==Yes, the prefix to use for the output files.";
           // message+="\nStatus: Output status of search to screen (and a log file, if open).";
           // message+="\nSteepest: Whether to look at all rearrangements and then take the best one or just take the first better one.";
			message+="\nSubsample: How extensively to try taxon reassignments on leaf splits. \n\tA value of 1 means try all of the possible reassignments, \n\ta value of 2 means try the square root of all the possible assignments,\n\t3 means the cube root, etc. A higher number means a faster but less effective search.\n\tThe program won't let you try fewer than 10 assignments on average.";
            PrintMessage();
            finishexecuting=false;
        }
        else if(token.Abbreviation("NReps")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            nreps=atoi( numbernexus.c_str() ); //convert to int
        }
        else if(token.Abbreviation("MAxnumspp")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            maxnumspecies=atoi( numbernexus.c_str() ); //convert to int
        }
		else if(token.Abbreviation("MINSamp")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            minsamplesperspecies=GSL_MIN(atoi( numbernexus.c_str() ),taxa->GetNumTaxonLabels()); //convert to int
			message="All species must now have at least ";
			message+=minsamplesperspecies;
			message+=" samples in them.";
			PrintMessage();
        }
		else if(token.Abbreviation("REArrlimit")) {
            nxsstring numbernexus;
			message="Rearrangement limit set to ";
            numbernexus = GetFileName(token);
			if (numbernexus[0]=='n' || numbernexus[0]=='N') {
					rearrlimit=-1;
				message+="None";
			}
			else {
				rearrlimit=atoi(numbernexus.c_str() );
				message+=rearrlimit;
			}
			PrintMessage();
        }		
        else if( token.Abbreviation("STEepest") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n' || yesnoreplace[0] == 'N') {
                steepest=false;
            }
            else {
                steepest=true;
                cout<<"%%%%%%%%%%%%%% Warning: steepest hasn't been tested everywhere %%%%%%%%%%%%%%%%%%%"<<endl;
            }
        }
		else if( token.Abbreviation("EXhaustive") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n' || yesnoreplace[0] == 'N') {
                exhaustive=false;
            }
            else {
                exhaustive=true;
                message="Wow, an exhaustive search.";
				PrintMessage();
            }
        }
        else if( token.Abbreviation("SHOwtries") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n' || yesnoreplace[0] == 'N') {
                showtries=false;
            }
            else {
                showtries=true;
            }
        }
        else if( token.Abbreviation("MoveFreq") ) {
            token.GetNextToken();
            token.GetNextToken(); //eat the equals sign
            vector<double> temporarymovefreqvector;
            if (!token.Equals("(")) {
                errormsg="Expecting next token to be a \'(\' but instead got ";
                errormsg+=token.GetToken();
                throw XNexus( errormsg);
            }
            int inputcount=0;
            while (!token.Equals(")")) {
                nxsstring numbernexus;
                numbernexus=GetNumberOnly(token);
                if (numbernexus!=")") {
                    temporarymovefreqvector.push_back(atof( numbernexus.c_str() ));
                    inputcount++;
                }
                else {
                    break;
                }
            }
            if (inputcount!=5) {
                errormsg="You should have entered five frequencies, you entered ";
                errormsg+=inputcount;
                throw XNexus( errormsg);
            }
            else {
                double sumoffreqs=temporarymovefreqvector[0]+temporarymovefreqvector[1]+temporarymovefreqvector[2]+temporarymovefreqvector[3]+temporarymovefreqvector[4];
                movefreqvector.clear();
                for (int i=0; i<temporarymovefreqvector.size(); i++) {
                    movefreqvector[i]=(temporarymovefreqvector[i])/sumoffreqs;
                }

            }
        }

        else if(token.Abbreviation("MInumspp")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            minnumspecies=atoi( numbernexus.c_str() ); //convert to int
        }
        else if(token.Abbreviation("FIlename")) {
            treefilename=GetFileName(token);
        }
        else if(token.Abbreviation("STRuctwt")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            double proposedwt=atof( numbernexus.c_str() ); //convert to double
            if ((proposedwt<0) || (proposedwt>1)) {
                errormsg="Error: The weight for population substructure must be between 0 and 1 (inclusive)";
                throw XNexus( errormsg);
            }
            else {
                structwt=proposedwt;
                message="You have chosen to have the gene duplication cost to make up ";
                double gtpwt=(1-structwt)*100;
                message+=gtpwt;
                message+="% and the substructure cost to make up the remaining ";
                message+=100-gtpwt;
                message+="%.";
                PrintMessage();
            }
        }
        else if(token.Abbreviation("PThreshold")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            double proposedthresh=atof( numbernexus.c_str() ); //convert to double
            if ((proposedthresh<0) || (proposedthresh>1)) {
                errormsg="Error: The p-value threshold must be between 0 and 1 (inclusive)";
                throw XNexus( errormsg);
            }
            else {
                pthreshold=proposedthresh;
                message="You have chosen to have structure p values > ";
                message+=pthreshold;
                message+=" ignored";
                PrintMessage();
            }
        }
		else if(token.Abbreviation("SUBsample")) {
			nxsstring numbernexus;
            numbernexus = GetNumber(token);
            double proposedsubsample=atof( numbernexus.c_str() ); //convert to double
			if (proposedsubsample<1) {
				errormsg="Error: The subsample threshold must be greater than one";
                throw XNexus( errormsg);
			}
			else {
				chosensubsampling=proposedsubsample;
				message="You have chosen a split leaf subsample value of ";
				message+=chosensubsampling;
				PrintMessage();
			}
			
			
		}
    }
}

vector<double> BROWNIE::GetCombinedScore(ContainingTree *SpeciesTreePtr)
{
    //We do this so we don't bother doing the GTP or triplet calculations if they're not needed.
	triplettoohigh=false;
	gtptoohigh=false;
    vector<double> scorevector;
	bool calculatescore=true;
    double combinedscore=0;
    double tripletscore=0;
    double gtpscore=0;
	if (minsamplesperspecies>1) {
		int maxspecies=0;
		for (int i=0;i<convertsamplestospecies.size();i++) {
			maxspecies=GSL_MAX(maxspecies,convertsamplestospecies[i]);
		}
		vector<int> speciesvector(maxspecies,0);
		for (int i=0;i<convertsamplestospecies.size();i++) {
			speciesvector[(convertsamplestospecies[i])-1]++;
		}
		for (int i=0;i<speciesvector.size();i++) {
			if (speciesvector[i]<minsamplesperspecies) {
				calculatescore=false;
				combinedscore=GSL_POSINF;
				tripletscore=GSL_POSINF;
				gtpscore=GSL_POSINF;
			}
		}		
	}
	if (calculatescore) {
		if (structwt>0) {
			tripletscore=double(GetTripletScore(SpeciesTreePtr));
        //	(*SpeciesTreePtr).Write(cout);
        //	cout<<endl;
        //	for (int i=0;i<convertsamplestospecies.size();i++) {
        //			cout<<convertsamplestospecies[i]<<" ";
        //	}
        //	cout<<endl;
        //	cout<<"\ntriplet score = "<<newscore<<endl<<endl;
        //  combinedscore+=newscore;
			
		}		
		if (structwt<1 && (triplettoohigh==false)) {
        //combinedscore+=(1.0-structwt)*GetGTPScore(SpeciesTreePtr);
        // cout<<"Mike gave "<<combinedscore/(1.0-structwt)<<"\nI gave "<<GetGTPScoreNew(SpeciesTreePtr)<<endl;
        //cout<<"SpeciesTree"<<endl;
        //(*SpeciesTreePtr).Write(cout);
        //cout<<endl;
			gtpscore=double(GetGTPScoreNew(SpeciesTreePtr));
		}
		combinedscore=((1.0-structwt)*gtpscore)+(structwt*tripletscore);
	}
    scorevector.push_back(combinedscore);
    scorevector.push_back(gtpscore);
    scorevector.push_back(tripletscore);
    return scorevector;
}

//checks to make sure that the convertsamplestospecies vector is not missing assignments to any species (i.e, doesn't consist only of species 1, 3 and 4). If fix==true, it'll fix this. It will return true if check
bool BROWNIE::CheckConvertSamplesToSpeciesVector(bool fix)
{
    int maxspecies=0;
    bool goodshape=true;
    int firstbadspecies=0;
    for (int i=0;i<convertsamplestospecies.size();i++) {
        maxspecies=GSL_MAX(maxspecies,convertsamplestospecies[i]);
    }
    vector<int> speciesvector(maxspecies,0);
    for (int i=0;i<convertsamplestospecies.size();i++) {
        speciesvector[(convertsamplestospecies[i])-1]++;
    }
    for (int i=0;i<speciesvector.size();i++) {
        if (speciesvector[i]==0) {
            goodshape=false;
            firstbadspecies=i+1;
            break;
        }
    }
    if (!goodshape && fix) {
        for (int i=0;i<convertsamplestospecies.size();i++) {
            if (convertsamplestospecies[i]>firstbadspecies) {
                convertsamplestospecies[i]=convertsamplestospecies[i]-1;
            }
        }
        goodshape=CheckConvertSamplesToSpeciesVector(true);

    }
    return goodshape;
}

bool BROWNIE::CombineSpeciesWithTooFewSamples(bool fix) { //The idea of this is to fix the assignment vector if there are too few samples per species by combining the too small species with the next species and then fixing the resulting hole
	bool goodsamplenumber=true;
	if (minsamplesperspecies>1) {
		int maxspecies=0;
		for (int i=0;i<convertsamplestospecies.size();i++) {
			maxspecies=GSL_MAX(maxspecies,convertsamplestospecies[i]);
		}
		vector<int> speciesvector(maxspecies,0);
		for (int i=0;i<convertsamplestospecies.size();i++) {
			speciesvector[(convertsamplestospecies[i])-1]++;
		}
		for (int i=0;i<speciesvector.size();i++) {
			if (speciesvector[i]<minsamplesperspecies) {
				goodsamplenumber=false;
				if (fix) {
					for (int j=0;j<convertsamplestospecies.size();j++) {
						if (convertsamplestospecies[j]==i+1) {
							if ((i+1)<maxspecies) { //so there's a species to the right to combine with.
								convertsamplestospecies[j]=convertsamplestospecies[j]+1;
							}
							else { 
								convertsamplestospecies[j]=convertsamplestospecies[j]-1;
							}
						}
					}
					bool goodshape=CheckConvertSamplesToSpeciesVector(true);
					goodsamplenumber=CombineSpeciesWithTooFewSamples(true);
				}
			}
		}		
	}
	return goodsamplenumber;
}

bool BROWNIE::CheckConvertSamplesToSpeciesVectorSpNum(int actualmaxspeciesnum)
{
    int maxspecies=0;
    bool goodshape=true;
    int firstbadspecies=0;
    for (int i=0;i<convertsamplestospecies.size();i++) {
        maxspecies=GSL_MAX(maxspecies,convertsamplestospecies[i]);
    }
	if (maxspecies!=actualmaxspeciesnum) {
		goodshape=false;
	}
	if (goodshape) {
		vector<int> speciesvector(maxspecies,0);
		for (int i=0;i<convertsamplestospecies.size();i++) {
			speciesvector[(convertsamplestospecies[i])-1]++;
		}
		for (int i=0;i<speciesvector.size();i++) {
			if (speciesvector[i]==0) {
				goodshape=false;
				firstbadspecies=i+1;
				break;
			}
		}
	}
	return goodshape;
}

bool BROWNIE::MoveSamples(vector<int> Originalconvertsamplestospecies)
{
    bool morereassignments=true;
    int sampletomove=SamplesToMove.back();
    int destination=SampleDestinations.back();
    for(int i=0;i<(CladeVector[sampletomove]).size();i++) {
        convertsamplestospecies[(CladeVector[sampletomove][i])]=destination;
    }
    SampleDestinations.pop_back();
    if (SampleDestinations.size()==0) { //we've finished trying to assign all destinations
        convertsamplestospecies.swap(Originalconvertsamplestospecies); //need to use the original vector; if this try results in a better vector, we'll throw out the Sample___ vectors anyway; if not, we'll go back to the Originalconvertsamplestospecies vector, so we have to make sure any moves are allowed under that vector
        int nspecies=0;
        for (int i=0;i<convertsamplestospecies.size();i++) {
            nspecies=GSL_MAX(convertsamplestospecies[i],nspecies);
        }
        vector<int> TempSampleDestinations;
        if (SamplesToMove.size()<2) { //we were already on the last sample, too
            morereassignments=false;
        }
        else {
            bool enoughdestinations=false;
            while (!enoughdestinations && SamplesToMove.size()>1) {
                SamplesToMove.pop_back(); //we're done trying to move that sample
                SampleDestinations.clear();
                for (int i=0;i<nspecies;i++) {
                    if(TestMoveSamples(SamplesToMove.back(),i+1)) {
                        TempSampleDestinations.push_back(i+1);
                        enoughdestinations=true;
                    }
                }
            }
            morereassignments=enoughdestinations;
            if(enoughdestinations) {
                gsl_permutation * v = gsl_permutation_alloc (TempSampleDestinations.size());
                gsl_permutation_init (v);
                gsl_ran_shuffle (r, v->data, TempSampleDestinations.size(), sizeof(size_t));
                for (int i=0; i<TempSampleDestinations.size(); i++) {
                    SampleDestinations.push_back(TempSampleDestinations[gsl_permutation_get (v,i)]);
                }
                gsl_permutation_free(v);
            }
        }
        convertsamplestospecies.swap(Originalconvertsamplestospecies);
    }
    return morereassignments;
}

bool BROWNIE::TestMoveSamples(int Sample, int Destination)
{
    bool validmove=true;
    int maxspeciesorig=0;
    for (int i=0; i<convertsamplestospecies.size();i++) {
        maxspeciesorig=GSL_MAX(maxspeciesorig,convertsamplestospecies[i]);
    }
    vector<int> origconvertsamplestospecies=convertsamplestospecies;
    for(int i=1;i<=(CladeVector[Sample]).size();i++) {
		//convertsamplestospecies[(CladeVector.at(Sample)).at(i)]=Destination;
		int convertsamplestospeciessize=convertsamplestospecies.size();
		int cladevectorsamplesize=(CladeVector[Sample]).size();
		int cladevectorsize=CladeVector.size();
		if(cladevectorsamplesize>convertsamplestospeciessize || cladevectorsamplesize<0) {
			nxsstring outputname=treefilename;
			outputname+=".CladeVector.txt";
			ofstream cladevector;
			cladevector.open( outputname.c_str(), ios::out | ios::app );
			for (int k=0;k<CladeVector.size();k++) {
				cladevector<<endl;
				for (int j=0;j<(CladeVector[k]).size();j++) {
					cladevector<<" "<<CladeVector[k][j];
				}
			}
			cladevector<<endl;
			cladevector<<endl<<"convertsamplestospeciessize = "<<convertsamplestospeciessize<<endl;
			cladevector<<"cladevectorsamplesize = "<<cladevectorsamplesize<<endl;
			cladevector<<"cladevectorsize = "<<cladevectorsize<<endl;
			cladevector<<"Sample = "<<Sample<<endl<<"Destination = "<<Destination<<endl<<"i = "<<i<<endl;
			cladevector.close();
		}
        convertsamplestospecies[(CladeVector[Sample][i-1])]=Destination; // so it doesn't run if CladeVector[Sample].size()=0
    }
    validmove=CheckConvertSamplesToSpeciesVector(false);
    if(validmove) {
        if(convertsamplestospecies==origconvertsamplestospecies) {
            validmove=false; //since this induces no change
        }
    }
    if(validmove) {
        int maxspeciesfinal=0;
        for (int i=0; i<convertsamplestospecies.size();i++) {
            maxspeciesfinal=GSL_MAX(maxspeciesfinal,convertsamplestospecies[i]);
        }
        if (maxspeciesfinal!=maxspeciesorig) {
            validmove=false; //need to have same number of species at start and end
        }
    }
    convertsamplestospecies.swap(origconvertsamplestospecies);
    return validmove;
}

void BROWNIE::DoExhaustiveSearch()
{
	message="This will do an exhaustive search for up to 6 species. This will take some time.";
	PrintMessage();
	bestscore=GSL_POSINF;
	vector<double> nextscorevector;
	double nextscore;
	message="Now starting with 1 species";
	PrintMessage();
	ContainingTree CurrentTree;
	CurrentTree.RandomTree(1);
	CurrentTree.ConvertTaxonNamesToRandomTaxonNumbers();
	convertsamplestospecies.clear();
	for (int i=0;i<taxa->GetNumTaxonLabels();i++) {
		convertsamplestospecies.push_back(1);
		cout<<" 1";
	}
	nextscorevector=GetCombinedScore(&CurrentTree);
	nextscore=nextscorevector[0]; 
	bestscore=nextscore; //this must be the best score, as it's the first tree
	//TotalScores.push_back(nextscorevector[0]);
	//GTPScores.push_back(nextscorevector[1]);
	//StructScores.push_back(nextscorevector[2]);
	FormatAndStoreBestTree(&CurrentTree,nextscorevector);
	
	message="Score: ";
	message+=nextscore;
	PrintMessage();
	message="Now starting with 2 species";
	ContainingTree TreeTwo;
	int maxspecies=2;
	TreeTwo.RandomTree(maxspecies);
	TreeTwo.ConvertTaxonNamesToRandomTaxonNumbers();
	bestscore=DoAllAssignments(bestscore,maxspecies,&TreeTwo);
	message="Now starting with 3 species";
	ContainingTree TreeThree;
	maxspecies=3;
	TreeThree.RandomTree(maxspecies);
	TreeThree.ConvertTaxonNamesToRandomTaxonNumbers();
	bestscore=DoAllAssignments(bestscore,maxspecies,&TreeThree);
	cout<<endl<<"Best trees overall"<<endl<<endl;
	for (int i=0; i<FormattedBestTrees.size(); i++) {
		(FormattedBestTrees[i]).Update();
		(FormattedBestTrees[i]).GetNodeDepths();
		(FormattedBestTrees[i]).Draw(cout);
		cout<<endl;
	}
	message="\n\n#nexus\nbegin trees;\n";
	for (int i=0; i<FormattedBestTrees.size(); i++) {
		message+="tree sptree";
		message+=i+1;
		message+=" = ";
		if (unrooted==1) {
			message+="[&U] ";
		}
		else {
			message+="[&R] ";
		}
		message+=ReturnFinalSpeciesTree(FormattedBestTrees[i]);
		message+="\n";
	}
	message+="end;\n\n";
	PrintMessage();
	message="Best score = ";
	message+=bestscore;
	PrintMessage();
	
}

double BROWNIE::DoAllAssignments(double bestscore, int maxspecies, ContainingTree *SpeciesTree )
{
	vector<double> nextscorevector;
	double nextscore;
	vector <int> countpertaxon(taxa->GetNumTaxonLabels(),1);
	for (int i=0;i<taxa->GetNumTaxonLabels();i++) {
		convertsamplestospecies[i]=1;
		cout<<" "<<convertsamplestospecies[i];
	}
	for (int assignmentrep=0;assignmentrep<pow(maxspecies,taxa->GetNumTaxonLabels());assignmentrep++) {
		if (CheckConvertSamplesToSpeciesVectorSpNum(maxspecies)) { //we have the right number of species
			nextscorevector=GetCombinedScore(SpeciesTree);
			nextscore=nextscorevector[0]; 
			cout<<"  Score: "<<nextscore<<" Best: "<<GSL_MIN(nextscore,bestscore);
			if (nextscore==bestscore) {
			//	TotalScores.push_back(nextscorevector[0]);
			//	GTPScores.push_back(nextscorevector[1]);
			//	StructScores.push_back(nextscorevector[2]);				
				FormatAndStoreBestTree(SpeciesTree,nextscorevector);
			}
			else if (nextscore<bestscore) {
				bestscore=nextscore;
				FormattedBestTrees.clear();
				TotalScores.clear();
				GTPScores.clear();
				StructScores.clear();
				BestConversions.clear();
				//TotalScores.push_back(nextscorevector[0]);
				//GTPScores.push_back(nextscorevector[1]);
				//StructScores.push_back(nextscorevector[2]);
				
                FormatAndStoreBestTree(SpeciesTree,nextscorevector);
			}
		}
		else {
				cout<<"  Invalid assignment for "<<maxspecies<<" species";
		}
		cout<<endl;
	//	for (int i=0;i<countpertaxon.size();i++) {
			//cout<<" "<<countpertaxon[i];
	//	}
		//cout<<endl;
		for (int i=0;i<countpertaxon.size();i++) {
			countpertaxon[i]++;
		//	cout<<" maxspecies = "<<maxspecies<<" i = "<<i<<" pow = "<<pow(maxspecies,i)<<" countpertaxon = "<<countpertaxon[i];
			if (countpertaxon[i]>pow(maxspecies,i)) {
			//	cout<<" must change ";
				convertsamplestospecies[i]++;
				countpertaxon[i]=1;
				if (convertsamplestospecies[i]>maxspecies) {
					convertsamplestospecies[i]=1;
				}
			}
			cout<<" "<<convertsamplestospecies[i];
		}

	}	
	return bestscore;
}


void BROWNIE::DoHeuristicSearch()
{
	if (jackknifesearch) {
			message="\n---------- Now starting jackknife search replicate ";
		message+=jackrep;
		message+=" ----------";
		PrintMessage();
		int includedtrees=0;
		while (includedtrees<2 || includedtrees==trees->GetNumTrees()) {
			jackknifevector.clear();
			geneidvector.clear();
			double currenttotalwt=0;
			int genenumber=1;
			for (int curtreenum=0; curtreenum<trees->GetNumTrees(); curtreenum++) {
				double newweight=trees->GetTreeWeight(curtreenum);
				currenttotalwt+=newweight;
				if (currenttotalwt>1) {
					currenttotalwt=newweight;
					genenumber++;
				}
				geneidvector.push_back(genenumber);
				jackknifevector.push_back(gsl_ran_bernoulli(r, newweight*(1-pctdelete))); //we randomly select trees
				if (jackknifevector.back()==1) {
					includedtrees++;
				}
			}
		}
		message="----------   ";
		message+=includedtrees;
		message+=" out of ";
		message+=trees->GetNumTrees();
		message+=" input trees are included.  ----------\n\n";
		PrintMessage();
	}
    message="Creating initial neighbor-joining tree for samples, based on triplet overlap. Please be patient.\n\nNow getting distances...";
    PrintMessage();

    GetTaxonTaxonTripletDistances();
    message="Computing tree...";
    PrintMessage();
    ContainingTree NJBrlenTree=ComputeTripletNJTree();
	ContainingTree TripletSupportBrlenTree=NJBrlenTree;
	TripletSupportBrlenTree.Update();
	//message="Tree of samples with triplet support brlen: ";
	//PrintMessage();
	//TripletSupportBrlenTree.Draw(cout);
	NodeIterator <Node> l (TripletSupportBrlenTree.GetRoot()); //goes from tips down
	NodePtr currentnode=l.begin();
	int nodecount=0;
	while (currentnode) {
		if (currentnode->IsLeaf()) {
			int leafnumtrans=(currentnode->GetLeafNumber())-1;
			currentnode->SetEdgeLength(CladeVectorTripletSupport[leafnumtrans]);
		}
		else if(currentnode!=TripletSupportBrlenTree.GetRoot()) {
			int nodenumber=nodecount+TripletSupportBrlenTree.GetNumLeaves();
			currentnode->SetEdgeLength(CladeVectorTripletSupport[nodenumber]);
			nodecount++;
		}
		currentnode=l.next();
	}
	if( logf_open ) {
		TripletSupportBrlenTree.Draw(logf);
		logf << endl;
		logf << "#nexus\nbegin trees;\ntree GuideTreeTripletSupportBrlen = ";
		TripletSupportBrlenTree.Write(logf);
		logf << "\ntree GuideTreeTripletSupportBrlen = ";
		NJBrlenTree.Write(logf);
		logf<<endl<<"end;"<<endl;
	}
    //This is just a rudimentary search: later, add options like time limits, changing the number of species, etc.
    bestscore=GSL_POSINF;
    //ADD STORAGE OF STARTING VALUE
    vector<int> intialconvertsamplestospeciesvector=convertsamplestospecies;
    double nextscore;
    double nextgtpscore;
    double nexttripletscore;
    vector<double> nextscorevector;
    nxsstring scoretype;
	message="Now starting the search proper.\nA \">\" before a score indicates that calculation of that score was aborted once the score for that move exceeded the best local score\n";
	if (!jackknifesearch) {
		message+="\nRep\tMoves\t#Spp\tType\tQual\tCombScore\t      GTP\t   Struct\t    Local\t   Global\tNTrees\tRemaining";
	}
	else {
	    message+="\nJackRep\tRep\tMoves\t#Spp\tType\tQual\tCombScore\t      GTP\t   Struct\t    Local\t   Global\tNTrees\tRemaining";
	}
    PrintMessage();

    for (int replicate=1;replicate<=nreps;replicate++) {
        convertsamplestospecies=intialconvertsamplestospeciesvector;
        ContainingTree StartingTree;
        //cout<<"Starting vector = "<<endl;
        //convertsamplestospecies=intialconvertsamplestospeciesvector;
        //for (int i=0;i<convertsamplestospecies.size();i++) {
        //		cout<<convertsamplestospecies[i]<<" ";
        //}
        //cout<<endl;
        bool assignmentbasedontriplettree=true;
        if (convertsamplestospecies.size()==0) {
            if (assignmentbasedontriplettree) { //use the triplet tree to get assignments
				if (gsl_ran_flat(r,0,1)<tripletdistthreshold) {
				//New method, based on splitting on longest internal branches in starting nj tree. Basically, split on internal branches with longer than average lengths
					int nsamples=taxa->GetNumTaxonLabels();
					int numnontrivialclades=nsamples-2;
					convertsamplestospecies.assign(nsamples,1);
					if (debugmode) {
						cout<<"assembling initial convertsamplestospeciesvector, using NJ tree distances"<<endl;
						for (int k=0;k<convertsamplestospecies.size();k++) {
							cout<<convertsamplestospecies[k]<<" ";
						}
						cout<<endl;
					}
					int currentspecies=2;
					for (int currentpos=CladeVectorNJBrlen.size()-1;currentpos>=nsamples-1;currentpos--) { //start at the root, work up (based on other order on way down) 
						if (CladeVectorNJBrlen[currentpos]>meaninternalbrlen && (CladeVector[currentpos]).size()>=minsamplesperspecies) {
							for(int j=0;j<(CladeVector[currentpos]).size();j++) {
								convertsamplestospecies[(CladeVector[currentpos][j])]=currentspecies;
							}
							currentspecies++;
							if (debugmode) {
								for (int taxon=0; taxon<nsamples; taxon++) {
									cout<<convertsamplestospecies[taxon]<<" "<<taxa->GetTaxonLabel(taxon)<<endl;
								}
								cout<<endl;
							}
							
						}
					}
					bool goodshape=CheckConvertSamplesToSpeciesVector(true);
					goodshape=CombineSpeciesWithTooFewSamples(true);
					if (debugmode) {
						for (int taxon=0; taxon<nsamples; taxon++) {
							cout<<convertsamplestospecies[taxon]<<" "<<taxa->GetTaxonLabel(taxon)<<endl;
						}
						cout<<endl;
					}
					
					/* //Old method for getting starting assignments: tended to split good clades too often, not pay attention to relative support
						int nsamples=taxa->GetNumTaxonLabels();
					int numnontrivialclades=nsamples-2;
					convertsamplestospecies.assign(nsamples,1);
					int currentpos=nsamples-1;
					int currentspecies=2;
					while (currentpos<CladeVector.size()) {
						currentpos+=1+gsl_ran_binomial(r,.5,4);
						if (currentpos<CladeVector.size()) {
							for(int j=0;j<(CladeVector[currentpos]).size();j++) {
								convertsamplestospecies[(CladeVector[currentpos][j])]=currentspecies;
							}
						}
						currentspecies++;
					}
					bool goodshape=CheckConvertSamplesToSpeciesVector(true);
					goodshape=CombineSpeciesWithTooFewSamples(true);
					*/ //Old method for getting starting assignments
				}
				else { //Use triplet support distances
					int nsamples=taxa->GetNumTaxonLabels();
					int numnontrivialclades=nsamples-2;
					convertsamplestospecies.assign(nsamples,1);
					if (debugmode) {
						cout<<"assembling initial convertsamplestospeciesvector, using triplet support distances"<<endl;
						for (int k=0;k<convertsamplestospecies.size();k++) {
							cout<<convertsamplestospecies[k]<<" ";
						}
						cout<<endl;
					}
					int currentspecies=2;
					for (int currentpos=CladeVectorTripletSupport.size()-1;currentpos>=nsamples-1;currentpos--) { //start at the root, work up (based on other order on way down) 
						double supportvalue=CladeVectorTripletSupport[currentpos];
						if ((supportvalue>gsl_ran_flat(r,0.6,1)) && (CladeVector[currentpos]).size()>=minsamplesperspecies) { //only split on branches with 60% support or more
							for(int j=0;j<(CladeVector[currentpos]).size();j++) {
								convertsamplestospecies[(CladeVector[currentpos][j])]=currentspecies;
							}
							currentspecies++;
							if (debugmode) {
								for (int taxon=0; taxon<nsamples; taxon++) {
									cout<<convertsamplestospecies[taxon]<<" "<<taxa->GetTaxonLabel(taxon)<<endl;
								}
								cout<<endl;
							}
							
						}
					}
					bool goodshape=CheckConvertSamplesToSpeciesVector(true);
					goodshape=CombineSpeciesWithTooFewSamples(true);
					if (debugmode) {
						for (int taxon=0; taxon<nsamples; taxon++) {
							cout<<convertsamplestospecies[taxon]<<" "<<taxa->GetTaxonLabel(taxon)<<endl;
						}
						cout<<endl;
					}
				}
            }
            else {
                //message="You didn't do an intial assignment of taxa to species, so we'll do a random assignment";
                //PrintMessage();
                int ntax=taxa->GetNumTaxonLabels();
                //cout<<"ntax is "<<ntax<<endl;
                int samplesperspecies=GSL_MAX(minsamplesperspecies,1+gsl_ran_binomial(r,0.5,7)); // can change this; set now for on average 4.5 samples per species
                convertsamplestospecies.clear();
                vector <int> tempconvertsamplestospecies;
                int speciesid=1;
                int assignmentcount=0;
                for (int i=0;i<ntax;i++) {
                    tempconvertsamplestospecies.push_back(speciesid);
                    assignmentcount++;
                    if (assignmentcount==samplesperspecies) {
                        assignmentcount=0;
                        speciesid++;
                        samplesperspecies=GSL_MAX(minsamplesperspecies,1+gsl_ran_binomial(r,0.5,7)); //currently set for an average of 4.5 samples per species
                                                                       //cout<<"speciesid is "<<speciesid<<endl;
                    }
                }
				if (assignmentcount<minsamplesperspecies && assignmentcount>0) { //Means the last species has too few samples in it, so we'll merge it with a random earlier species
					for (int i=0;i<tempconvertsamplestospecies.size();i++) {
						if(tempconvertsamplestospecies[i]==speciesid) {
							tempconvertsamplestospecies[i]=1+gsl_ran_binomial(r,0.5,speciesid-2);
						}
					}
				}
                gsl_permutation * c = gsl_permutation_alloc (tempconvertsamplestospecies.size());
                gsl_permutation_init (c);
                gsl_ran_shuffle (r, c->data, tempconvertsamplestospecies.size(), sizeof(size_t));
                for (int i=0; i<tempconvertsamplestospecies.size(); i++) {
                    convertsamplestospecies.push_back(tempconvertsamplestospecies[gsl_permutation_get (c,i)]);
                    //cout<<tempconvertsamplestospecies[gsl_permutation_get (c,i)]<<endl;
                }
                gsl_permutation_free(c);
            }
        }
        int CurrentSppNum=0;
        for (int vectorpos=0;vectorpos<convertsamplestospecies.size();vectorpos++) {
            if (convertsamplestospecies[vectorpos]>CurrentSppNum) {
                CurrentSppNum=convertsamplestospecies[vectorpos];
            }
        }
        //cout<<"Starting vector = "<<endl;
        //for (int i=0;i<convertsamplestospecies.size();i++) {
        //	cout<<convertsamplestospecies[i]<<" ";
        //}
        //cout<<endl;
bestscorelocal=GSL_POSINF;
vector<ContainingTree> BestTreesThisRep;

//ContainingTree BestTree;
StartingTree.RandomTree(CurrentSppNum);
StartingTree.ConvertTaxonNamesToRandomTaxonNumbers();
//cout<<"Starting tree"<<endl;
//StartingTree.Write(cout);
//cout<<endl;
//StartingTree.Draw(cout);
//cout<<endl;
//cout<<"Root "<<StartingTree.GetRoot()<<endl;
//cout<<"Root child = "<<(StartingTree.GetRoot())->GetChild()<<endl;
//cout<<"Root child is leaf? "<<((StartingTree.GetRoot())->GetChild())->IsLeaf()<<endl;
//cout<<"\nTree health\n";
//StartingTree.ReportTreeHealth();
//cout<<"\n";
// char* inputforgtp=OutputForGTP(&CurrentTree);
// cout<<OutputForGTP(&StartingTree)<<endl;
// bestscorelocal=ReturnScore(OutputForGTP(&StartingTree));
//cout<<"currentsppnum = "<<CurrentSppNum<<endl;
//cout<<"starting tree health:\n";
//StartingTree.ReportTreeHealth();
vector<double> bestscorelocalvector=GetCombinedScore(&StartingTree);
bestscorelocal=bestscorelocalvector[0];
if (bestscorelocal==bestscore) {
    RawBestTrees.push_back(StartingTree);
	//TotalScores.push_back(bestscorelocalvector[0]);
	//GTPScores.push_back(bestscorelocalvector[1]);
	//StructScores.push_back(bestscorelocalvector[2]);	
    FormatAndStoreBestTree(&StartingTree,bestscorelocalvector);
    scoretype="*G\t";
}
else if (bestscorelocal<bestscore) {
    RawBestTrees.clear();
    RawBestTrees.push_back(StartingTree);
    FormattedBestTrees.clear();
	TotalScores.clear();
	GTPScores.clear();
	StructScores.clear();
	BestConversions.clear();
//	TotalScores.push_back(bestscorelocalvector[0]);
//	GTPScores.push_back(bestscorelocalvector[1]);
//	StructScores.push_back(bestscorelocalvector[2]);
    FormatAndStoreBestTree(&StartingTree,bestscorelocalvector);
    bestscore=bestscorelocal;
    scoretype="=G\t";
}
else {
    scoretype="*L\t";
}

//while (bestscorelocal<0) { //due to error in GTP
//    bestscorelocal=ReturnScore(OutputForGTP(&StartingTree));
//}
BestTreesThisRep.push_back(StartingTree);

//BestTree=StartingTree;
int movecount=0;
if (status) {
    //cout<<"\n\n"<<OutputForGTP(&CurrentTree)<<"\n\n";
    // cout<<"Starting tree: \n\n"; //Rewrite the draw function to allow output to a file
    // CurrentTree.Draw(cout);
    // cout<<"\nScore is "<<bestscore<<"\n";
}
// if (BestTrees.size()==0) {
//    BestTrees.push_back(CurrentTree);
// }
bool improvement=true;
while (improvement && (rearrlimit<0 || movecount<rearrlimit)) {
    //cout<<"\nimprovement, restarting\n";
    improvement=false;
    bool moreswaps=true;
    bool morereassignments=true;
    bool moreincreases=true;
    bool moredecreases=true;
    bool morererootings=true;
    // cout<<"moreswaps = "<<moreswaps<<" morereassignments = "<<morereassignments<<" moreincreases = "<<moreincreases<<" moredecreases = "<<moredecreases<<" morererootings = "<<morererootings<<endl;
    assert(BestTreesThisRep.size()>0);
    //for(int i=0;i<BestTreesThisRep.size();i++) {
    //    BestTreesThisRep[i].Write(cout);
    //    cout<<endl;
    // }
    // (BestTreesThisRep.back()).Draw(cout);
    // (BestTreesThisRep.back()).ReportTreeHealth();
    (BestTreesThisRep.back()).Update();
    //(BestTreesThisRep.back()).ReportTreeHealth();
    // cout<<"GetRoot: "<<(BestTreesThisRep.back()).GetRoot()<<endl;
    assert(BestTreesThisRep.size()>0);
    ContainingTree CurrentTree=BestTreesThisRep.back();
    CurrentTree.Update();
    CurrentTree.ResetBreakVector();
    CurrentTree.UpdateCherries();
    CurrentTree.SetLeafNumbers();
    //cout<<"GetNumLeaves = "<<CurrentTree.GetNumLeaves()<<endl;
    if ((sppnumfixed==true) || CurrentTree.GetNumLeaves()<=minnumspecies) {
        moredecreases=false;
    }
    if ((sppnumfixed==true) || CurrentTree.GetNumLeaves()>=maxnumspecies) {
        moreincreases=false;
    }
    if (movefreqvector[0]==0 || CurrentTree.GetNumLeaves()<3) {
        moreswaps=false;
    }
    if (movefreqvector[1]==0 || CurrentTree.GetNumLeaves()==1) {
        morereassignments=false;
    }
    if (movefreqvector[2]==0) {
        moreincreases=false;
    }
    if (movefreqvector[3]==0 || CurrentTree.GetNumLeaves()==1) {
        moredecreases=false;
    }
    if (movefreqvector[4]==0 || CurrentTree.GetNumLeaves()<3) {
        morererootings=false;
    }
    //Get list of cherries to collapse; we do this at the start so that we try each cherry at random but only once.
    vector<int> TempCherriesToMash;
    vector<int> CherriesToMash;
    for (int i=0;i<CurrentTree.GetNumCherries(); i++) {
        TempCherriesToMash.push_back(i);
    }
    if (CurrentTree.GetNumCherries()>0) {
        gsl_permutation * p = gsl_permutation_alloc (TempCherriesToMash.size());
        gsl_permutation_init (p);
        gsl_ran_shuffle (r, p->data, TempCherriesToMash.size(), sizeof(size_t));
        for (int i=0; i<TempCherriesToMash.size(); i++) {
            CherriesToMash.push_back(TempCherriesToMash[gsl_permutation_get (p,i)]);
        }
        gsl_permutation_free(p);
    }
    else {
        moredecreases=false;
    }

    //Get list of nodes to reroot on
    vector<int> NodesToReRootOn=CurrentTree.GetPotentialNewRoots();
    if (NodesToReRootOn.size()==0) {
        morererootings=false;
    }
    //cout<<"There are potentially "<<NodesToReRootOn.size()<<" new roots\n";


    //Get list of leaves to split; we do this at the start so that we try each possible leaf (leaves with at least two samples) at random but only once.
    vector<int> TempLeavesToSplit;
    vector<int> LeavesToSplit;
    for (int i=0;i<CurrentTree.GetNumLeaves(); i++) {
        int numsamples=0;
        for (int j=0;j<convertsamplestospecies.size();j++) {
            if(convertsamplestospecies[j]==i+1) {
                numsamples++;
            }
        }
        if (numsamples>minsamplesperspecies) {
            TempLeavesToSplit.push_back(i+1);
        }
    }
	if (TempLeavesToSplit.size()>0) {
		gsl_permutation * q = gsl_permutation_alloc (TempLeavesToSplit.size());
		gsl_permutation_init (q);
		gsl_ran_shuffle (r, q->data, TempLeavesToSplit.size(), sizeof(size_t));
		for (int i=0; i<TempLeavesToSplit.size(); i++) {
			LeavesToSplit.push_back(TempLeavesToSplit[gsl_permutation_get (q,i)]);
		}
		gsl_permutation_free(q);
		if (showtries) {
			cout<<"Made LeavesToSplitVector of size "<<LeavesToSplit.size()<<endl<<"contents: ";
			for (int k=0;k<LeavesToSplit.size();k++) {
					cout<<LeavesToSplit[k]<<"\t";
			}
			cout<<endl;
		}
	}
	else {
		moreincreases=false;
	}

    SamplesToMove.clear();
    vector<int> TempSamplesToMove;
    int maxsamplesperspecies=0;
    vector<int> SamplesPerSpecies(1+CurrentTree.GetNumLeaves(),0); //so SamplesPerSpecies[0] is empty but then SamplesPerSpecies[X] is the number of samples for species X
    for (int i=0;i<convertsamplestospecies.size();i++) {
        SamplesPerSpecies[(convertsamplestospecies[i])]++;
        maxsamplesperspecies=GSL_MAX(maxsamplesperspecies,SamplesPerSpecies[(convertsamplestospecies[i])]);
    }
    if (maxsamplesperspecies<=minsamplesperspecies || movefreqvector[1]==0) {
        morereassignments=false;
    }
    else {
        vector<int> TempSampleDestinations;
        for (int i=0;i<CladeVector.size();i++) {
            TempSamplesToMove.push_back(i);
        }

        gsl_permutation * u = gsl_permutation_alloc (TempSamplesToMove.size());
        gsl_permutation_init (u);
        gsl_ran_shuffle (r, u->data, TempSamplesToMove.size(), sizeof(size_t));
        for (int i=0; i<TempSamplesToMove.size(); i++) {
            SamplesToMove.push_back(TempSamplesToMove[gsl_permutation_get (u,i)]);
        }
        gsl_permutation_free(u);

        SampleDestinations.clear();
        bool enoughdestinations=false;
        while (!enoughdestinations && SamplesToMove.size()>0) {
            SampleDestinations.clear();

            for (int i=0;i<CurrentTree.GetNumLeaves();i++) {
                if(TestMoveSamples(SamplesToMove.back(),i+1)) {
                    TempSampleDestinations.push_back(i+1);
                    enoughdestinations=true;
                }
            }
            if (!enoughdestinations) {
                SamplesToMove.pop_back(); //we're done trying to move that sample
            }
        }
        morereassignments=enoughdestinations;
        if(enoughdestinations) {
            gsl_permutation * v = gsl_permutation_alloc (TempSampleDestinations.size());
            gsl_permutation_init (v);
            gsl_ran_shuffle (r, v->data, TempSampleDestinations.size(), sizeof(size_t));
            for (int i=0; i<TempSampleDestinations.size(); i++) {
                SampleDestinations.push_back(TempSampleDestinations[gsl_permutation_get (v,i)]);
            }
            gsl_permutation_free(v);
        }
    }
    //cout<<"moreswaps = "<<moreswaps<<" morereassignments = "<<morereassignments<<" moreincreases = "<<moreincreases<<" moredecreases = "<<moredecreases<<" morererootings = "<<morererootings<<endl;

    if (movecount==0) {

        message="";
		if (jackknifesearch) {
			message+=jackrep;
			message+="\t";
		}
        message+=replicate;
        message+="\t";
        message+=movecount;
        message+="\t";
        message+=CurrentTree.GetNumLeaves();
        message+="\t\t";
        message+=scoretype;
        char outputstring[9];
        sprintf(outputstring,"%9.3f",bestscorelocal);
        message+=outputstring;
        message+="\t";
        sprintf(outputstring,"%9.3f",bestscorelocalvector[1]);
        message+=outputstring;
        message+="\t";
        sprintf(outputstring,"%9.3f",bestscorelocalvector[2]);
        message+=outputstring;
        message+="\t";
        sprintf(outputstring,"%9.3f",bestscorelocal);
        message+=outputstring;
        message+="\t";
        sprintf(outputstring,"%9.3f",GSL_MIN(bestscore,bestscorelocal));
        message+=outputstring;
        message+="\t";
        message+=int(FormattedBestTrees.size());
        if (moreswaps) {
            message+="\ts";
        }
        else {
            message+="\t_";
        }
        if (morereassignments) {
            message+="a";
        }
        else {
            message+="_";
        }
        if (moreincreases) {
            message+="i";
        }
        else {
            message+="_";
        }
        if (moredecreases) {
            message+="d";
        }
        else {
            message+="_";
        }
        if (morererootings) {
            message+="r";
        }
        else {
            message+="_";
        }
        if (status) {
            PrintMessage();
        }
    }

    while ((moreswaps || morereassignments || moreincreases || moredecreases || morererootings) && (rearrlimit<0 || movecount<rearrlimit)) {
        bool somethinghappened=true;
        ContainingTree NextTree=CurrentTree;
        //   cout<<"\n\nOldTree\n"<<ReturnFinalSpeciesTree(CurrentTree)<<endl;
        // for (int i=0;i<convertsamplestospecies.size();i++) {
        //      cout<<convertsamplestospecies[i]<<" ";
        //  }
        // cout<<endl;
        NextTree.UpdateCherries();
        NextTree.Update();
        NextTree.GetNodeDepths();
        vector<int> Originalconvertsamplestospecies=convertsamplestospecies;

        //decide chosen move
        if (CurrentTree.GetNumLeaves()==1) {
            moreswaps=false;
            morereassignments=false;
            moredecreases=false;
            morererootings=false;
        }
        double randomvalue=double(gsl_ran_flat (r,0,1));
        int chosenmove=0;
        nxsstring chosenmovestring="?";
        vector<double> possiblemovefreqvector;
        vector<int> possiblemovechoicevector;
        vector<nxsstring> possiblemoveabbrevvector;
        if (moreswaps) {
            possiblemovefreqvector.push_back(movefreqvector[0]);
            possiblemovechoicevector.push_back(1);
            possiblemoveabbrevvector.push_back("s");
        }
        if (morereassignments) {
            possiblemovefreqvector.push_back(movefreqvector[1]);
            possiblemovechoicevector.push_back(2);
            possiblemoveabbrevvector.push_back("a");
        }
        if (moreincreases) {
            possiblemovefreqvector.push_back(movefreqvector[2]);
            possiblemovechoicevector.push_back(3);
            possiblemoveabbrevvector.push_back("i");
        }
        if (moredecreases) {
            possiblemovefreqvector.push_back(movefreqvector[3]);
            possiblemovechoicevector.push_back(4);
            possiblemoveabbrevvector.push_back("d");
        }
        if (morererootings) {
            possiblemovefreqvector.push_back(movefreqvector[4]);
            possiblemovechoicevector.push_back(5);
            possiblemoveabbrevvector.push_back("r");
        }
        double sumofpossiblemovefreqs=0;
        for (int k=0; k<possiblemovefreqvector.size(); k++) {
            sumofpossiblemovefreqs+=possiblemovefreqvector[k];
        }
        for (int k=0; k<possiblemovefreqvector.size(); k++) {
            possiblemovefreqvector[k]=(possiblemovefreqvector[k])/sumofpossiblemovefreqs;
            //  cout<<"possiblemovefreqvector["<<k<<"] = "<<possiblemovefreqvector[k]<<"\tpossiblemovechoicevector["<<k<<"] = "<<possiblemovechoicevector[k]<<endl;
        }
        double runningtotal=0;
        for (int k=0; k<possiblemovefreqvector.size(); k++) {
            runningtotal+=possiblemovefreqvector[k];
            // cout<<"randomvalue = "<<randomvalue<<" runningtotal = "<<runningtotal;
            if (randomvalue<=runningtotal) {
                chosenmove=possiblemovechoicevector[k];
                chosenmovestring=possiblemoveabbrevvector[k];
                //cout<<" chosenmove is "<<chosenmove;
                break;
            }
            //cout<<endl;
        }
        //cout<<"randomvalue is "<<randomvalue<<" chosenmove is "<<chosenmove<<" moreswaps = "<<moreswaps<<" morereassignments = "<<morereassignments<<" moreincreases = "<<moreincreases<<" moredecreases = "<<moredecreases<<" morererootings = "<<morererootings<<endl;
        movecount++;
        //  cout<<"convertsamplestospecies\n";
        //  for (int m=0;m<convertsamplestospecies.size();m++) {
        //      cout<<convertsamplestospecies[m]<<" ";
        //  }
        //  cout<<endl;
        //cout<<"chosenmove = "<<chosenmove<<endl;
        if (chosenmove==1) { //Try branch swap
            if (showtries) {
                cout<<"Trying branch swap"<<endl;
                cout<<"Start tree = \n";
                NextTree.Draw(cout);
            }
            NextTree.SetBreakVector(CurrentTree.GetBreakVector());
            NextTree.SetAttachVector(CurrentTree.GetAttachVector());
            NextTree.FindAndSetRoot();
            NextTree.Update();
            moreswaps=NextTree.NextSPR();
            if (showtries) {
                cout<<"Swap tree = \n";
                NextTree.Draw(cout);
            }
            CurrentTree.SetBreakVector(NextTree.GetBreakVector()); //due to how the vectors are updated during a swap.
            CurrentTree.SetAttachVector(NextTree.GetAttachVector());
            nextscorevector=GetCombinedScore(&NextTree);
            nextscore=nextscorevector[0];
        }
        else if(chosenmove==2) {
            if (showtries) {
                cout<<"Moving a sample from one species to another\n";
            }
            if (showtries) {
                cout<<"Start assignment = (";
                for (int i=0;i<(CladeVector[SamplesToMove.back()]).size();i++) {
                    cout<<" "<<CladeVector[SamplesToMove.back()][i];
                }
                cout<<" )\n";
                for (int i=0; i<convertsamplestospecies.size();i++) {
                    cout<<convertsamplestospecies[i]<<" ";
                }
                cout<<"\n";
            }
            morereassignments=MoveSamples(Originalconvertsamplestospecies);
            if (showtries) {
                for (int i=0; i<convertsamplestospecies.size();i++) {
                    cout<<convertsamplestospecies[i]<<" ";
                }
                cout<<"\n";
            }
            nextscorevector=GetCombinedScore(&NextTree);
            nextscore=nextscorevector[0];
        }
        else if(chosenmove==3) {
            //cout<<"LeavesToSplitVect\n";
            //for (int i=0;i<LeavesToSplit.size();i++) {
            //     cout<<" "<<LeavesToSplit[i];
            // }
            //increase the number of species
            if (showtries) {
				cout<<"LeavesToSplitVect\n";
				cout<<"vector size is "<<LeavesToSplit.size()<<endl;
				for (int i=0;i<LeavesToSplit.size();i++) {
					cout<<" "<<LeavesToSplit[i];
				}
                cout<<"Splitting a leaf\n";
            }
            int ChosenLeaf=LeavesToSplit.back();
            //cout<<"\nChosenLeaf = "<<ChosenLeaf<<endl;
            LeavesToSplit.pop_back();
            // cout<<"\nVector size now "<<LeavesToSplit.size();
            if (LeavesToSplit.size()==0) {
                moreincreases=false;
            }
            // cout<<"\nmoreincreases value = "<<moreincreases<<endl;
            vector<int>changevector;
            if (showtries) {
                cout<<"Split leaf start tree = \n";
                NextTree.Draw(cout);
            }
            changevector=NextTree.SplitLeaf(ChosenLeaf); //first element is split taxon, second element is new taxon
            if (showtries) {
                cout<<"Final tree = \n";
                NextTree.Draw(cout);
            }
            //now try optimizing the new assignments (which descendant the samples go with) before actually getting the score
            vector<int> samplestomove;
            int sampletostay;
            for (int i=0;i<convertsamplestospecies.size();i++) {
                //cout<<convertsamplestospecies[i]<<" ";
                if(convertsamplestospecies[i]==changevector[0]) {
                    samplestomove.push_back(i);
                }
            }
            // cout<<endl;
            sampletostay=samplestomove.back();
            samplestomove.pop_back();
            //so idea here is to try all combinations, with one sample fixed in the old species and the others allowed to be in either species
            vector<int> Startingconvertsamplestospecies=convertsamplestospecies;
            vector<int> Bestconvertsamplestospecies=convertsamplestospecies;
            bool bestscorefound=false;
            double bestscoreforcombination=GSL_POSINF;
			vector<double> lastscorevector;
            //  cout<<"leaf split starting assignment"<<endl;
            //  for (int i=0;i<convertsamplestospecies.size();i++) {
            //      cout<<convertsamplestospecies[i]<<" ";
            //   }
            //   cout<<endl;
            size_t j;
            double numberofcomparisons=0;
            for(int l=1;l<=samplestomove.size();l++){
                numberofcomparisons+=gsl_sf_choose(samplestomove.size(),l);
				if (showtries) {
					cout<<"numberofcomparisons now "<<numberofcomparisons<<endl;
				}
				if (numberofcomparisons<1) { 
					cout<<"Number of comparisons was "<<numberofcomparisons<<" and samplestomove.size() was "<<samplestomove.size()<<endl;
				}
				assert(numberofcomparisons>0);
            }
            double probofacomb=(pow((1.0*numberofcomparisons),1.0/chosensubsampling))/(1.0*numberofcomparisons); //a way to reduce the search effort
			if (probofacomb*numberofcomparisons<10) { //so that if we choose a ridiculous number we expect to do at least ten swaps
				probofacomb=1;
			}
			else if (probofacomb!=probofacomb) {
				if (showtries) {
					cout<<"prob of a comb is "<<probofacomb<<" so we're adjusting it"<<endl;
				}
				probofacomb=GSL_MIN(100.0/numberofcomparisons,1);
			}
			if (showtries) {
				cout<<"Prob of a comb is "<<probofacomb<<" expected number of assignments to examine is "<<probofacomb*numberofcomparisons<<endl;
			}
            int combinationmoves=0;
            //while(bestscorefound==false) {
			gsl_combination *c;
			for (j=1;j<=samplestomove.size();j++) { //always move
				c=gsl_combination_calloc(samplestomove.size(),j);
				do
				{
					combinationmoves++;
					if (gsl_ran_bernoulli(r,probofacomb)==1 || steepest || exhaustive) {
						if (showtries) {
							
							cout<<combinationmoves<<"/"<<numberofcomparisons<<" = ";
							cout<<(1.0*combinationmoves)/(1.0*numberofcomparisons)<<endl;
                                ProgressBar(0);
							
						}
						convertsamplestospecies=Startingconvertsamplestospecies;
						for (int k=0;k<j;k++) {
							convertsamplestospecies[samplestomove[int(gsl_combination_get(c,k))]]=changevector[1]; //assign this taxon to the new species
						}
                            //  for (int m=0;m<convertsamplestospecies.size();m++) {
                            //      cout<<convertsamplestospecies[m]<<" ";
                            //   }
                            //for (int i=0;i<convertsamplestospecies.size();i++) {
                            //     cout<<convertsamplestospecies[i]<<" ";
                            // }
                            // cout<<endl;
                            //    for (int i=0;i<convertsamplestospecies.size();i++) {
                            //      cout<<convertsamplestospecies[i]<<" ";
                            //   }
                            //  cout<<endl;
						vector<double> newscoreforcombinationvector=GetCombinedScore(&NextTree);
						double newscoreforcombination=newscoreforcombinationvector[0];
						lastscorevector.swap(newscoreforcombinationvector);
						if (newscoreforcombination<bestscoreforcombination) {
							bestscoreforcombination=newscoreforcombination;
							Bestconvertsamplestospecies=convertsamplestospecies;
							if (!steepest && !exhaustive) {
								bestscorefound=true;
							}
						}
                            // cout<<"\t"<<newscoreforcombination<<endl;
					}
				}
				while ((gsl_combination_next (c) == GSL_SUCCESS) && (bestscorefound==false));
				gsl_combination_free(c);
			}
			//gsl_combination_free(c); /moved up to stop leak
           // }
			
			
            convertsamplestospecies.swap(Bestconvertsamplestospecies);
            //  cout<<"Final assignment after split leaf"<<endl;
            //  for (int i=0;i<convertsamplestospecies.size();i++) {
            //      cout<<convertsamplestospecies[i]<<" ";
            //  }
            //  cout<<endl;
            //  cout<<"Bestscorefound = "<<bestscorefound<<endl;
            //convertsamplestospecies=Originalconvertsamplestospecies;
			if(isinf(bestscoreforcombination)==0) {//so the best score is NOT infinity
				nextscorevector=GetCombinedScore(&NextTree);
				nextscore=nextscorevector[0];
			}
			else {
				nextscorevector.swap(lastscorevector);
				nextscore=nextscorevector[0];
			}
        }
        else if(chosenmove==4) {  //reduce the number of species, if possible
                                  //cout<<"\nCherry vector"<<endl;
                                  //for (int i=0;i<CherriesToMash.size(); i++) {
                                  //    cout<<" "<<CherriesToMash[i];
                                  //}
                                  //cout<<endl;
            if (showtries) {
                cout<<"Collapsing a cherry\n";
            }
            int CherryToMash=CherriesToMash.back();
            //cout<<"Cherry to mash = "<<CherryToMash<<endl;
            CherriesToMash.pop_back(); //this is so we look at each cherry once
            if (CherriesToMash.size()==0) {
                moredecreases=false;
            }
            vector<int>changevector;
            if (showtries) {
                cout<<"Cherry collapse start tree = \n";
                NextTree.Draw(cout);
            }
            changevector=NextTree.CollapseCherry(CherryToMash);
            for (int i=0; i<convertsamplestospecies.size(); i++) {
                if(convertsamplestospecies[i]==changevector[1]) {
                    convertsamplestospecies[i]=changevector[0];
                }
                if(convertsamplestospecies[i]>changevector[1]) {
                    convertsamplestospecies[i]--; //so if we have taxa 1-8, and delete taxon 6, taxon 7 becomes the new 6 and taxon 8 becomes the new 7
                }
            }
            if (showtries) {
                cout<<"Next tree = \n";
                NextTree.Draw(cout);
            }
            nextscorevector=GetCombinedScore(&NextTree);
            nextscore=nextscorevector[0];
        }
		
        else if (chosenmove==5) {
            if (showtries) {
                cout<<"Rerooting"<<endl;
            }
            int NodeToReRootOnNum=NodesToReRootOn.back();
            NodesToReRootOn.pop_back(); //this is so we look at each potential position once
            if (NodesToReRootOn.size()==0) {
                morererootings=false;
            }
            vector<int>changevector;
            if (showtries) {
                cout<<"Rerooting start tree = \n";
                NextTree.Draw(cout);
            }
            NextTree.ReRootTree(NextTree.SelectNodeToReRootOn(NodeToReRootOnNum));
            if (showtries) {
                cout<<"Next tree = \n";
                NextTree.Draw(cout);
            }
            nextscorevector=GetCombinedScore(&NextTree);
            nextscore=nextscorevector[0];
        }
        else {
            somethinghappened=false;
            movecount--;
        }
        if (somethinghappened) {
            if (NextTree.GetNumLeaves()==NextTree.GetNumInternals()) {
                errormsg="Error: num leaves = num internals\nLast move chosen was";
                errormsg+=chosenmove;
                NextTree.ReportTreeHealth();
                throw XNexus( errormsg);
				
            }
            assert(CheckConvertSamplesToSpeciesVector(false));
            //int maxspnum=0;
            //for (int i=0;i<convertsamplestospecies.size();i++) {
            //    cout<<convertsamplestospecies[i]<<" ";
            //  if (convertsamplestospecies[i]>maxspnum) {
            //    maxspnum=convertsamplestospecies[i];
            //  }
            //}
            //  cout<<endl;
            //  cout<<"maxspnum="<<maxspnum<<endl;
            // for (int j=1;j<=maxspnum;j++) {
            //     int samplecount=0;
            //      for (int i=0;i<convertsamplestospecies.size();i++) {
            //  if (convertsamplestospecies[i]==j) {
            //        samplecount++;
            //      }
            //    }
            //  cout<<j<<"\t"<<samplecount<<endl;
            //      assert(samplecount>0);
            //    }
			
            //  cout<<"\nNextTree\n"<<ReturnFinalSpeciesTree(NextTree)<<endl;
            scoretype="\t";
            bool modifiedscoretype=false;
            if (nextscore<bestscore) {
                scoretype="*G\t";
                modifiedscoretype=true;
                improvement=true;
                RawBestTrees.clear();
                RawBestTrees.push_back(NextTree);
                FormattedBestTrees.clear();
				TotalScores.clear();
				GTPScores.clear();
				StructScores.clear();
				BestConversions.clear();
			//	TotalScores.push_back(nextscorevector[0]);
			//	GTPScores.push_back(nextscorevector[1]);
			//	StructScores.push_back(nextscorevector[2]);				
                FormatAndStoreBestTree(&NextTree,nextscorevector);
                //  (FormattedBestTrees.back()).Update();
                //  (FormattedBestTrees.back()).GetNodeDepths();
                //  (FormattedBestTrees.back()).Draw(cout);
                bestscore=nextscore;
                if (showtries) {
                    cout<<"GOT BETTER TREE"<<endl<<endl;
                    NextTree.Draw(cout);
                }
            }
            else if (nextscore==bestscore) {
                scoretype="=G\t";
                modifiedscoretype=true;
                NextTree.Update();
                RawBestTrees.push_back(NextTree);
				//TotalScores.push_back(nextscorevector[0]);
				//GTPScores.push_back(nextscorevector[1]);
				//StructScores.push_back(nextscorevector[2]);				
                FormatAndStoreBestTree(&NextTree,nextscorevector);
                //    (FormattedBestTrees.back()).Update();
                //    (FormattedBestTrees.back()).GetNodeDepths();
                //    (FormattedBestTrees.back()).Draw(cout);
            }
			
            if (nextscore<bestscorelocal) {
                if (!modifiedscoretype) {
                    scoretype="*L\t";
                }
                bestscorelocal=nextscore;
                improvement=true;
                BestTreesThisRep.clear();
                NextTree.FindAndSetRoot();
                NextTree.Update();
                BestTreesThisRep.push_back(NextTree);
                CurrentTree.ResetBreakVector(); ///////figure out when to  reset this: any time you move to a new optimum
            }
            else if (nextscore==bestscorelocal) {
                if (!modifiedscoretype) {
                    scoretype="=L\t";
                }
                NextTree.FindAndSetRoot();
                NextTree.Update();
                BestTreesThisRep.push_back(NextTree);
                //cout<<"Swapping back the convertsamplestospecies vector\n";
                //convertsamplestospecies.swap(Originalconvertsamplestospecies); //need to reassign the original one
            }
            else {
                //cout<<"Swapping back the convertsamplestospecies vector\n";
                // convertsamplestospecies.swap(Originalconvertsamplestospecies); //need to reassign the original one
            }
            message="";
			if (jackknifesearch) {
				message+=jackrep;
				message+="\t";
			}
            message+=replicate;
            message+="\t";
            message+=movecount;
            message+="\t";
            message+=CurrentTree.GetNumLeaves();
            message+="->";
            message+=NextTree.GetNumLeaves();
            message+="\t";
            message+=chosenmovestring;
            message+="\t";
            message+=scoretype;
			if (gtptoohigh || triplettoohigh) {
				message+=">";
			}			
            char outputstring[9];
            sprintf(outputstring,"%9.3f",nextscore);
            message+=outputstring;
            message+="\t";
			if (gtptoohigh) {
				message+=">";
			}						
            sprintf(outputstring,"%9.3f",nextscorevector[1]);
            message+=outputstring;
            message+="\t";
			if (gtptoohigh || triplettoohigh) { //since we abort gtp calculations if the triplet cost is already too high
				message+=">";
			}									
            sprintf(outputstring,"%9.3f",nextscorevector[2]);
            message+=outputstring;
            message+="\t";
            sprintf(outputstring,"%9.3f",bestscorelocal);
            message+=outputstring;
            message+="\t";
            sprintf(outputstring,"%9.3f",GSL_MIN(bestscore,bestscorelocal));
            message+=outputstring;
            message+="\t";
            message+=int(FormattedBestTrees.size());
            if (moreswaps) {
                message+="\ts";
            }
            else {
                message+="\t_";
            }
            if (morereassignments) {
                message+="a";
            }
            else {
                message+="_";
            }
            if (moreincreases) {
                message+="i";
            }
            else {
                message+="_";
            }
            if (moredecreases) {
                message+="d";
            }
            else {
                message+="_";
            }
            if (morererootings) {
                message+="r";
            }
            else {
                message+="_";
            }
            if (badgtpcount>0) {
                message+="\t!!!";
                message+=badgtpcount;
                message+="!!!";
            }
            if (status) {
                PrintMessage();
            }
            if ((steepest==false) && improvement) {
                if (showtries) {
                    cout<<"NOW BREAKING..."<<endl;
                }
                break;
            }
            else {
				convertsamplestospecies.assign( Originalconvertsamplestospecies.begin(), Originalconvertsamplestospecies.end() );
                //convertsamplestospecies.swap(Originalconvertsamplestospecies);
            }
        } //if something happened
    } //while (moreswaps || morereassignments || moreincreases || moredecreases )
}//while improvement
 //DelDupes();
 //message="Replicate finished, now removing duplicate trees and saving best to file";
 //PrintMessage();
 //	ofstream outtreef;
 //	outtreef.open(treefilename.c_str());
 //outtreef<<"#nexus\nbegin trees;\n";
 //outtreef<<"[heuristic search, best results after replicate "<<replicate<<"\nSearch options: \n]\n";

//	for (int i=0; i<FormattedBestTrees.size(); i++) {
//		outtreef<<"tree sptre"<<i+1<<" = [&R] ";
//		outtreef<<ReturnFinalSpeciesTree(FormattedBestTrees[i]);
//	       outtreef<<endl;
//	}
//outtreef<<"end;";
//outtreef.close();
//  message="Best trees by the end of rep ";
//    message+=replicate;
//   message+="\n";
//    PrintMessage();
//    for (int i=0; i<FormattedBestTrees.size(); i++) {
//       (FormattedBestTrees[i]).Draw(cout);
//       cout<<endl;
//}
    }//nrep
     //DelDupes();
cout<<endl<<"Best trees overall"<<endl<<endl;
for (int i=0; i<FormattedBestTrees.size(); i++) {
    (FormattedBestTrees[i]).Update();
    (FormattedBestTrees[i]).GetNodeDepths();
    (FormattedBestTrees[i]).Draw(cout);
    cout<<endl;
}
message="\n\n#nexus\nbegin trees;\n";
for (int i=0; i<FormattedBestTrees.size(); i++) {
    message+="tree sptree";
    message+=i+1;
    message+=" = ";
    if (unrooted==1) {
        message+="[&U] ";
    }
    else {
        message+="[&R] ";
    }
    message+=ReturnFinalSpeciesTree(FormattedBestTrees[i]);
    message+="\n";
}
message+="end;\n\n";
PrintMessage();
gsl_matrix_free(TaxonDistance);
gsl_matrix_free(TaxonProportDistance);
}

//gives total score for structure within each putative species
double BROWNIE::GetTripletScore(ContainingTree *SpeciesTreePtr) { //Note that excess structure requires four or more individuals per species
                                                                  //cut each gene tree where it crosses a species tree boundary and reroot on the node connecting to the deleted edge
                                                                  //Watch out: don't compare two trees for the same gene (weighted trees)
                                                                  //Have matrix of expected random triplet scores (3 tax vs 4 tax, etc. ) so that you can subtract this from the observed overlap
    int oldchosentree=chosentree;
    int nspecies=SpeciesTreePtr->GetNumLeaves();
    double totalscore=0; //this is right now the default score if there are no triplets
    for (int i=1; i<=nspecies; i++) {
		if (triplettoohigh) {
			break;
		}
        int nsamplesinspecies=0;
        vector<nxsstring> taxatoexclude;
        for (int j=0; j<convertsamplestospecies.size();j++) {
            if (convertsamplestospecies[j]==i) {
                nsamplesinspecies++;
            }
            else {
                taxatoexclude.push_back(taxa->GetTaxonLabel(j));
                //	cout<<"taxatoexclude "<<taxa->GetTaxonLabel(j)<<endl;
            }
        }
        if (nsamplesinspecies>=3) { //so there's actually a triplet; otherwise, score is ____default____.
                                    //cout<<"at least three species"<<endl;
            double totalweight=0;
            vector<vector <ContainingTree> > GeneTreesVector;
            vector<vector <double> > GeneTreesWeights;
            int numberofgenes=0;
            vector<int> TreesPerGene;
            vector<ContainingTree> OneGeneTreeVector;
            vector<double> OneGeneTreeWeights;
            for (int chosentreenum=0; chosentreenum<trees->GetNumTrees(); chosentreenum++) {
				if (triplettoohigh) {
					break;
				}				
				double newweight=trees->GetTreeWeight(chosentreenum);
				if (!jackknifesearch) {
										//cout<<"Tree weight is "<<newweight<<endl;
					totalweight+=newweight;
					if (totalweight>=1) { //new gene (remember, if we have several bootstrap trees from one gene, we don't want triplet scores against them. Weight of this means that we're starting a new gene
						if (OneGeneTreeVector.size()>0) {
							GeneTreesVector.push_back(OneGeneTreeVector);
							GeneTreesWeights.push_back(OneGeneTreeWeights);
							numberofgenes++;
							TreesPerGene.push_back(OneGeneTreeVector.size());
							OneGeneTreeVector.clear();
						}
						totalweight=newweight;
					}
				}
				else {
					bool usethistree=true;
					if (jackknifevector[i]==0) {
						usethistree=false;
					}
					if (usethistree) {
						if (geneidvector[chosentreenum]>numberofgenes+1) {
							if (OneGeneTreeVector.size()>0) {
								GeneTreesVector.push_back(OneGeneTreeVector);
								GeneTreesWeights.push_back(OneGeneTreeWeights);
								numberofgenes++;
								TreesPerGene.push_back(OneGeneTreeVector.size());
								OneGeneTreeVector.clear();
							}
						}
					}
				}
				bool usethistree=true;
				if (jackknifesearch) {
				if (jackknifevector[i]==0) {
					usethistree=false;
				}
				}
				if (usethistree) {
					//Now we can deal with the gene tree, having re-initialized the vector if need be
					Tree CurrentGeneTreeTreeFmt=intrees.GetIthTree(chosentreenum);
					//	cout<<"Gene tree is "<<endl;
					//	CurrentGeneTreeTreeFmt.Draw(cout);
					ContainingTree CurrentGeneTree;
					CurrentGeneTree.SetRoot(CurrentGeneTreeTreeFmt.CopyOfSubtree(CurrentGeneTreeTreeFmt.GetRoot()));
					vector<ContainingTree> SplitTreeVector;
					SplitTreeVector=CurrentGeneTree.SplitOnTaxon(taxatoexclude);
					for (int k=0; k<SplitTreeVector.size(); k++) {
						OneGeneTreeVector.push_back(SplitTreeVector[k]);
						OneGeneTreeWeights.push_back(newweight);
					}
				}
            }
            if (OneGeneTreeVector.size()>0) { //add the last gene
                GeneTreesVector.push_back(OneGeneTreeVector);
                GeneTreesWeights.push_back(OneGeneTreeWeights);
                numberofgenes++;
                TreesPerGene.push_back(OneGeneTreeVector.size());
                OneGeneTreeVector.clear();
            }
            //now compare all genes to all other genes
            for (int chosengene1=0;chosengene1<numberofgenes-1; chosengene1++) {
				if (triplettoohigh) {
					break;
				}				
                for (int chosengene2=chosengene1+1;chosengene2<numberofgenes; chosengene2++) {
					if (triplettoohigh) {
						break;
					}					
                    for (int chosentreenum1=0; chosentreenum1<TreesPerGene[chosengene1]; chosentreenum1++) {
						if (triplettoohigh) {
							break;
						}						
                        //ContainingTree Tree1=(GeneTreesVector[chosengene1][chosentreenum1]);
                        assert(GeneTreesVector.size()>chosengene1);
                        assert((GeneTreesVector[chosengene1]).size()>chosentreenum1);
                        double Tree1Wt=GeneTreesWeights[chosengene1][chosentreenum1];
                        int Tree1Ntax=(GeneTreesVector[chosengene1][chosentreenum1]).GetNumLeaves();
                        for (int chosentreenum2=0; chosentreenum2<TreesPerGene[chosengene2]; chosentreenum2++) {
							if (triplettoohigh) {
								break;
							}							
                            double Tree2Wt=GeneTreesWeights[chosengene2][chosentreenum2];
                            assert(GeneTreesVector.size()>chosengene2);
                            assert((GeneTreesVector[chosengene2]).size()>chosentreenum2);
                            //	(GeneTreesVector[chosengene2][chosentreenum2]).Write(cout);
                            //	cout<<endl;
                            (GeneTreesVector[chosengene2][chosentreenum2]).Update();
                            int Tree2Ntax=(GeneTreesVector[chosengene2][chosentreenum2]).GetNumLeaves();
                            if ((Tree1Ntax>2) && (Tree2Ntax>2)) { //so, at least a triplet in each
                                ContainingTree t1=GeneTreesVector[chosengene1][chosentreenum1]; //Need to store copies of trees because they're modified (leaves deleted)
                                ContainingTree t2=GeneTreesVector[chosengene2][chosentreenum2];
                                int taxaincommon=PrepareTreesForTriplet(&t1,&t2);
                                if (taxaincommon>2) {
									//cout<<"t1=\n";
									//t1.Draw(cout);
									//cout<<"\nt2=\n";
									//t2.Draw(cout);
                                    vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,taxaincommon);
                                    int maxnumber=tripletoverlapoutput[0];
                                    int numberdisagree=tripletoverlapoutput[1];
                                    int numberunresolved=tripletoverlapoutput[2];
                                    int numberagree=maxnumber-numberdisagree-numberunresolved; //Note that this is the number of triplets resolved IN BOTH TREES that agree
                                    double newscore=ComputeTripletCost(numberagree,maxnumber,taxaincommon,Tree1Wt,Tree1Ntax,Tree2Wt,Tree2Ntax,numberofgenes);
                                    // cout<<"maxnum="<<maxnumber<<" dis="<<numberdisagree<<" un="<<numberunresolved<<" agr="<<numberagree<<" newscore="<<newscore<<endl;
                                    totalscore+=newscore;
									if ((totalscore*structwt)>bestscorelocal) {
										totalscore=(0.0001+bestscorelocal)/structwt;
										triplettoohigh=true;
										break;
									}
                                }
                            }
                        }
                    }
                }
            }
            
        }
        
    }
    chosentree=oldchosentree;
    return totalscore;
}

//computes the TaxonDistance matrix. each entry (i,j) is the number of times taxon i and taxon j are each others' closest relatives in a triplet
void BROWNIE::GetTaxonTaxonTripletDistances() {
	TripletCounts.clear(); //Clear the triplet counts map
	TripleCounts.clear(); //Clears the triples counts map (number of time each set of three taxa co-occur on a tree)
    int nsamples=taxa->GetNumTaxonLabels();
    //	cout<<"nsamples = "<<nsamples<<endl;
    TaxonDistance=gsl_matrix_calloc(nsamples,nsamples); //on diagonal is number of triplets containing the taxon, other elements are as above.
                                                        //	gsl_matrix *TaxonTripletsNotOnSameSide=gsl_matrix_calloc(,);
                                                        //	gsl_matrix *TaxonTripletsTotal=gsl_matrix_calloc(intrees.GetNumTrees(),intrees.GetNumTrees()); //at end, do element by element division, watch for division by zero
    for (int i=0;i<intrees.GetNumTrees();i++) {
			message="now getting triplets on tree ";
			message+=i+1;
			message+=" of ";
			message+=intrees.GetNumTrees();
			PrintMessage();
		bool usethistree=true;
		if (jackknifesearch) {
			if (jackknifevector[i]==0) {
				usethistree=false;
			}
		}
		if (usethistree) {
			Tree t1=intrees.GetIthTree(i);
			gsl_combination *combo;
			combo=gsl_combination_calloc(t1.GetNumLeaves(),3);
			vector<nxsstring> LeafLabelVect; //this section taken from containingtree.cpp (no easy way to convert tree to containing tree)
			NodeIterator <Node> n (t1.GetRoot());
			NodePtr currentnode = n.begin();
			t1.MakeNodeList(); //moved from below
			while (currentnode)
			{
				if(currentnode->IsLeaf()) {
					LeafLabelVect.push_back(currentnode->GetLabel());
				}
				currentnode = n.next();
			}
			do
			{
				NodePtr p;
				NodePtr q;
				nxsstring a=LeafLabelVect[gsl_combination_get(combo,0)];
				nxsstring b=LeafLabelVect[gsl_combination_get(combo,1)];
				nxsstring c=LeafLabelVect[gsl_combination_get(combo,2)];
				std::map<Node *, int, std::less<Node *> > depth;
				PreorderIterator <Node> n (t1.GetRoot());
				int count = 0;
				q = n.begin();
				while (q)
				{
					depth[q] = count++;
					q = n.next();
				}
				
            //Now get depths for each pair
				vector<int> LCADepthVectorT1;
				//t1.MakeNodeList(); //Move this up a level -- calculate it once, see if it goes faster
				p = t1.GetLeafWithLabel(a.c_str());
				q = t1.GetLeafWithLabel(b.c_str());
				while (depth[p] != depth[q])
				{
					if (depth[p] < depth[q])
						q = q->GetAnc();
					else
						p = p->GetAnc();
				}
				LCADepthVectorT1.push_back(depth[p]);
				p = t1.GetLeafWithLabel(b.c_str());
				q = t1.GetLeafWithLabel(c.c_str());
				while (depth[p] != depth[q])
				{
					if (depth[p] < depth[q])
						q = q->GetAnc();
					else
						p = p->GetAnc();
				}
				LCADepthVectorT1.push_back(depth[p]);
				p = t1.GetLeafWithLabel(a.c_str());
				q = t1.GetLeafWithLabel(c.c_str());
				while (depth[p] != depth[q])
				{
					if (depth[p] < depth[q])
						q = q->GetAnc();
					else
						p = p->GetAnc();
				}
				LCADepthVectorT1.push_back(depth[p]);
				int T1DepthMax=0; //root has depth 0, others have higher depths.
				int T1DepthMaxIndex=0;
				int anum=taxa->FindTaxon(a);
				int bnum=taxa->FindTaxon(b);
				int cnum=taxa->FindTaxon(c);
            //cout<<anum<<": "<<a<<" "<<bnum<<": "<<b<<" "<<cnum<<": "<<c<<endl;
				for (int j=0;j<3;j++) {
					if(LCADepthVectorT1[j]>T1DepthMax) {
						T1DepthMax=LCADepthVectorT1[j];
						T1DepthMaxIndex=j;
					}
				}
            //j=0->ab, j=1->bc j=2->ac
				gsl_matrix_set(TaxonDistance,anum,anum,1+gsl_matrix_get(TaxonDistance,anum,anum));
				gsl_matrix_set(TaxonDistance,bnum,bnum,1+gsl_matrix_get(TaxonDistance,bnum,bnum));
				gsl_matrix_set(TaxonDistance,cnum,cnum,1+gsl_matrix_get(TaxonDistance,cnum,cnum));
				nxsstring tripletlabel="";
				if (T1DepthMaxIndex==0) {
					gsl_matrix_set(TaxonDistance,anum,bnum,1+gsl_matrix_get(TaxonDistance,anum,bnum));
					gsl_matrix_set(TaxonDistance,bnum,anum,gsl_matrix_get(TaxonDistance,anum,bnum));
					tripletlabel+=GSL_MIN(anum,bnum);
					tripletlabel+="_";
					tripletlabel+=GSL_MAX(anum,bnum);
					tripletlabel+="_";
					tripletlabel+=cnum;
					TripletCounts[tripletlabel]++;
				}
				else if (T1DepthMaxIndex==1) {
					gsl_matrix_set(TaxonDistance,bnum,cnum,1+gsl_matrix_get(TaxonDistance,bnum,cnum));
					gsl_matrix_set(TaxonDistance,cnum,bnum,gsl_matrix_get(TaxonDistance,bnum,cnum));
					tripletlabel+=GSL_MIN(cnum,bnum);
					tripletlabel+="_";
					tripletlabel+=GSL_MAX(cnum,bnum);
					tripletlabel+="_";
					tripletlabel+=anum;
					TripletCounts[tripletlabel]++;
				}
				else if (T1DepthMaxIndex==2) {
					gsl_matrix_set(TaxonDistance,anum,cnum,1+gsl_matrix_get(TaxonDistance,anum,cnum));
					gsl_matrix_set(TaxonDistance,cnum,anum,gsl_matrix_get(TaxonDistance,anum,cnum));
					tripletlabel+=GSL_MIN(cnum,anum);
					tripletlabel+="_";
					tripletlabel+=GSL_MAX(cnum,anum);
					tripletlabel+="_";
					tripletlabel+=bnum;
					TripletCounts[tripletlabel]++;
				}
				int taxonarray[3]={anum, bnum, cnum};
				sort(taxonarray,taxonarray+3);
				nxsstring triplelabel="";
				triplelabel+=taxonarray[0];
				triplelabel+="_";
				triplelabel+=taxonarray[1];
				triplelabel+="_";
				triplelabel+=taxonarray[2];
				TripleCounts[triplelabel]++;
				//if (debugmode) {
				//	cout<<"just got triplet for taxa "<<anum<<", "<<bnum<<", and "<<cnum<<endl;
				//}
			}
			while (gsl_combination_next(combo) == GSL_SUCCESS);
			gsl_combination_free(combo);
		}
    }
    TaxonProportDistance=gsl_matrix_calloc(nsamples,nsamples);
	double maxnumbertriplets=0;
for (int i=0;i<nsamples;i++) {
	for (int j=0;j<nsamples;j++) {
		if (i!=j) {
			maxnumbertriplets=GSL_MAX(gsl_matrix_get(TaxonDistance,i,j),maxnumbertriplets);
		}
	}
}
maxnumbertriplets++; //This is to ensure that taxa are closer to themselves than to other taxa
    //double maxnumbertriplets=gsl_matrix_max(TaxonDistance);
if (debugmode) {
	cout<<"Original similarity matrix\n";
}
for (int i=0;i<nsamples;i++) {
	if (debugmode) {
		cout<<endl<<i<<"\t";
	}
	for (int j=0;j<nsamples;j++) {
		if (debugmode) {
			cout<<gsl_matrix_get(TaxonDistance,i,j)<<"\t";
		}
	}
}
if (debugmode) {
	cout<<endl<<"\nNew distance matrix\n";
}



    for (int i=0;i<nsamples;i++) {
		if (debugmode) {
				cout<<endl<<i<<"\t";
		}
        for (int j=0;j<nsamples;j++) {
			if (i==j) {
				gsl_matrix_set(TaxonProportDistance,i,j,0.0);
			}
			else {
				gsl_matrix_set(TaxonProportDistance,i,j,1-((1.0*gsl_matrix_get(TaxonDistance,i,j))/maxnumbertriplets));
			}
            //gsl_matrix_set(TaxonProportDistance,i,j,1-((1.0*gsl_matrix_get(TaxonDistance,i,j))/maxnumbertriplets));
			if (debugmode) {
					cout<<gsl_matrix_get(TaxonProportDistance,i,j)<<"\t";
			}
        }
    }
	if (debugmode) {
		cout<<endl;
	}
    //for (int i=0;i<nsamples;i++) {
    //		cout<<endl;
    //		cout<<i<<": "<<taxa->GetTaxonLabel(i)<<"\t";
    //		for (int j=0;j<nsamples;j++) {
    //			cout<<"\t"<<gsl_matrix_get(TaxonProportDistance,i,j);
    //		}
    //	}
    //for (int i=0;i<nsamples;i++) {
    //	cout<<endl;
    //	cout<<i<<": "<<taxa->GetTaxonLabel(i)<<"\t";
    //	for (int j=0;j<nsamples;j++) {
    //		cout<<"\t"<<gsl_matrix_get(TaxonDistance,i,j);
    //	}
    //}
    //	cout<<"\nTaxondistance dimensions = "<<TaxonDistance->size1<<","<<TaxonDistance->size2<<endl;
    //	cout<<endl<<"last element = ";
    //	cout<<gsl_matrix_get(TaxonDistance,-1+(TaxonDistance->size1),-1+(TaxonDistance->size1));
}

ContainingTree BROWNIE::ComputeTripletNJTree() {
    ContainingTree NJresult;
    vector<NodePtr> NewNodeVector;
    vector<bool> UnattachedNodeVector;
    vector<int> MatrixToNodeVector;
    int nsamples=taxa->GetNumTaxonLabels();
    int nsamplesremaining=nsamples;
    gsl_matrix *dmatrix=gsl_matrix_calloc(nsamples,nsamples);
    for (int i=0;i<nsamples;i++) {
        NodePtr p=NJresult.NewNode();
        p->SetLabel(taxa->GetTaxonLabel(i));
        p->SetLeaf(true);
        p->SetLeafNumber(i+1);
        p->SetWeight(0);
        NJresult.SetEdgeLengths(true);
        NJresult.SetRooted(true);
        NewNodeVector.push_back(p);
        UnattachedNodeVector.push_back(true); //true means unattached
        MatrixToNodeVector.push_back(i);
        for (int j=0;j<nsamples;j++) {
            gsl_matrix_set(dmatrix,i,j,gsl_matrix_get(TaxonProportDistance,i,j));
        }
    }
    while (nsamplesremaining>2) {
        //cout<<"\n\nMatrixToNodeVector=";
        ////for (int i=0;i<nsamplesremaining;i++) {
        //	cout<<" "<<MatrixToNodeVector[i];
        //}
        //cout<<"\n\nNewNodeVector=";
        //for (int i=0;i<NewNodeVector.size();i++) {
        //	cout<<" "<<NewNodeVector[i]<<"("<<(NewNodeVector[i])->GetLabel()<<")";
        //}
        //for (int i=0;i<nsamplesremaining;i++) {
        //	cout<<endl<<i<<": ";
        //	for (int j=0;j<nsamplesremaining;j++) {
        //		cout<<"\t"<<gsl_matrix_get(dmatrix,i,j);
        //	}
        //}
        //cout<<endl;
gsl_vector *rvector=gsl_vector_calloc(nsamplesremaining);
for (int i=0;i<nsamplesremaining;i++) {
    for (int j=0;j<nsamplesremaining;j++) {
        if(i!=j) {
            gsl_vector_set(rvector,i,(gsl_vector_get(rvector,i))+gsl_matrix_get(dmatrix,i,j));
        }
    }
}
int imin;
int jmin;
double Mijmin=GSL_POSINF;
for (int i=0;i<nsamplesremaining;i++) {
    for (int j=i+1;j<nsamplesremaining;j++) {
        double Mij=gsl_matrix_get(dmatrix,i,j)-((gsl_vector_get(rvector,i)+gsl_vector_get(rvector,j))/(nsamplesremaining-2.0));
        if (Mij<Mijmin) {
            Mijmin=Mij;
            imin=i;
            jmin=j;
        }
    }
}
NodePtr u=NJresult.NewNode();
u->SetWeight(0);
NodePtr i=NewNodeVector[MatrixToNodeVector[imin]];
NodePtr j=NewNodeVector[MatrixToNodeVector[jmin]];
//cout<<"Now joining imin: "<<imin<<": "<<i<<"("<<i->GetLabel()<<") and jmin: "<<jmin<<": "<<j<<"("<<j->GetLabel()<<") to make "<<u<<endl;
u->SetLeaf(false);
u->SetChild(i);
i->SetAnc(u);
i->SetSibling(j);
j->SetAnc(u);
gsl_matrix *dmatrixnew=gsl_matrix_calloc(nsamplesremaining-1,nsamplesremaining-1);
gsl_matrix *dmatrixorig=gsl_matrix_calloc(nsamplesremaining,nsamplesremaining);
gsl_matrix_memcpy(dmatrixorig,dmatrix);
i->SetEdgeLength(gsl_matrix_get(dmatrix,imin,jmin)/2.0+(1.0*gsl_vector_get(rvector,imin)-1.0*gsl_vector_get(rvector,jmin))/(2.0*(1.0*nsamplesremaining-2.0)));
j->SetEdgeLength(gsl_matrix_get(dmatrix,imin,jmin)-(i->GetEdgeLength()));
NewNodeVector.push_back(u);
UnattachedNodeVector.push_back(true);
UnattachedNodeVector[MatrixToNodeVector[imin]]=false;
UnattachedNodeVector[MatrixToNodeVector[jmin]]=false;
bool seenj=false;
for (int i=0;i<MatrixToNodeVector.size();i++) {
    if (seenj==false) {
        if (i==jmin) {
            seenj=true;
        }
        else if (i==imin) {
            MatrixToNodeVector[i]=-1+NewNodeVector.size();
        }
    }
    else {
        if (i==imin) {
            MatrixToNodeVector[i-1]=-1+NewNodeVector.size();
        }
        else {
            MatrixToNodeVector[i-1]=MatrixToNodeVector[i];
        }
    }
}
MatrixToNodeVector.pop_back();
for (int i=0;i<nsamplesremaining;i++) {
    for (int j=0;j<nsamplesremaining;j++) {
        if(i==imin) {
            double newvalue=0.0;
            newvalue+=gsl_matrix_get(dmatrixorig,j,imin);
            newvalue+=gsl_matrix_get(dmatrixorig,j,jmin);
            newvalue-=gsl_matrix_get(dmatrixorig,imin,jmin);
            newvalue*=0.5;
            //			cout<<"changed pos ("<<i<<","<<j<<") from value "<<gsl_matrix_get(dmatrix,i,j)<<" or "<<gsl_matrix_get(dmatrixorig,i,j);
            gsl_matrix_set(dmatrix,i,j,newvalue);
            //			cout<<" to have value "<<gsl_matrix_get(dmatrixorig,j,imin)<<"+"<<gsl_matrix_get(dmatrixorig,j,jmin)<<"-"<<gsl_matrix_get(dmatrixorig,imin,jmin)<<"="<<gsl_matrix_get(dmatrix,i,j)<<endl;
        }
        if(j==imin) {
            double newvalue=0.0;
            newvalue+=gsl_matrix_get(dmatrixorig,imin,i);
            newvalue+=gsl_matrix_get(dmatrixorig,jmin,i);
            newvalue-=gsl_matrix_get(dmatrixorig,imin,jmin);
            newvalue*=0.5;
            //			cout<<"changed pos ("<<i<<","<<j<<") from value "<<gsl_matrix_get(dmatrix,i,j)<<" or "<<gsl_matrix_get(dmatrixorig,i,j);
            gsl_matrix_set(dmatrix,i,j,newvalue);
            //			cout<<" to have value "<<gsl_matrix_get(dmatrixorig,imin,i)<<"+"<<gsl_matrix_get(dmatrixorig,jmin,i)<<"-"<<gsl_matrix_get(dmatrixorig,imin,jmin)<<"="<<gsl_matrix_get(dmatrix,i,j)<<endl;
        }
    }
}
gsl_matrix_set(dmatrix,imin,imin,0);
gsl_matrix_set(dmatrix,jmin,jmin,0);
gsl_matrix_set(dmatrix,imin,jmin,0);
gsl_matrix_set(dmatrix,jmin,imin,0);
int rowoffset=0;
gsl_matrix *dmatrixintermed=gsl_matrix_calloc(nsamplesremaining-1,nsamplesremaining);
for (int i=0;i<nsamplesremaining;i++) {
    if (i==jmin) {
        rowoffset=1;
    }
    else {
        for (int j=0;j<nsamplesremaining;j++) {
            gsl_matrix_set(dmatrixintermed,i-rowoffset,j,gsl_matrix_get(dmatrix,i,j));
        }
    }
}
//cout<<"\ndintermed\n";
//for (int i=0;i<nsamplesremaining-1;i++) {
//	cout<<endl<<i<<":\t";
//	for (int j=0;j<nsamplesremaining;j++) {
//		cout<<"\t"<<gsl_matrix_get(dmatrixintermed,i,j);
//	}
//}
for (int i=0;i<nsamplesremaining-1;i++) {
    int coloffset=0;
    for (int j=0;j<nsamplesremaining;j++) {
        if (j==jmin) {
            coloffset=1;
        }
        else {
            gsl_matrix_set(dmatrixnew,i,j-coloffset,gsl_matrix_get(dmatrixintermed,i,j));
        }
    }
}

nsamplesremaining--;
gsl_vector_free(rvector);
gsl_matrix_free(dmatrix);
dmatrix=gsl_matrix_calloc(nsamplesremaining,nsamplesremaining);
gsl_matrix_memcpy(dmatrix,dmatrixnew);
gsl_matrix_free(dmatrixintermed);
gsl_matrix_free(dmatrixnew);
gsl_matrix_free(dmatrixorig);
    }
//now two nodes left to attach
NodePtr u=NJresult.NewNode();
u->SetWeight(0);
vector<NodePtr> NodesToConnect;
for(int i=0;i<UnattachedNodeVector.size();i++) {
    if(UnattachedNodeVector[i]) {
        NodesToConnect.push_back(NewNodeVector[i]);
    }
}
assert(NodesToConnect.size()==2);
NodePtr i=NodesToConnect[0];
NodePtr j=NodesToConnect[1];
u->SetLeaf(false);
u->SetChild(i);
i->SetAnc(u);
i->SetSibling(j);
j->SetAnc(u);
i->SetEdgeLength(0.5*gsl_matrix_get(dmatrix,0,0));
j->SetEdgeLength(0.5*gsl_matrix_get(dmatrix,0,0));
for (int i=0;i<NewNodeVector.size();i++) {
    NodePtr p=NewNodeVector[i];
    //cout<<"p = "<<p;
    //cout<<" "<<p->GetLabel();
    //cout<<" Anc "<<p->GetAnc()<<endl;
}
NJresult.SetRoot(u);
//cout<<"Root is "<<NJresult.GetRoot()<<endl;
//NJresult.ReportTreeHealth();
//NJresult.GetNodeDepths();
NJresult.Update();
//	NJresult.Draw(cout);
//	cout<<endl;
//	NJresult.Write(cout);
//	NJresult.SetPathLengths();
//	double minpathlength=GSL_POSINF;
//	double maxpathlength=GSL_NEGINF;
//	NodePtr MaxLengthNode;
//	NodePtr MinLengthNode;
//	NodeIterator <Node> n (NJresult.GetRoot());
//	NodePtr currentnode = n.begin();
//	while (currentnode) {
//		if (currentnode->IsLeaf()) {
//			if (currentnode->GetPathLength()<minpathlength) {
//				minpathlength=currentnode->GetPathLength();
//				MinLengthNode=currentnode;
//			}
//			else if (currentnode->GetPathLength()>maxpathlength) {
//				maxpathlength=currentnode->GetPathLength();
//				MaxLengthNode=currentnode;
//			}
//		}
//		currentnode = n.next();
//	}
//	//now, midpoint root the tree
//	double desiredpathlength=(maxpathlength+minpathlength)/2.0;
//	currentnode=MaxLengthNode;
//	while (currentnode->GetPathLength()>desiredpathlength) {
//		if (((currentnode->GetAnc())->GetPathLength())<=desiredpathlength) {
//			NodePtr newroot=NJresult.ReRootTree(currentnode);
//			break;
//		}
//		currentnode=currentnode->GetAnc();
//	}
vector<int> origconvertsamplestospecies;
origconvertsamplestospecies.swap(convertsamplestospecies);
convertsamplestospecies.clear();
for (int i=0;i<nsamples;i++) {
    convertsamplestospecies.push_back(i+1);
}
ContainingTree NJrenamed=NJresult;
nxsstring NewLeafLabel;
NodeIterator <Node> n (NJrenamed.GetRoot());
NodePtr currentnode = n.begin();
while (currentnode)
{
    if (currentnode->IsLeaf())
    {
        NewLeafLabel="taxon";
        NewLeafLabel+=currentnode->GetLeafNumber();
        currentnode->SetLabel(NewLeafLabel);
    }
    currentnode = n.next();
}
//cout<<"Renamed tree"<<endl;
NJrenamed.Update();
//NJrenamed.Draw(cout);
//cout<<endl;
//cout<<"Now geting initial GTP score, it =";
NodePtr OldRoot=(NJrenamed.GetRoot())->GetChild(); //since the previous root itself will be deleted, need to reroot on its descendant.
NodePtr tmp=NJrenamed.ReRootTree(OldRoot); //so that the child its sib don't change when we reset
NJrenamed.Update();
double bestGTPscore=GetGTPScoreNew(&NJrenamed);
//cout<<bestGTPscore;
//cout<<endl;
int NodeToReRootOn=0;
//cout<<"Now starting nodeIterator m\n";
NodeIterator <Node> m (NJrenamed.GetRoot());
//cout<<"NOw setting currentnode=m.begin()\n";
NodePtr mcurrentnode = m.begin();
//cout<<"Now in while(currentnode)\n";
vector<NodePtr> nodestorerootNJrenamedon;
while (mcurrentnode) {
    if (mcurrentnode!=(NJrenamed.GetRoot())) {
        nodestorerootNJrenamedon.push_back(mcurrentnode);
    }
    mcurrentnode=m.next();
}
for (int i=0;i<nodestorerootNJrenamedon.size();i++) {
    NodePtr NewRoot=NJrenamed.ReRootTree(nodestorerootNJrenamedon[i]);
    NJrenamed.Update();
    double newGTPscore=GetGTPScoreNew(&NJrenamed);
    if (newGTPscore<bestGTPscore) {
        NodeToReRootOn=i+1;
        bestGTPscore=newGTPscore;
    }
    //  cout<<"FOr rerooting "<<i+1<<endl;
    //  NJrenamed.Draw(cout);
    // cout<<"\nscore is "<<newGTPscore<<endl;
    NewRoot=NJrenamed.ReRootTree(OldRoot); //Reset it to initial tree
}

//cout<<"OPtimal rerooting on node "<<NodeToReRootOn<<endl;

int currentnodenum=0;
convertsamplestospecies.swap(origconvertsamplestospecies);
if (NodeToReRootOn!=0) {
    NodeIterator <Node> o (NJresult.GetRoot());
    currentnode=o.begin();
    while (currentnode) {
        if (currentnode!=(NJresult.GetRoot())) {
            currentnodenum++;
            if (currentnodenum==NodeToReRootOn) {
                NodePtr NewRoot=NJresult.ReRootTree(currentnode);
            }
        }
        currentnode=o.next();
    }
}



NJresult.Update();
message="Tree of samples: ";
PrintMessage();
NJresult.Draw(cout);
message="";
PrintMessage();

//now modify CladeVector;
CladeVector.clear();
CladeVectorNJBrlen.clear(); //This is used to record internal brlen for the triplet tree, which is used to help do the initial assignments
CladeVectorTripletSupport.clear(); //stores the mean internal brlen of the starting nj tree, used for deciding on initial splits
for (int i=0;i<(2*nsamples)-2;i++) {
    vector<int> tempvector;
    CladeVector.push_back(tempvector);
	CladeVectorNJBrlen.push_back(0.0); //Just initializing
	CladeVectorTripletSupport.push_back(0.0); //Just initializing
}
if (debugmode) {
NJresult.ReportTreeHealth();
}
NodeIterator <Node> l (NJresult.GetRoot()); //goes from tips down
currentnode=l.begin();
int nodecount=0;
double brlensum=0.0;
while (currentnode) {
    vector<int> tempvector; //stores the list of terminal taxa as we go down the tree
    if (currentnode->IsLeaf()) {
        int leafnumtrans=(currentnode->GetLeafNumber())-1;
		if (debugmode) {
			cout<<"looking at leaf with node number "<<leafnumtrans<<" and address "<<currentnode<<endl;
		}
        currentnode->SetIndex(leafnumtrans);
        tempvector.push_back(leafnumtrans);
        CladeVector[leafnumtrans]=tempvector;
		CladeVectorNJBrlen[leafnumtrans]=-1.0*(currentnode->GetEdgeLength()); //terminal branches are stored as negative lengths so we can distinguish them later
    }
    else if(currentnode!=NJresult.GetRoot()) {
        int nodenumber=nodecount+NJresult.GetNumLeaves();
		if (debugmode) {
			cout<<"looking at internal with node number "<<nodenumber<<" and address "<<currentnode<<endl;
		}		
        currentnode->SetIndex(nodenumber);
        int desc1=(currentnode->GetChild())->GetIndex();
        int desc2=((currentnode->GetChild())->GetSibling())->GetIndex();
        for (int i=0;i<(CladeVector[desc1]).size();i++) {
            tempvector.push_back(CladeVector[desc1][i]);
        }
        for (int i=0;i<(CladeVector[desc2]).size();i++) {
            tempvector.push_back(CladeVector[desc2][i]);
        }
        CladeVector[nodenumber]=tempvector;
		CladeVectorNJBrlen[nodenumber]=currentnode->GetEdgeLength();
		brlensum+=currentnode->GetEdgeLength();
		

		
        nodecount++;
    }
	
    currentnode=l.next();
}
meaninternalbrlen=brlensum/(1.0*nodecount);
if (debugmode) {
	cout<<endl<<"average = "<<meaninternalbrlen<<endl;
}

//Get triplet support on edges: have to loop through the tree again (first had to assign clade vectors to tree)
NodeIterator <Node> q (NJresult.GetRoot()); //goes from tips down
currentnode=q.begin();
nodecount=0;
while (currentnode) {
    if (currentnode->IsLeaf()) {
        int leafnumtrans=(currentnode->GetLeafNumber())-1;
		if (debugmode) {
			cout<<"looking at leaf with node number "<<leafnumtrans<<" and address "<<currentnode<<endl;
		}
		CladeVectorTripletSupport[leafnumtrans]=0.0; //triplet support for an OTU is meaningless
    }
    else if(currentnode!=NJresult.GetRoot()) {
        int nodenumber=nodecount+NJresult.GetNumLeaves();
		if (debugmode) {
			cout<<"looking at internal with node number "<<nodenumber<<" and address "<<currentnode<<endl;
		}		
		
		//Get triplet support on edges
		/*
		 If this is the tree, and we want to get support for the edge below the MRCA of taxon p1samp2 and p1samp3, 
		 we need the proportion of relevant triplets with (samp2 or samp6) forming a clade with samp3 relative to (samp1 or samp7)
		 (we want to know the support ONLY at that edge, not using triplets relevant for branches higher or lower in the tree).
		 
		 +- p1samp1
							  +--|
		 |  +- p1samp7
		 +-|
		 | |  +- p1samp2
		 | | +|
		 | | |+- p1samp6
		 | +-|
		 |   +-- p1samp3
		 ---------|
		 |    +- p1samp4
		 +----|
		 +- p1samp5
		 
		 basically, for node with two descendant clades A and B and with sister node (clade) C, want the proportion of triplets matching (A,B),C, where the 
		 taxa are pulled from those clades
		 */		 
		
		//Remember, this tree must be binary (thank goodness!)
		vector<int> avector=CladeVector[(currentnode->GetChild())->GetIndex()]; //Get the child of the current node, get its node number, get the list of taxa at that node
		vector<int> bvector=CladeVector[((currentnode->GetChild())->GetSibling())->GetIndex()]; //Get the sibling of the child of the current node, get its node number, get the list of taxa at that node
		vector<int> cvector;
		if ((currentnode->GetSibling())!=NULL) {
			cvector=CladeVector[(currentnode->GetSibling())->GetIndex()]; //get list of taxa in sister node
			if (debugmode) {
				cout<<"Sister node = "<<currentnode->GetSibling()<<endl;
			}
			
		}
		else { //This node is the sibling, so have to go down one step and back up to get its actual sib
			cvector=CladeVector[((currentnode->GetAnc())->GetChild())->GetIndex()];
			if (debugmode) {
				cout<<"Sister node = "<<(currentnode->GetAnc())->GetChild()<<endl;
			}
		}
		int numberofmatchingtriplets=0;
		int numberoftriplespossible=0;
		if (debugmode) {
			cout<<"vector sizes (a, b, c): "<<avector.size()<<" "<<bvector.size()<<" "<<cvector.size()<<endl;
		}
		for (int aindex=0;aindex<avector.size();aindex++) {
			for (int bindex=0;bindex<bvector.size();bindex++) {
				nxsstring tripletstringroot="";
				int anum=avector[aindex];
				int bnum=bvector[bindex];
				tripletstringroot+=GSL_MIN(anum,bnum);
				tripletstringroot+="_";
				tripletstringroot+=GSL_MAX(anum,bnum);
				tripletstringroot+="_";
				for (int cindex=0;cindex<cvector.size();cindex++) {
					nxsstring tripletstring=tripletstringroot;
					int cnum=cvector[cindex];
					tripletstring+=cnum;
					int taxonarray[3]={anum, bnum, cnum};
					sort(taxonarray,taxonarray+3);
					nxsstring triplelabel="";
					triplelabel+=taxonarray[0];
					triplelabel+="_";
					triplelabel+=taxonarray[1];
					triplelabel+="_";
					triplelabel+=taxonarray[2];
					map<nxsstring,int>::iterator iter = TripletCounts.find(tripletstring);
					if( iter != TripletCounts.end() ) {
						numberofmatchingtriplets+=iter->second;
					}
					map<nxsstring,int>::iterator iter2 = TripleCounts.find(triplelabel);
					if( iter != TripleCounts.end() ) {
						numberoftriplespossible+=iter2->second;
					}					
					if (debugmode) {
						cout<<"Triplet = "<<tripletstring<<" triple = "<<triplelabel<<endl;
					}
				}
			}
		}
		if (debugmode) {
			cout<<"number of matching triplets = "<<numberofmatchingtriplets<<" number of triplets possible = "<<numberoftriplespossible;
		}
		CladeVectorTripletSupport[nodenumber]=(1.0*numberofmatchingtriplets)/(1.0*numberoftriplespossible);
		assert (CladeVectorTripletSupport[nodenumber]<=1 && CladeVectorTripletSupport[nodenumber]>=0);
		if (debugmode) {
			cout<<" proport = "<<CladeVectorTripletSupport[nodenumber]<<endl;
		}
        nodecount++;
    }
	
    currentnode=q.next();
}






if (debugmode) {
	
	for( map<nxsstring, int>::iterator iter = TripletCounts.begin(); iter != TripletCounts.end(); iter++ ) {
		cout << "Triplet "<<(*iter).first << " found " << (*iter).second << " times" << endl;
	}
	
	cout<<endl<<endl<<endl;
	for( map<nxsstring, int>::iterator iter = TripleCounts.begin(); iter != TripleCounts.end(); iter++ ) {
		cout << "Triple "<<(*iter).first << " found " << (*iter).second << " times" << endl;
	}
}
//for (int i=0;i<CladeVector.size();i++) {
//    cout<<endl;
//    for (int j=0;j<(CladeVector[i]).size();j++) {
//        cout<<" "<<CladeVector[i][j];
//   }
//}
//cout<<endl;
gsl_matrix_free(dmatrix);
return NJresult;
}


vector<nxsstring> BROWNIE::ReturnClade(Tree *T, nxsstring a, nxsstring b, nxsstring c) {
	NodePtr p;
	NodePtr q;
	std::map<Node *, int, std::less<Node *> > depth;
	PreorderIterator <Node> n ((*T).GetRoot());
	int count = 0;
	q = n.begin();
	while (q)
	{
		depth[q] = count++;
		q = n.next();
	}
	
	//Now get depths for each pair
	vector<int> LCADepthVectorT1;
	(*T).MakeNodeList();
	p = (*T).GetLeafWithLabel(a.c_str());
	q = (*T).GetLeafWithLabel(b.c_str());
	while (depth[p] != depth[q])
	{
		if (depth[p] < depth[q])
			q = q->GetAnc();
		else
			p = p->GetAnc();
	}
	LCADepthVectorT1.push_back(depth[p]);
	p = (*T).GetLeafWithLabel(b.c_str());
	q = (*T).GetLeafWithLabel(c.c_str());
	while (depth[p] != depth[q])
	{
		if (depth[p] < depth[q])
			q = q->GetAnc();
		else
			p = p->GetAnc();
	}
	LCADepthVectorT1.push_back(depth[p]);
	p = (*T).GetLeafWithLabel(a.c_str());
	q = (*T).GetLeafWithLabel(c.c_str());
	while (depth[p] != depth[q])
	{
		if (depth[p] < depth[q])
			q = q->GetAnc();
		else
			p = p->GetAnc();
	}
	LCADepthVectorT1.push_back(depth[p]);
	int T1DepthMax=0; //root has depth 0, others have higher depths.
	int T1DepthMaxIndex=-1;
	//int anum=taxa->FindTaxon(a);
	//int bnum=taxa->FindTaxon(b);
	//int cnum=taxa->FindTaxon(c);
	//cout<<anum<<": "<<a<<" "<<bnum<<": "<<b<<" "<<cnum<<": "<<c<<endl;
	for (int j=0;j<3;j++) {
		if(LCADepthVectorT1[j]>T1DepthMax) {
			T1DepthMax=LCADepthVectorT1[j];
			T1DepthMaxIndex=j;
		}
	}
	if (LCADepthVectorT1[0]==LCADepthVectorT1[1] && LCADepthVectorT1[0]==LCADepthVectorT1[2]) {
		T1DepthMaxIndex=-1;
	}
	//j=0->ab, j=1->bc j=2->ac
	//gsl_matrix_set(TaxonDistance,anum,anum,1+gsl_matrix_get(TaxonDistance,anum,anum));
	//gsl_matrix_set(TaxonDistance,bnum,bnum,1+gsl_matrix_get(TaxonDistance,bnum,bnum));
	//gsl_matrix_set(TaxonDistance,cnum,cnum,1+gsl_matrix_get(TaxonDistance,cnum,cnum));
	vector<nxsstring> returnvector;
	if (T1DepthMaxIndex==0) {
		returnvector.push_back(a);
		returnvector.push_back(b);
	}
	else if (T1DepthMaxIndex==1) {
		returnvector.push_back(b);
		returnvector.push_back(c);
	}
	else if (T1DepthMaxIndex==2) {
		returnvector.push_back(a);
		returnvector.push_back(c);
	}
	else if (T1DepthMaxIndex==-1) {
		returnvector.push_back("nope");
	}
	return returnvector;
}


void BROWNIE::InitializeQuartetCounts() {
	ProgressBar(trees->GetNumTrees());
	for (int treenum=0;treenum<trees->GetNumTrees(); treenum++) {
		Tree T=intrees.GetIthTree(treenum);
		vector<nxsstring> LeafLabelVector;
		NodeIterator <Node> n (T.GetRoot());
		NodePtr currentnode = n.begin();
		while (currentnode)
		{
			if(currentnode->IsLeaf()) {
				LeafLabelVector.push_back(currentnode->GetLabel());
			}
			currentnode = n.next();
		}
		gsl_combination *comb;
		comb=gsl_combination_calloc(T.GetNumLeaves(),4);
		do {
			nxsstring a=LeafLabelVector[gsl_combination_get(comb,0)];
			nxsstring b=LeafLabelVector[gsl_combination_get(comb,1)];
			nxsstring c=LeafLabelVector[gsl_combination_get(comb,2)];
			nxsstring d=LeafLabelVector[gsl_combination_get(comb,3)];
			map< nxsstring, int> LabelToNum;
			map< int, nxsstring > NumToLabel;
			vector<int> LabelNums;
			vector<int> OutputVector;
			LabelToNum[a]=taxa->FindTaxon(a);
			NumToLabel[LabelToNum[a]]=a;
			LabelToNum[b]=taxa->FindTaxon(b);
			NumToLabel[LabelToNum[b]]=b;
			LabelToNum[c]=taxa->FindTaxon(c);
			NumToLabel[LabelToNum[c]]=c;
			LabelToNum[d]=taxa->FindTaxon(d);
			NumToLabel[LabelToNum[d]]=d;
			int minlabel=-1;
			int maxlabel=100000000;
			map<nxsstring, int>::const_iterator itr;
			
			for(itr = LabelToNum.begin(); itr != LabelToNum.end(); ++itr){
				minlabel=GSL_MIN((*itr).second,minlabel);
				maxlabel=GSL_MAX((*itr).second,maxlabel);
				LabelNums.push_back((*itr).second);
			}
			
			vector<nxsstring> taxainclade;
			bool validquartet=false;
			taxainclade=ReturnClade(&T,LeafLabelVector[gsl_combination_get(comb,0)],LeafLabelVector[gsl_combination_get(comb,1)],LeafLabelVector[gsl_combination_get(comb,2)]);
			if(taxainclade.size()==2) {
				taxainclade=ReturnClade(&T,taxainclade[0],taxainclade[1],LeafLabelVector[gsl_combination_get(comb,3)]);	
			}
			if (taxainclade.size()==2) {
				OutputVector.push_back(minlabel);
				if (LabelToNum[taxainclade[0]]==minlabel) {
					
					OutputVector.push_back(LabelToNum[taxainclade[1]]);
					vector<int>remainingoutput;
					for (int i=0;i<4;i++) {
						if(LabelNums[i]!=OutputVector[0] && LabelNums[i]!=OutputVector[1]) {
							remainingoutput.push_back(LabelNums[i]);
						}
					}
					OutputVector.push_back(GSL_MIN(remainingoutput[0],remainingoutput[1]));
					OutputVector.push_back(GSL_MAX(remainingoutput[0],remainingoutput[1]));
				}
				else if (LabelToNum[taxainclade[1]]==minlabel) {
					OutputVector.push_back(LabelToNum[taxainclade[0]]);
					vector<int> remainingoutput;
					for (int i=0;i<4;i++) {
						if(LabelNums[i]!=OutputVector[0] && LabelNums[i]!=OutputVector[1]) {
							remainingoutput.push_back(LabelNums[i]);
						}
					}
					OutputVector.push_back(GSL_MIN(remainingoutput[0],remainingoutput[1]));
					OutputVector.push_back(GSL_MAX(remainingoutput[0],remainingoutput[1]));
				}
				else {
					vector<int> remainingoutput;
					for (int i=0;i<4;i++) {
						if(LabelNums[i]!=OutputVector[0] && LabelNums[i]!=OutputVector[1] && LabelNums[i]!=minlabel) {
							OutputVector.push_back(LabelNums[i]);
						}
					}
					OutputVector.push_back(GSL_MIN(LabelToNum[taxainclade[0]],LabelToNum[taxainclade[1]]));
					OutputVector.push_back(GSL_MAX(LabelToNum[taxainclade[0]],LabelToNum[taxainclade[1]]));
				}
				
			}
			else {
				OutputVector.push_back(maxlabel);
				vector<int> remainingoutput;
				int newmaxlabel=-1;
				for (int i=0;i<4;i++) {
					if(LabelNums[i]!=maxlabel) {
						newmaxlabel=GSL_MAX(LabelNums[i],newmaxlabel);
					}
				}
				OutputVector.push_back(newmaxlabel);
				int newmaxlabel2=-1;
				int newminlabel2=100000000;
				for (int i=0;i<4;i++) {
					if(LabelNums[i]!=maxlabel && LabelNums[i]!=newmaxlabel) {
						newmaxlabel2=GSL_MAX(LabelNums[i],newmaxlabel2);
						newminlabel2=GSL_MAX(LabelNums[i],newminlabel2);
					}
				}
				OutputVector.push_back(newmaxlabel2);
				OutputVector.push_back(newminlabel2);
			}
			quartetcounts[OutputVector]++;
		}
		while (gsl_combination_next(comb) == GSL_SUCCESS);
		gsl_combination_free(comb);
		ProgressBar(0);
	}
	
}

void BROWNIE::HandleNast( NexusToken& token ) {
	bool finishexecuting=true;
	for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (finishexecuting) {
                DoNast();
            }
            break;
        }
        else if( token.Abbreviation("?") ) {
			finishexecuting=false;
		}
        else if(token.Abbreviation("NPercent")) {
            nxsstring numbernexus;
			double initialnumbernexus;
            numbernexus = GetNumber(token);
            initialnumbernexus=atof( numbernexus.c_str() );
			if (initialnumbernexus>1 || initialnumbernexus<0) {
				message="Error: npercent should be between 0 and 1, inclusive";
				PrintMessage();
				finishexecuting=false;
			}
			else {
				npercent=initialnumbernexus;
			}
        }
	}
	
}

void BROWNIE::DoNast() {
	if (quartetcounts.size()==0) {
		message="Now initializing quartet counts";
		PrintMessage();
		InitializeQuartetCounts();
	}
	
}



vector<int> BROWNIE::GetTripletOverlap(ContainingTree *t1, ContainingTree *t2,int taxaincommon) {
    int maxnumber=0;
    int numberdisagree=0;
    int numunresolved=0; //number unresolved in one or the other trees
	int numunresolvedinboth=0; //number unresolved in both trees
	int numunresolvedinT1only=0;
	int numunresolvedinT2only=0;
	(*t1).MakeNodeList(); //added here, took out from containingtree:GetLCADepthVectorCommand
	(*t2).MakeNodeList();
    gsl_combination *c;
    c=gsl_combination_calloc(taxaincommon,3); //this is for triplets, change 3 to 4 (and make some other mods below) for quartets
    vector<nxsstring> LeafLabelVect=(*t1).GetLeafLabelVector();
    do
    {
        maxnumber++;
        vector<int> LCADepthVectorT1=(*t1).GetLCADepthVector(LeafLabelVect[gsl_combination_get(c,0)],LeafLabelVect[gsl_combination_get(c,1)],LeafLabelVect[gsl_combination_get(c,2)]);
        vector<int> LCADepthVectorT2=(*t2).GetLCADepthVector(LeafLabelVect[gsl_combination_get(c,0)],LeafLabelVect[gsl_combination_get(c,1)],LeafLabelVect[gsl_combination_get(c,2)]);
        int T1DepthMax=0; //root has depth 0, others have higher depths.
        int T2DepthMax=0;
        int T1DepthMaxIndex=0;
        int T2DepthMaxIndex=0;
        if ((LCADepthVectorT1[0]==LCADepthVectorT1[1] && LCADepthVectorT1[0]==LCADepthVectorT1[2]) || (LCADepthVectorT2[0]==LCADepthVectorT2[1] && LCADepthVectorT2[0]==LCADepthVectorT2[2])) {
            numunresolved++;
			if ((LCADepthVectorT1[0]==LCADepthVectorT1[1] && LCADepthVectorT1[0]==LCADepthVectorT1[2]) && (LCADepthVectorT2[0]==LCADepthVectorT2[1] && LCADepthVectorT2[0]==LCADepthVectorT2[2])) {
				numunresolvedinboth++;
			}
			else if ((LCADepthVectorT1[0]==LCADepthVectorT1[1] && LCADepthVectorT1[0]==LCADepthVectorT1[2])) {
				numunresolvedinT1only++;
			}
			else if ((LCADepthVectorT2[0]==LCADepthVectorT2[1] && LCADepthVectorT2[0]==LCADepthVectorT2[2])) {
				numunresolvedinT2only++;
			}
        }
        else {
            for (int j=0;j<3;j++) {
                if(LCADepthVectorT1[j]>T1DepthMax) {
                    T1DepthMax=LCADepthVectorT1[j];
                    T1DepthMaxIndex=j;
                }
                if(LCADepthVectorT2[j]>T2DepthMax) {
                    T2DepthMax=LCADepthVectorT2[j];
                    T2DepthMaxIndex=j;
                }
            }
            if (T2DepthMaxIndex!=T1DepthMaxIndex) {
                numberdisagree++;
            }
        }
    }
    while (gsl_combination_next(c) == GSL_SUCCESS);
    gsl_combination_free(c);
    vector<int> tripletoverlapoutput;
    tripletoverlapoutput.push_back(maxnumber);
    tripletoverlapoutput.push_back(numberdisagree);
    tripletoverlapoutput.push_back(numunresolved);
	tripletoverlapoutput.push_back(numunresolvedinboth);
	tripletoverlapoutput.push_back(numunresolvedinT1only);
	tripletoverlapoutput.push_back(numunresolvedinT2only);
    return tripletoverlapoutput;
}


void BROWNIE::HandleAccuracy( NexusToken& token )
{
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
			ComputeAccuracy();
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="Usage: Accuracy\n\n";
            message+="Outputs accuracy stats, assuming the first tree is the true tree.\n";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading ACCURACY command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}


void BROWNIE::ComputeAccuracy() {
	Tree TrueTreeGeneTreeFmt=intrees.GetIthTree(0); //we assume the first tree in the file is the accurate one
	ContainingTree TrueTree;
	TrueTree.SetRoot(TrueTreeGeneTreeFmt.CopyOfSubtree(TrueTreeGeneTreeFmt.GetRoot()));
	//message="TreeNum\tTreeName\tNumTripletsProperlyResolved\tNumTripletsProperlyUnresolved\tNumTripletsImproperlyResolved\tNumTripletsImproperlyUnresolved";
	message="\t\tCorrect\t\tIncorrect\nNum\tName\tRes\tUn\tRes\tUn";
	PrintMessage();
	for (int selectedtree=1;selectedtree<trees->GetNumTrees();selectedtree++) {
		ContainingTree ModifiedTrueTree=TrueTree;
		Tree CurrentGeneTreeTreeFmt=intrees.GetIthTree(selectedtree);
		ContainingTree CurrentGeneTree;
		CurrentGeneTree.SetRoot(CurrentGeneTreeTreeFmt.CopyOfSubtree(CurrentGeneTreeTreeFmt.GetRoot()));
		int ntaxincommon=PrepareTreesForTriplet(&ModifiedTrueTree,&CurrentGeneTree);
		vector<int> tripletoverlapoutput=GetTripletOverlap(&ModifiedTrueTree,&CurrentGeneTree,ntaxincommon);
		message="";
		message+=selectedtree+1;
		message+="\t";
		//if (selectedtree==0) {
		//	message+="TRUTH";
		//}
		//else {
			message+=trees->GetTreeName(selectedtree-1);
		//}
		message+="\t";
		
		double maxnumber=tripletoverlapoutput[0];
		double numberdisagree=tripletoverlapoutput[1];
		double numunresolved=tripletoverlapoutput[2];
		double numunresolvedinboth=tripletoverlapoutput[3];
		double numunresolvedinTTrueonly=tripletoverlapoutput[4];
		double numunresolvedinTOtheronly=tripletoverlapoutput[5];
		//message+=double((tripletoverlapoutput[0]-tripletoverlapoutput[1]-tripletoverlapoutput[4]-tripletoverlapoutput[5])/tripletoverlapoutput[0]);
		message+=(maxnumber-numunresolved)/maxnumber;
		message+="\t";
		//message+=double(tripletoverlapoutput[3]/tripletoverlapoutput[0]);
		message+=numunresolvedinboth/maxnumber;
		message+="\t";
		//message+=double(tripletoverlapoutput[4]/tripletoverlapoutput[0]);
		message+=numunresolvedinTTrueonly/maxnumber;
		message+="\t";
		//message+=double(tripletoverlapoutput[5]/tripletoverlapoutput[0]);
		message+=numunresolvedinTOtheronly/maxnumber;
		PrintMessage();
	}
}


double BROWNIE::ComputeTripletCost(int numberagree,int maxnumber,int ntaxincommon, double Tree1Wt,int Tree1Ntax,double Tree2Wt,int Tree2Ntax, int numberofgenes) {
    double initialscorealgorithm;
    double cdfundermodel;
    if (ntaxincommon<=50) { //change this to use exact values
                            //cdfundermodel=gsl_cdf_binomial_P(numberagree, 2.0/3.0, maxnumber);
        cdfundermodel=CDFvector[ntaxincommon-3][maxnumber-numberagree]; //assumes gene trees are resolved
                                                                        //	cout<<"ntax="<<ntaxincommon<<" cdfundermodelVector="<<cdfundermodel<<" ";
                                                                        //	if(ntaxincommon>5) {
                                                                        //		double x=(1.0*(maxnumber-numberagree))-(2.0*maxnumber/3.0);
                                                                        //		double calculatedvar=sqrt(((1.0/30.0)*double(gsl_sf_choose(ntaxincommon,5)))+((10.0/27.0)*double(gsl_sf_choose(ntaxincommon,4)))+((2.0/9.0)*double(gsl_sf_choose(ntaxincommon,3))));
                                                                        //		assert(calculatedvar>0);
                                                                        //		cout<<"cdfundermodelApprox="<<gsl_cdf_gaussian_P(x,calculatedvar)<<" ";
                                                                        //	}
    }
    else { //use approximation from D.E. Critchlow, D.K. Pearl, and C. Qian 1996. "The triples distance for rooted bifurcating phylogenetic trees" Syst Biol 45(3): 323-334.
        double x=(1.0*(maxnumber-numberagree))-(2.0*maxnumber/3.0);
        double calculatedvar=sqrt(((1.0/30.0)*double(gsl_sf_choose(ntaxincommon,5)))+((10.0/27.0)*double(gsl_sf_choose(ntaxincommon,4)))+((2.0/9.0)*double(gsl_sf_choose(ntaxincommon,3))));
        assert(calculatedvar>0);
        cdfundermodel=gsl_cdf_gaussian_P(x,calculatedvar);
        //	cout<<"ntax="<<ntaxincommon<<" cdfundermodelApprox="<<cdfundermodel<<" ";
        //	for (int i=0;i<=35;i++) {
        //		 x=(1.0*i)-(2.0*35/3.0);
        //		 calculatedvar=sqrt(((1.0/30.0)*double(gsl_sf_choose(7,5)))+((10.0/27.0)*double(gsl_sf_choose(7,4)))+((2.0/9.0)*double(gsl_sf_choose(7,3))));
        //		cout<<"i="<<i<<" cdf="<<gsl_cdf_gaussian_P(x,calculatedvar)<<endl;
        //	}
    }
    if (cdfundermodel>pthreshold) {
        initialscorealgorithm=0; //so we don't penalize random structure
    }
    else {
        initialscorealgorithm=Tree1Wt*Tree2Wt*((1.0/cdfundermodel)-1.0)/(1.0*numberofgenes); //Motivation: treeweighting obvious; use number of identical triplets as score; divide by number of genes since we do O(numberofgenes*numberofgenes) comparisons for structure but only O(numberofgenes) costs for GTP and want the costs to be similar regardless of gene number
		if (debugmode) {
			if (isinf(initialscorealgorithm)==1) {
				cout<<"Error: initialscorealgorithm = "<<initialscorealgorithm<<" cdfundermodel = "<<cdfundermodel;
				if (ntaxincommon<=32) {
					cout<<" calculated using cdfundermodelVector"<<endl;
				}
				else {
					cout<<" calculated using Critchlow et al approximation"<<endl;
				}
			}
		}
    }
                                                                                         //initialscorealgorithm=Tree1Wt*Tree2Wt*(maxnumber-numberdisagree)/numberofgenes;

    return initialscorealgorithm;
}

//Gets trees ready to compare using triplet or quartet scores. After this, just use
//QTValues Q;
//CompareTriplets(t1,t2, Q);
//SummaryStats(Q);
//and get the results you want from Q.
int BROWNIE::PrepareTreesForTriplet(ContainingTree *t1, ContainingTree *t2) {
    int taxaincommon=0;
    (*t1).Update();
    (*t2).Update();
    (*t1).SetLeafNumbers();
    (*t2).SetLeafNumbers();
    //cout<<"starting to compare "<<endl;
    //(*t1).Write(cout);
    //cout<<endl;
    //(*t2).Write(cout);
    //cout<<endl;
    taxaincommon=PruneToOverlappingLeaves(t1, t2);
    if (taxaincommon>2) {
        //   (*t1).BuildLeafClusters();
        //   (*t1).BuildLabelClusters();
        //   (*t1).MakeNodeList();
        //   (*t2).BuildLeafClusters();
        //   (*t2).BuildLabelClusters();
        //   (*t2).MakeNodeList();
        (*t1).SetLeafNumbers();
        (*t2).SetLeafNumbers();
		assert((*t1).GetNumLeaves()==(*t2).GetNumLeaves());

    }
    //cout<<"continuing to compare "<<endl;
    //(*t1).Write(cout);
    //cout<<endl;
    //(*t2).Write(cout);
    //cout<<endl<<taxaincommon<<endl<<endl;
    return taxaincommon;
}

//Uses some code modified from Rod Page's supertrees program to prune leaves from two trees so they have the same leaf set
int BROWNIE::PruneToOverlappingLeaves(ContainingTree *t1, ContainingTree *t2)
{
	//cout<<"\n\n\nStarting Prune\n\n";
	//(*t1).ReportTreeHealth();
	//(*t2).ReportTreeHealth();
    (*t1).MakeNodeList();
    (*t2).MakeNodeList();
    int maxoverlap=(*t1).GetNumLeaves();
    int leavesincommon=0;
    vector<NodePtr> NodesToDelete;
    NodeIterator <Node> n ((*t1).GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if (currentnode->IsLeaf()) {
            NodePtr matchingLeaf = (*t2).GetLeafWithLabel(currentnode->GetLabel());
            //cout<<"Current node is "<<currentnode<<" "<<currentnode->GetLabel()<<" and matching leaf is "<<matchingLeaf<<endl;
            if (matchingLeaf != NULL)
            {
                // This leaf is in t1 AND t2. Set the LeafNumber of this leaf
                // in t1 to match that in t2
                leavesincommon++;
                currentnode->SetLeafNumber(leavesincommon); //We need to do this & the following step so with MakeNodeList, the leaf number is 1:#Leaves
                matchingLeaf->SetLeafNumber(leavesincommon);

            }
            else
            {
                // This leaf is not in t2 so we will prune it from the tree
                NodesToDelete.push_back(currentnode);
                //cout<<"Going to delete "<<currentnode<<" with label "<<currentnode->GetLabel()<<" from t1"<<endl;
                maxoverlap--;

            }
        }
        currentnode=n.next();

    }
    if (maxoverlap>0) {
        for (int i=0;i<NodesToDelete.size();i++) {
            (*t1).DeleteSubtendingEdge(NodesToDelete[i]);
            delete (NodesToDelete[i]);
        }
        NodesToDelete.clear();
        (*t1).SuppressInternalNodesWithOneDescendant();
        (*t1).Update();
        (*t1).GetNodeDepths();
		(*t1).SetLeafNumbers();
        (*t2).SetLeafNumbers();
        //now do the same for the reverse comparison

        (*t1).MakeNodeList();
        (*t2).MakeNodeList();
		//(*t1).ReportTreeHealth();
		//(*t2).ReportTreeHealth();
        leavesincommon=0;
        NodeIterator <Node> m ((*t2).GetRoot());
        currentnode = m.begin();
        while (currentnode)
        {
            if (currentnode->IsLeaf()) {
                NodePtr matchingLeaf = (*t1).GetLeafWithLabel (currentnode->GetLabel());
                if (matchingLeaf != NULL)
                {
                    leavesincommon++;
                    currentnode->SetLeafNumber(leavesincommon); //We need to do this & the following step so with MakeNodeList, the leaf number is 1:#Leaves
                    matchingLeaf->SetLeafNumber(leavesincommon);
                }
                else
                {
                    NodesToDelete.push_back(currentnode);
                    //cout<<"Going to delete "<<currentnode<<" with label "<<currentnode->GetLabel()<<" from t2"<<endl;
                }
            }
            currentnode=m.next();
        }
        for (int i=0;i<NodesToDelete.size();i++) {
			//cout<<"\n------------------- \nGoing to delete "<<NodesToDelete[i]<<endl<<endl;
            (*t2).DeleteSubtendingEdge(NodesToDelete[i]);
            delete (NodesToDelete[i]);
			
			//(*t2).ReportTreeHealth();
			//cout<<"\nJust deleted "<<NodesToDelete[i]<<endl<<endl;
        }
        if (leavesincommon>2) {
            (*t2).SuppressInternalNodesWithOneDescendant();
            (*t2).Update();
            (*t2).GetNodeDepths();
        }
    }
			//(*t1).ReportTreeHealth();
			//(*t2).ReportTreeHealth();
    return (GSL_MIN(leavesincommon,maxoverlap));
}



//Eliminate duplicate trees in the FormattedBestTrees vector
void BROWNIE::DelDupes()
{
    vector<ContainingTree> FormattedBestTreesBkup;
    FormattedBestTreesBkup.swap(FormattedBestTrees); //empties besttrees into best treesbkup
    assert(FormattedBestTreesBkup.size()>0);
    FormattedBestTrees.push_back(FormattedBestTreesBkup.back());
    (FormattedBestTrees.back()).FindAndSetRoot();
    (FormattedBestTrees.back()).Update();
    FormattedBestTreesBkup.pop_back(); // transfer a tree back into best trees bkup
    for (int i=0; i<FormattedBestTreesBkup.size(); i++) {
        (FormattedBestTreesBkup[i]).FindAndSetRoot();
        (FormattedBestTreesBkup[i]).Update();
        bool newtree=true;
        for (int j=0; j<FormattedBestTrees.size(); j++) {
            ContainingTree t1=FormattedBestTreesBkup[i];
            ContainingTree t2=FormattedBestTrees[j];
            // cout<<"Compare \n\t";
            // (FormattedBestTreesBkup[i]).Write(cout);
            //cout<<"\n\t";
            //  (FormattedBestTrees[j]).Write(cout);
            int ntaxt1=t1.GetNumLeaves();
            int ntaxt2=t2.GetNumLeaves();
            int ntaxincommon=PrepareTreesForTriplet(&t1,&t2);
            // cout<<"\n\tT1 leaves before = "<<ntaxt1<<" after "<<t1.GetNumLeaves()<<"\n\tT2 leaves before = "<<ntaxt2<<" after "<<t2.GetNumLeaves()<<"\n\tntaxincommon = "<<ntaxincommon<<endl;
            if (ntaxincommon==ntaxt1 && ntaxincommon==ntaxt2) {
                if (ntaxincommon>=3) { //so, all leaves are in common, but is the topology the same? If there are only two or one species, yes, so newtree=false. If there are at least 3 taxa, there are different possible topologies, so maybe, thus check below
                    vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,ntaxincommon);
                    int numberdisagree=tripletoverlapoutput[1];
                    if (numberdisagree==0) {
                        newtree=false;
                    }
                    //      cout<<"\tNumDisagree = "<<numberdisagree<<endl;
                }
                else {
                    newtree=false; //all leaves in common, must have same topology, so same tree
                }
            }
            //ShowHeader(cout);
            //ShowQTRecord(cout,Q);
            //cout<<endl<<endl;
        }
        if (newtree) {
            FormattedBestTrees.push_back(FormattedBestTreesBkup[i]); //BestTreesBkup is going out of scope soon, anyway, so don't need to bother pop_backing
        }
    }
    for (int k=0; k<FormattedBestTrees.size(); k++) {
        (FormattedBestTrees[k]).Update();
        (FormattedBestTrees[k]).GetNodeDepths();
        (FormattedBestTrees[k]).ClearInternalLabels();
    }
}

/**
* @method HandleCitation [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the CITATION command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleCitation( NexusToken& token )
{
    nxsstring numbernexus;
    bool donenothing=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (donenothing) {
				PrintCitations();
            }
            break;
        }
        else if( token.Abbreviation("?") ) {
            donenothing=false;
            message="Usage: Citation\n\n";
            message+="Prints out list of papers that you should read and cite for the analyses\nrun this session.";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading CHOOSE command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

void BROWNIE::PrintCitations()
{
	message="Papers whose methods you have used so far in this session:\n[you should read and probably cite them]";
	if (citationarray[0]) {
		message+="\n\nO'Meara, B.C., C. Ane, M.J. Sanderson, and P.C. Wainwright. 2006. \"Testing for different rates of evolution using likelihood.\" Evolution 60(5): 922-933\n---Citation for this program and for rate comparison methods.";
	}
	if (citationarray[1]) {
		message+="\n\nBlomberg, S.P., T. Garland, Jr., and A.R. Ives. 2003. \"Testing for phylogenetic signal in comparative data: Behavioral traits are more labile.\" Evolution 57(4) 717-745.\n---Citation for constant mean, constant pull OU and ACDC transformations (d and g parameters, respectively).";
	}
				if (citationarray[2]) {
					message+="\n\nButler, M.A., King, A.A. 2004. \"Phylogenetic comparative analysis: a modeling approach for adaptive evolution.\" American Naturalist. 164(6):683-695.";
					message+="\n\nHansen, T.F., 1997. \"Stabilizing selection and the comparative analysis of adaptation.\" Evolution, 51:1341-1351.";
					message+="\n\nO'Meara, B.C. Brownie v2.0b7. Distributed by the author at http://www.brianomeara.info/brownie";
					message+="\n---Citations for Ornstein-Uhlenbeck model with multiple means but one attraction and rate parameter";
				}
	if (citationarray[3]) {
			message+="\n\nO'Meara, B.C. MS in prep \"Species delimitation using multiple gene trees\"";
			message+="\n---Citation for species delimitation approach";
	}
	PrintMessage();
}

/**
* @method HandleChoose [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the CHOOSE command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleChoose( NexusToken& token )
{
    nxsstring numbernexus;
    bool donenothing=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (donenothing) {
                message="Usage: Choose [tree=<integer>] [char=<integer>]\n\n";
                PrintMessage();
            }
            break;
        }
        else if( token.Abbreviation("Tree") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            chosentree=atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen tree number ";
            message+=chosentree;
            PrintMessage();
            if (chosentree<1) {
                errormsg = "Error: must select a number greater than zero";
                chosentree=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            if (chosentree>trees->GetNumTrees()) {
                errormsg = "Error: you chose tree number ";
                errormsg += chosentree;
                errormsg += " but there are only ";
                errormsg += trees->GetNumTrees();
                errormsg += " trees loaded.\n";
                errormsg += "Tree 1 has been selected by default.";
                chosentree=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("Char") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            chosenchar=atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen char number ";
            message+=chosenchar;
            PrintMessage();
            if (chosenchar<1) {
                errormsg = "Error: must select a number greater than zero";
                chosenchar=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            if (chosenchar>characters->GetNChar()) {
                errormsg = "Error: you chose char number ";
                errormsg += chosenchar;
                errormsg += " but there are only ";
                errormsg += characters->GetNChar();
                errormsg += " characters loaded.\n";
                errormsg += "Character 1 has been selected by default.";
                chosenchar=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("?") ) {
            donenothing=false;
            message="Usage: Choose [tree=<integer>] [char=<integer>]\n\n";
            message+="Selects one tree and/or character to use for subsequent analyses.\n";
            message+="The remaining trees and characters are still stored in memory and can be chosen later.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --\n";
            message+="Tree         <integer-value>                      ";
            if (trees->GetNumTrees()>0) {
                message+=chosentree;
            }
            else {
                message+="[no trees loaded]";
            }
            message+="\nChar         <integer-value>                      ";
            if (characters->GetNChar()>0) {
                message+=chosenchar;
            }
            else {
                message+="[no characters loaded]";
            }
            message+="\n\n";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading CHOOSE command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

void BROWNIE::HandleCompareRandomTrees( NexusToken& token)
{
    bool donenothing=false;
    int ltax=3;
    int rtax=3;
    int nreps=1000;
    bool loopmode=false;
    bool automode=false;
    int requiredminimum=5;
    nxsstring numbernexus;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (donenothing) {
                message="Usage: Compare [ltax=<integer> rtax=<integer> nreps=<integer> [loop]]\n\nIf you use loop, ltax=taxmax and rtax=taxmin";
                PrintMessage();
            }
            else {
                if (loopmode && !automode) {
                    message="loopmode=yes automode=no\nltax\trtax\tnreps\tdiffs\tobs\tprobability\tcum prob";
                    PrintMessage();
                    int taxmax=ltax;
					int taxmin=rtax;
					nxsstring filename="comparisoncodeMax";
					filename+=taxmax;
					filename+="Min";
					filename+=taxmin;
					filename+=".cpp";
					ofstream comparisonfile;
					comparisonfile.open( filename.c_str() );
                    for (int ntax=taxmin;ntax<=taxmax;ntax++) {
                        ltax=ntax;
                        rtax=ntax;
                        int minnumleaves=GSL_MIN(ltax,rtax);
                        int maxnumsametriplets=(minnumleaves*(minnumleaves-1)*(minnumleaves-2))/6; //minnumleaves choose 3
						vector<int> differencevector(maxnumsametriplets+1,0); // the +1 is so we can have zero in the vector, too
						for (int i=0;i<nreps;i++) {
							ContainingTree t1;
							ContainingTree t2;
							t1.RandomTree(ltax);
							t2.RandomTree(rtax);
							t1.ConvertTaxonNamesToRandomTaxonNumbers();
							t2.ConvertTaxonNamesToRandomTaxonNumbers();
							int taxaincommon=PrepareTreesForTriplet(&t1,&t2);
							// QTValues Q;
							// CompareTriplets(t1,t2, Q);
							// SummaryStats(Q);
							// int numberdisagree=Q.d;
							vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,taxaincommon);
							int numberdisagree=tripletoverlapoutput[1];
							if (debugmode) {
								cout<<"rep "<<i+1<<" "<<numberdisagree<<endl;
							}
							differencevector[numberdisagree]++;
						}
						
						double cumprob=0;
						comparisonfile<<"contentsofrow.clear();\n";
						for (int i=0; i<differencevector.size(); i++) {
							message="";
							message+=ltax;
							message+="\t";
							message+=rtax;
							message+="\t";
							message+=nreps;
							message+="\t";
							message+=i;
							message+="\t";
							message+=differencevector[i];
							message+="\t";
							char outputstring[11];
							double ratio=(1.0*differencevector[i])/nreps;
							cumprob+=ratio;
							sprintf(outputstring,"%11.10f",ratio);
							message+=outputstring;
							message+="\t";
							sprintf(outputstring,"%11.10f",cumprob);
							message+=outputstring;
							PrintMessage();
							comparisonfile<<"contentsofrow.push_back("<<cumprob<<");\n";
						}
						comparisonfile<<"CDFvector.push_back(contentsofrow); //comparison where ntax="<<ntax<<endl;
                    }
					comparisonfile.close();
                }
                else if (automode) {
                    int minnumleaves=GSL_MIN(ltax,rtax);
                    int maxnumsametriplets=(minnumleaves*(minnumleaves-1)*(minnumleaves-2))/6; //minnumleaves choose 3
					vector<int> differencevector(maxnumsametriplets+1,0); // the +1 is so we can have zero in the vector, too
					int vectorminimum=0;
					int loopcounttotal=0;
					int loopcount=0;
					message="Loops\tNumber of samples in the 4 most similar & 4 least similar scores\tTime remaining";
					PrintMessage();
					nxsstring filename="comparisonof";
					filename+=GSL_MIN(ltax,rtax);
					filename+=".txt";
					ofstream comparisonfile;
					comparisonfile.open( filename.c_str() );
					time_t start,end;
					time (&start);
					int oldminimum=0;
					double secondstotal=-1;
					while (vectorminimum<requiredminimum) {
						ContainingTree t1;
						ContainingTree t2;
						t1.RandomTree(ltax);
						t2.RandomTree(rtax);
						t1.ConvertTaxonNamesToRandomTaxonNumbers();
						t2.ConvertTaxonNamesToRandomTaxonNumbers();
						int taxaincommon=PrepareTreesForTriplet(&t1,&t2);
						QTValues Q;
						CompareTriplets(t1,t2, Q);
						SummaryStats(Q);
						int numberdisagree=Q.d;
						comparisonfile<<1.0*numberdisagree/maxnumsametriplets<<endl;
						differencevector[numberdisagree]++;
						vectorminimum=requiredminimum*10;
						for (int i=0; i<differencevector.size(); i++) {
							if (differencevector[i]<vectorminimum) {
								vectorminimum=differencevector[i];
							}
						}
						if(vectorminimum>oldminimum) {
							oldminimum=vectorminimum;
							time(&end);
							double timedifference=difftime(end,start);
							secondstotal=(1.0*timedifference)*(1.0*requiredminimum)/(1.0*vectorminimum); // all the 1.0 to make them treated as doubles.
																										 // cout<<"timedif is "<<timedifference<<" sectotal is "<<secondstotal<<" vectorminimum is "<<vectorminimum<<" requiredminimum is "<<requiredminimum<<endl;
						}
						loopcount++;
						loopcounttotal++;
						if (loopcount==5000) {
							message="";
							message+=loopcounttotal;
							message+="\t";
							for (int k=0;k<GSL_MIN(4,differencevector.size());k++) {
								message+="\t";
								message+=differencevector[k];
							}
							message+="\t...";
							for (int k=GSL_MAX(0,differencevector.size()-4);k<differencevector.size();k++) {
								message+="\t";
								message+=differencevector[k];
							}
							message+="\t";
							time(&end);
							if (secondstotal>3600) {
								int hours=int(floor((secondstotal-difftime(end,start))/3600.0));
								int minutes=int(floor(((secondstotal-difftime(end,start))/60.0)-(hours*60)));
								message+=hours;
								message+="h ";
								message+=minutes;
								message+=" min";
							}
							else if (secondstotal>120) {
								int minutes=int(ceil((secondstotal-difftime(end,start))/60.0));
								if (minutes<1) {
									message+="almost done";
								}
								else {
									message+=minutes;
									message+=" min";
								}
							}
							else if (secondstotal>0) {
								int seconds=int(ceil((secondstotal-difftime(end,start))));
								if (seconds<1) {
									message+="almost done";
								}
								else {
									message+=seconds;
									message+=" sec";
								}
							}
							else {
								message+="Too long";
							}
							PrintMessage();
							loopcount=0;
						}
					}
					comparisonfile.close();
					message="";
					message+=loopcounttotal;
					message+="\t";
					for (int k=0;k<GSL_MIN(4,differencevector.size());k++) {
						message+="\t";
						message+=differencevector[k];
					}
					message+="\t...";
					for (int k=GSL_MAX(0,differencevector.size()-4);k<differencevector.size();k++) {
						message+="\t";
						message+=differencevector[k];
					}
					message+="\n\t\t\tElapsed time: ";
					time(&end);
					message+=floor(difftime(end,start)/60);
					message+=" min";
					message+="\n";
					PrintMessage();
					message="ltax\trtax\tnreps\tdiffs\tobs\tproportion\tprobability\tcum prob";
					PrintMessage();
					double cumprob=0;
					int maxcount=0;
					for (int i=0; i<differencevector.size(); i++) {
						if (differencevector[i]>maxcount) {
							maxcount=differencevector[i];
						}
					}
					for (int i=0; i<differencevector.size(); i++) {
						message="";
						message+=ltax;
						message+="\t";
						message+=rtax;
						message+="\t";
						message+=loopcounttotal;
						message+="\t";
						message+=i;
						message+="\t";
						message+=differencevector[i];
						message+="\t";
						char outputstring[11];
						double proportion=(1.0*i/maxnumsametriplets);
						sprintf(outputstring,"%11.10f",maxnumsametriplets);
						message+=outputstring;
						message+="\t";
						double ratio=(1.0*differencevector[i])/loopcounttotal;
						cumprob+=ratio;
						sprintf(outputstring,"%11.10f",ratio);
						message+=outputstring;
						message+="\t";
						sprintf(outputstring,"%11.10f",cumprob);
						message+=outputstring;
						message+="\t";
						for (int k=0;k<floor(10.0*(differencevector[i])/maxcount);k++) {
							message+="*";
						}
						PrintMessage();
					}
					
					gsl_vector * rawdata = gsl_vector_calloc(loopcounttotal);
					int insertionposition=0;
					for (int i=0; i<differencevector.size(); i++) {
						for (int j=0;j<differencevector[i];j++) {
							gsl_vector_set(rawdata,insertionposition,1.0*i/maxnumsametriplets);
							insertionposition++;
						}
					}
					
					message="\n";
					gsl_sort_vector(rawdata);
					int cdfresolution=100;
					for (int k=1;k<=cdfresolution;k++) {
						double proportiondesired=1.0*k/cdfresolution;
						message+="CDF percentile = ";
						message+=proportiondesired;
						message+=" proportion mismatch = ";
						message+=gsl_stats_quantile_from_sorted_data(rawdata->data,1,loopcounttotal,proportiondesired);
						message+="\n";
					}
					
					message+="\nmean is ";
					message+=gsl_stats_mean(rawdata->data, 1, loopcounttotal);
					message+="\nvariance is ";
					message+=gsl_stats_variance(rawdata->data, 1, loopcounttotal);
					PrintMessage();
                }
				
                else {
                    int minnumleaves=GSL_MIN(ltax,rtax);
                    int maxnumsametriplets=(minnumleaves*(minnumleaves-1)*(minnumleaves-2))/6; //minnumleaves choose 3
					vector<int> differencevector(maxnumsametriplets+1,0); // the +1 is so we can have zero in the vector, too
					ProgressBar(nreps);
					for (int i=0;i<nreps;i++) {
						ProgressBar(0);
						ContainingTree t1;
						ContainingTree t2;
						t1.RandomTree(ltax);
						t2.RandomTree(rtax);
						t1.ConvertTaxonNamesToRandomTaxonNumbers();
						t2.ConvertTaxonNamesToRandomTaxonNumbers();
						int taxaincommon=PrepareTreesForTriplet(&t1,&t2);
						QTValues Q;
						CompareTriplets(t1,t2, Q);
						SummaryStats(Q);
						int numberdisagree=Q.d;
						differencevector[numberdisagree]++;
					}
					message="ltax\trtax\tnreps\tdiffs\tobs\tprobability\tcum prob";
					PrintMessage();
					double cumprob=0;
					for (int i=0; i<differencevector.size(); i++) {
						message="";
						message+=ltax;
						message+="\t";
						message+=rtax;
						message+="\t";
						message+=nreps;
						message+="\t";
						message+=i;
						message+="\t";
						message+=differencevector[i];
						message+="\t";
						char outputstring[11];
						double ratio=(1.0*differencevector[i])/nreps;
						cumprob+=ratio;
						sprintf(outputstring,"%11.10f",ratio);
						message+=outputstring;
						message+="\t";
						sprintf(outputstring,"%11.10f",cumprob);
						message+=outputstring;
						PrintMessage();
					}
					nxsstring filename="comparisonof";
					filename+=GSL_MIN(ltax,rtax);
					filename+=".txt";
					ofstream comparisonfile;
					comparisonfile.open( filename.c_str() );
					gsl_vector * rawdata = gsl_vector_calloc(nreps);
					int insertionposition=0;
					for (int i=0; i<differencevector.size(); i++) {
						for (int j=0;j<differencevector[i];j++) {
							comparisonfile<<i<<endl;
							gsl_vector_set(rawdata,insertionposition,1.0*i/maxnumsametriplets);
							insertionposition++;
						}
					}
					comparisonfile.close();
					gsl_sort_vector(rawdata);
					for (int k=1;k<=10;k++) {
						double proportiondesired=k/10.0;
						cout<<"CDF percentile = "<<proportiondesired<<" proportion mismatch = "<<gsl_stats_quantile_from_sorted_data(rawdata->data,1,nreps,proportiondesired)<<endl;
					}
					cout<<"mean is "<<gsl_stats_mean(rawdata->data, 1, nreps)<<endl;
					cout<<"variance is "<<gsl_stats_variance(rawdata->data, 1, nreps)<<endl;
                }
            }
            break;
        }
        else if( token.Abbreviation("Ltax") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            ltax=atoi( numbernexus.c_str() ); //convert to int
            if (ltax<3) {
                errormsg = "Error: must select a number greater than 2";
                ltax=3;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("Rtax") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            rtax=atoi( numbernexus.c_str() ); //convert to int
            if (rtax<3) {
                errormsg = "Error: must select a number greater than 2";
                rtax=3;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("LOop") ) {
            loopmode=true;
        }
        else if ( token.Abbreviation("AUto") ) {
            automode=true;
        }
        else if( token.Abbreviation("Nreps") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            nreps=atoi( numbernexus.c_str() ); //convert to int
            if (nreps<1) {
                errormsg = "Error: must select a number greater than 0";
                nreps=1000;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
		else if( token.Abbreviation("Jreps") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            jreps=atoi( numbernexus.c_str() ); //convert to int
            if (jreps<1) {
                errormsg = "Error: must select a number greater than 0";
                jreps=0;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
		else if( token.Abbreviation("PCtdelete") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            pctdelete=atof( numbernexus.c_str() ); //convert to int
            if (pctdelete>=1 || pctdelete<=0) {
                errormsg = "Error: must select a number between 0 and 1 (and not equal to either of these)";
                pctdelete=1.0/3.0;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("?") ) {
            donenothing=false;
            message="Usage: Compare [ltax=<integer> rtax=<integer> nreps=<integer>]\n\n";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading COMPARE command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}



/**
* @method HandleSet [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the SET command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleSet( NexusToken& token )
{
    nxsstring numbernexus;
    bool donenothing=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (donenothing) {
                message="Usage: Set [maxspecies=<integer>]\n\n";
                PrintMessage();
            }
            break;
        }
        else if( token.Abbreviation("Maxspecies") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            maxnumspecies=atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen to limit the search to no more than ";
            message+=maxnumspecies;
            PrintMessage();
            if (maxnumspecies<1) {
                errormsg = "Error: must select a number greater than zero";
                maxnumspecies=taxa->GetNumTaxonLabels();
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            if (maxnumspecies>taxa->GetNumTaxonLabels()) {
                errormsg = "Error: you chose a maximum of ";
                errormsg += maxnumspecies;
                errormsg += " but there are only ";
                errormsg += taxa->GetNumTaxonLabels();
                errormsg += " samples loaded.\n";
                maxnumspecies=taxa->GetNumTaxonLabels();
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("Compare")) {
            ContainingTree t1;
            ContainingTree t2;
            t1.RandomTree(6);
            t2.RandomTree(4);
            t1.FindAndSetRoot();
            t1.Update();
            t1.GetNodeDepths();
            t1.ConvertTaxonNamesToRandomTaxonNumbers();
            t2.ConvertTaxonNamesToRandomTaxonNumbers();
            t1.SetLeafNumbers();
            t2.SetLeafNumbers();
            int taxaincommon=PruneToOverlappingLeaves(&t1,&t2);
            cout<<"Taxa in common are "<<taxaincommon<<endl;
            t1.BuildLeafClusters();
            t1.BuildLabelClusters();
            t1.MakeNodeList();
            t2.FindAndSetRoot();
            t2.Update();
            t2.GetNodeDepths();
            t2.BuildLeafClusters();
            t2.BuildLabelClusters();
            t2.MakeNodeList();

            cout<<"clusters tree 1\n";
            t1.ShowClusters();
            cout<<"clusters tree 2\n";
            t2.ShowClusters();
            cout<<"Now comparing\n";
            t1.Draw(cout);
            cout<<endl;
            t2.Draw(cout);
            cout<<endl;
            QTValues Q;
            CompareTriplets(t1,t2, Q);
            SummaryStats(Q);
            int numberdisagree=Q.d;
            int maxnumber=Q.n;
            cout<<"Out of "<<maxnumber<<" only "<<numberdisagree<<" disagreed\n";
        }
        else if( token.Abbreviation("Randtree")) {
            donenothing=false;
            ContainingTree testtree;
            testtree.RandomTree(6);
            cout<<endl;
            // testtree.Write(cout);
            testtree.ConvertTaxonNamesToOrderedTaxonNumbers();
            cout<<endl;
            //testtree.Write(cout);
            cout<<endl;
            testtree.Update();
            //testtree.Draw(cout);
            cout<<endl<<"Tree has "<<testtree.GetNumLeaves()<<" leaves and "<<testtree.GetNumInternals()<<" internals"<<endl;
            testtree.TestRerooting();
            // testtree.Update();
            //  testtree.Draw(cout);
            cout<<endl<<"Tree has "<<testtree.GetNumLeaves()<<" leaves and "<<testtree.GetNumInternals()<<" internals"<<endl;
        }
        else if( token.Abbreviation("Write")) {
            donenothing=false;
            ContainingTree testtree;
            int nspecies=0;
            for (int vectorpos=0;vectorpos<convertsamplestospecies.size();vectorpos++) {
                if (convertsamplestospecies[vectorpos]>nspecies) {
                    nspecies=convertsamplestospecies[vectorpos];
                }
            }
            testtree.RandomTree(nspecies);
            testtree.Draw(cout);
            testtree.ConvertTaxonNamesToRandomTaxonNumbers();
            testtree.Draw(cout);
            cout<<OutputForGTP(&testtree);
            //printUsage();
            //cout<<"ReturnScore = "<<ReturnScore(OutputForGTP(&testtree),unrooted);
        }
        else if( token.Abbreviation("?") ) {
            donenothing=false;
            message="Usage: Set [maxspecies=<integer>]\n\n";
            message+="Sets the maximum number of species to test.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --\n";
            message+="MaxSpecies   <integer-value>                      ";
            message+=maxnumspecies;
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading CHOOSE command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}



void BROWNIE::HandleModel( NexusToken& token )
{
	// Retrieve all tokens for this command, stopping only in the event
	// of a semicolon or an unrecognized keyword
	//
    for(;;)
    {
        token.GetNextToken();
		
        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
			//message="Usage: Model type= [BM | OU | ACDC] assign=[1 | CS | CC | TI | TC] ou=[M | MV]\n";
        //    message="Usage: Model type = [ BM1 | OU1 | ACDC1 | BMS | BMC | BMAN | BMAO | OUSMVA | OUCMVA | OUSMA | OUCMA | OUSM | OUCM | TSBMI |TSBMC]\n\n"; //add this once everything's been tested
            message="Usage: Model type = [ BM1 | BMS | BMC | OU1 | ... ] states = ( vector )\n";
            message+="";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------------------------------- Current setting -----";
            message+="\n Type        BM1 | BMS | BMC | OU1 | ACDC                                  ";

		//	message+="\n Type        BM | OU | ACDC                                                           ";
		//	message+="\n Assign      1 | CS | CC | TI | TC | RBS | RBO                                                    ";
		//	message+="\n OU          M | MV                                                                   ";
         //   message+="\n Type        BM1 | OU1 | ACDC1 | BM[S|C] | OU[SMVA|CMVA|SMA|CMA|CM|SM] | TS[BMS|BMC]  ";
            if (chosenmodel==1 || chosenmodel==2) {
                message+="BM1";
            }
            else if(chosenmodel==3) {
                message+="OU1";
            }
            else if(chosenmodel==4) {
                message+="ACDC";
            }
            else if(chosenmodel==5) {
                message+="BMS";
            }
            else if(chosenmodel==6) {
                message+="BMC";
            }
			else if(chosenmodel==20) {
				message+="Pagel's Kappa";
			}
			else if(chosenmodel==21) {
				message+="Pagel's Delta";
			}
			else if(chosenmodel==22) {
				message+="Pagel's Lambda";
			}
            else {
                message+=chosenmodel;
            }
			message+="\n             OUSM | OUCM";
            message+="\n States      vector                                    (";
            for (int i=0; i<staterestrictionvector.size();i++) {
                message+=" ";
                message+=staterestrictionvector[i];
            }
            message+=" )";
           // message+="\n\nType\nBM1    = Brownian motion, one rate parameter";
			//ROse suggests OUCMA approach would be better
			
           // message+="\nOU: Ornstein-Uhlenbeck";
		// message+="\nACDC: Acceleration/deceleration model of Blomberg et al (2003)";
           // message+="\nAssign=1: Assign one model over the whole tree";
          //  message+="\nAssign=CS: Assign one model for each discrete Character State mapped on the tree";
           // message+="\nAssign=CC: Assign one model on branches where the discrete Character Changes, and a different model to the other branches";
         //   message+="\nAssign=TI: Assign one model for each Time Interval";
         //   message+="\nAssign=TC: Assign one model on branches that Cross a Time interval, and a different model to the other branches";
         //   message+="\nOU=M: Allow means (optimal values) to vary across OU models, but keep variance and attraction parameters the same";
          //  message+="\nOU-MV: Allow different means and variances across OU models, but keep attraction parameter the same";
			message+="\n\nBM1     = Brownian motion, one rate parameter";
			message+="\nBMS     = Brownian motion, with different rate parameters for each state on a tree";
			message+="\nBMC     = Brownian motion, with one rate parameter for branches with state changes and another for branches without changes";
			message+="\nOU1     = Blomberg et al. Ornstein-Uhlenbeck (one attraction parameter (d), one mean)";
			message+="\nACDC    = Blomberg et al. Acceleration/Deceleration (g parameter)";
			message+="\nOUSM    = Ornstein-Uhlenbeck with one mean per discrete state (attraction and rate parameters constant across tree)";
			message+="\nOUCM    = Ornstein-Uhlenbeck with independent means on branches with and without changes in a discrete character";
			message+="\n         (attraction and rate parameters constant across tree)";
			//message+="\nKappa   = Pagel's Kappa parameter (continuous char)";
			//message+="\nDelta   = Pagel's Delta parameter (continuous char)"; //NOTE: Calculation of this seems off
			//message+="\nLambda  = Pagel's Lambda parameter (continuous char)"; //NOTE: Calculation of this seems off

		//	message+="\nOUSMVA = Ornstein-Uhlenbeck, with different mean, variance, and attraction parameters for each state on a tree";
		//	message+="\nOUCMVA = Ornstein-Uhlenbeck, with different mean, variance, and attraction parameters on branches with and without state changes";
		//	message+="\nOUSMA  = Ornstein-Uhlenbeck, with different mean and attraction parameters for each state, but one variance parameter";
		//	message+="\nOUCMA  = Ornstein-Uhlenbeck, with different mean and attraction parameters on branches with and without state changes, but one variance parameter";
		//	message+="\nTSMBI  = Brownian motion, with different rate parameters for each time interval (set using timeslice)";
		//	message+="\nTSMBC  = Brownian motion, with different rate parameters for branches that cross and do not cross time interval boundaries";
			
			message+="\n\nState vector allows restrictions, so that character states 0 and 2, for example, may be\n     viewed by the program as identical. To do this, you'd enter:\n\n       states=(0 1 0 2 3 4 5 6 7 8)";
            PrintMessage();
        }
        else if( token.Abbreviation("Type") ) {
			citationarray[0]=true;
            nxsstring chosenmodelinput=GetFileName(token);
            if (token.Abbreviation("BM1")) {
                if (tipvariancetype==1) {
                    chosenmodel=2;
                }
                else {
                    chosenmodel=1;
                }
                message="You have chosen a single rate Brownian motion model";
                PrintMessage();
            }
            else if(token.Abbreviation("OU1")) {
                chosenmodel=3;
                citationarray[1]=true;
                message="You have chosen the OU1 model";
                PrintMessage();
            }
            else if(token.Abbreviation("ACDC1")) {
                chosenmodel=4;
                citationarray[1]=true;
                message="You have chosen the ACDC model";
                PrintMessage();
            }
            else if(token.Abbreviation("BMS")) {
                chosenmodel=5;
				citationarray[0]=true;
                message="You have chosen Brownian motion with one rate parameter per state";
                PrintMessage();
            }
            else if(token.Abbreviation("BMC")) {
                chosenmodel=6;
				citationarray[0]=true;
                message="You have chosen Brownian motion with different rates on branches with and without state changes.";
                PrintMessage();
            }
            else if(token.Abbreviation("OUSMVA")) {
                chosenmodel=7;
                message="You have chosen the OU model with different means, brownian rate parameters, and attraction parameters for each state.";
                PrintMessage();
            }
            else if(token.Abbreviation("OUCMVA")) {
                chosenmodel=8;
                message="You have chosen the OU model with different means, brownian rate parameters, and attraction parameters on branches with and without state changes.";
                PrintMessage();
            }
            else if(token.Abbreviation("OUSMA")) {
                chosenmodel=10;
                message="You have chosen the OU model with different means and attraction parameters for each state.";
                PrintMessage();
            }
            else if(token.Abbreviation("OUCMA")) {
                chosenmodel=11;
                message="You have chosen the OU model with different means and attraction parameters on branches with and without state changes.";
                PrintMessage();
            }
            else if(token.Abbreviation("OUSM")) {
                chosenmodel=12;
                message="You have chosen the OU model with different means for each state.";
				citationarray[2]=true;
                PrintMessage();
            }
            else if(token.Abbreviation("OUCM")) {
                chosenmodel=13;
                message="You have chosen the OU model with different means on branches with and without state changes.";
				citationarray[2]=true;
                PrintMessage();
            }
            else if(token.Abbreviation("TSBMI")) {
                chosenmodel=14;
                message="You have chosen to have a different rate parameter in each time interval.";
                PrintMessage();
            }
            else if(token.Abbreviation("TSBMC")) {
                chosenmodel=15;
                message="You have chosen to have a different rate parameter on branches that do and do not cross a time interval.";
                PrintMessage();
            }
            else if (token.Abbreviation("BMAN")) {
                chosenmodel=16;
                message="You have chosen to try all possible assignments of zero rates to branches, leaving the other brlen unchanged";
                PrintMessage();
            }
            else if (token.Abbreviation("BMAO")) {
                chosenmodel=17;
                message="You have chosen to try all possible assignments of zero rates to branches, assigning a branch length of one to the other branches";
                PrintMessage();
            }
            else if (token.Abbreviation("BMPN")) {
                chosenmodel=18;
                message="You have chosen to try, at each internal node, each assignment of length 0 to one descendant branch and the original brlen to the other";
                PrintMessage();
            }
            else if (token.Abbreviation("BMPO")) {
                chosenmodel=19;
                message="You have chosen to try, at each internal node, each assignment of length 0 to one descendant branch and a length of 1 to the other";
                PrintMessage();
            }
			else if (token.Abbreviation("KAppa")) {
                chosenmodel=20;
                message="You have chosen to try Pagel's Kappa parameter model";
                PrintMessage();
            }
			else if (token.Abbreviation("DElta")) {
                chosenmodel=21;
                message="You have chosen to try Pagel's Delta parameter model";
                PrintMessage();
            }
			else if (token.Abbreviation("LAmbda") || token.Abbreviation("LAmda")) {
                chosenmodel=22;
                message="You have chosen to try Pagel's Lambda parameter model";
                PrintMessage();
            }
            else {
                errormsg = "Unexpected option (";
                errormsg += chosenmodelinput;
                errormsg += ") encountered reading Model command";
                throw XNexus( errormsg);
            }
        }
        else if( token.Abbreviation("States") ) {
            token.GetNextToken();
            token.GetNextToken(); //eat the equals sign
            vector<int> temporarystatevector;
            if (!token.Equals("(")) {
                errormsg="Expecting next token to be a \'(\' but instead got ";
                errormsg+=token.GetToken();
                throw XNexus( errormsg);
            }
            int inputcount=0;
            while (!token.Equals(")")) {
                nxsstring numbernexus;
                numbernexus=GetNumberOnly(token);
                if (numbernexus!=")") {
                    temporarystatevector.push_back(atoi( numbernexus.c_str() ));
                    inputcount++;
                }
                else {
                    break;
                }
            }
            if (inputcount<1) {
                errormsg="You should have entered at least one number";
                throw XNexus( errormsg);
            }
            else {
                int maxstate=0;
                for (int i=0;i<GSL_MIN(temporarystatevector.size(),staterestrictionvector.size());i++) {
                    staterestrictionvector[i]=temporarystatevector[i];
                    
                }
                for (int i=0;i<staterestrictionvector.size();i++) {
                    maxstate=GSL_MAX(staterestrictionvector[i],maxstate);
                }
                bool complete=false;
                while (!complete) {
                    for (int j=0;j<=maxstate;j++) {
                        int count=0;
                        for (int i=0;i<staterestrictionvector.size();i++) {
                            if(staterestrictionvector[i]==j) {
                                count++;
                            }
                        }
                        if (count==0) {
                            for (int i=0;i<staterestrictionvector.size();i++) {
                                if(staterestrictionvector[i]>=j) {
                                    staterestrictionvector[i]=staterestrictionvector[i]-1;
                                }
                            }
                            maxstate--;
                            complete=false;
                            break;
                        }
                        else {
                            complete=true;
                        }
                    }
                }
                //int sizediff=staterestrictionvector.size()-temporarystatevector.size();
               //// if (sizediff>0) {
               //     for (int i=temporarystatevector.size();i<staterestrictionvector.size();i++) {
               //         staterestrictionvector[i]=i-1-(-1+temporarystatevector.size()-maxstate);
               //     }
              //  }
            }
			
            ///below here is old way of reading vector
          //  nxsstring numbernexus;
           // numbernexus=GetNumber(token);
          //  staterestrictionvector[0]=atoi( numbernexus.c_str() );
         //   for (int position=1;position<=9;position++) {
         //       numbernexus=GetNumberOnly(token);
         //       if (debugmode) {
         //           message="Now reading for position ";
          //          message+=position;
          //          message+=" state is ";
          //          message+=numbernexus;
          //          PrintMessage();
          //      }
           //     staterestrictionvector[position]=atoi( numbernexus.c_str() );
           // }
        }
		
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Model command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

void BROWNIE::HandleLoss ( NexusToken& token )
{
	// Retrieve all tokens for this command, stopping only in the event
	// of a semicolon or an unrecognized keyword
	//
	bool donesomething=false;
    for(;;)
    {
        token.GetNextToken();
		
        if( token.Equals(";") ) {
			if (donesomething==false) {
			//double rateA=0.0000000000000000001;
			//while (rateA<100) {
			//	double neglnL=CalculateDiscreteLindy1(rateA);
			//	message="rateA=";
			//	message+=rateA;
			//	message+=" -lnL=";
			//	message+=neglnL;
			//	PrintMessage();
			//	rateA*=100.0;
			//}
			//double neglnL=CalculateDiscreteLindy1(0.022322);
			//message="rateA=";
			//message+=0.022322;
			//message+=" -lnL=";
			//message+=neglnL;
			//PrintMessage();
				
			//BROWNIE my_fnA(maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput,trees, taxa, assumptions, characters);
				gsl_vector *optimalrateA=gsl_vector_calloc(3);
				//gsl_vector_memcpy(optimalrate,my_fn.OptimizeRateWithGivenTipVariance());
				gsl_vector_memcpy(optimalrateA,LindyGeneralOptimization(1));
			//for(int i=0;i<optimalrateA->size;i++){
			//	message="i=";
			//	message+=i;
			//	message+=" val is ";
			//	message+=gsl_vector_get(optimalrateA,i);
			//	PrintMessage();
			//}
			//BROWNIE my_fnB(maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput,trees, taxa, assumptions, characters);
				gsl_vector *optimalrateB=gsl_vector_calloc(5);
				gsl_vector_memcpy(optimalrateB,LindyGeneralOptimization(2));
		//	for(int i=0;i<optimalrateB->size;i++){
		//		message="i=";
		//		message+=i;
		//		message+=" val is ";
		//		message+=gsl_vector_get(optimalrateB,i);
		//		PrintMessage();
		//	}
				double modelAaicc=(2.0*gsl_vector_get(optimalrateA,2))+2+4.0/((characters->GetNChar())-2); //AICc, n=1;
				double modelBaicc=(2.0*gsl_vector_get(optimalrateB,4))+4+12.0/((characters->GetNChar())-3);
				double modelAaic=(2.0*gsl_vector_get(optimalrateA,2))+2*1;
				double modelBaic=(2.0*gsl_vector_get(optimalrateB,4))+2*2;				
				message="ModelA\n\trate = ";
				message+=gsl_vector_get(optimalrateA,0);
				message+=" +/- ";
				message+=gsl_vector_get(optimalrateA,1);
				message+="\n\tlnL = ";
				message+=gsl_vector_get(optimalrateA,2);
				message+="\n\tAIC = ";
				message+=modelAaic-GSL_MIN(modelAaic,modelBaic);
				message+="\n\tAICc = ";
				message+=modelAaicc-GSL_MIN(modelAaicc,modelBaicc);;				
				message+="\nModelB\n\trate0 = ";
				message+=gsl_vector_get(optimalrateB,0);
				message+=" +/- ";
				message+=gsl_vector_get(optimalrateB,2);
				message+="\n\trate1 = ";
				message+=gsl_vector_get(optimalrateB,1);
				message+=" +/- ";
				message+=gsl_vector_get(optimalrateB,3);
				message+="\n\tlnL = ";
				message+=gsl_vector_get(optimalrateB,4);
				message+="\n\tAIC = ";
				message+=modelBaic-GSL_MIN(modelAaic,modelBaic);
				message+="\n\tAICc = ";
				message+=modelBaicc-GSL_MIN(modelAaicc,modelBaicc);	
				PrintMessage();
			}
            break;
        }
        else if( token.Abbreviation("?") ) {
			donesomething=true;
            message="Usage: Loss \n\n";
            message+="";
            message+="Available options: none\n\n";
		message+="This evaluates two models. For each model, the ancestral state is assumed to be 1 (presence) and the model is applied across all characters. \n  ModelA: Same loss rate on all branches\n  ModelB: Independent loss rates on branches with label 0 and label 1.";
            PrintMessage();
        }
		
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Loss command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}


void BROWNIE::HandleDiscrete( NexusToken& token )
{
	// Retrieve all tokens for this command, stopping only in the event
	// of a semicolon or an unrecognized keyword
	//
    for(;;)
    {
        token.GetNextToken();
		
        if( token.Equals(";") ) {
			double rateA=0.0000000000000000001;
			while (rateA<100) {
				double neglnL=CalculateDiscreteLindy1(rateA);
				message="rateA=";
				message+=rateA;
				message+=" -lnL=";
				message+=neglnL;
				PrintMessage();
				rateA*=100.0;
			}
			double neglnL=CalculateDiscreteLindy1(0.022322);
			message="rateA=";
			message+=0.022322;
			message+=" -lnL=";
			message+=neglnL;
			PrintMessage();
			
			//BROWNIE my_fnA(maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput,trees, taxa, assumptions, characters);
			gsl_vector *optimalrateA=gsl_vector_calloc(3);
				//gsl_vector_memcpy(optimalrate,my_fn.OptimizeRateWithGivenTipVariance());
			gsl_vector_memcpy(optimalrateA,LindyGeneralOptimization(1));
			for(int i=0;i<optimalrateA->size;i++){
					message="i=";
				message+=i;
				message+=" val is ";
				message+=gsl_vector_get(optimalrateA,i);
				PrintMessage();
			}
			//BROWNIE my_fnB(maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput,trees, taxa, assumptions, characters);
			gsl_vector *optimalrateB=gsl_vector_calloc(5);
			gsl_vector_memcpy(optimalrateB,LindyGeneralOptimization(2));
			for(int i=0;i<optimalrateB->size;i++){
				message="i=";
				message+=i;
				message+=" val is ";
				message+=gsl_vector_get(optimalrateB,i);
				PrintMessage();
			}
			
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="Usage: Discrete type = [ set | optimize ] states = [empirical | set | optimize] treeloop=[yes|no] allchar=[yes|no]\n\n";
            message+="";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ----------------------------- Current setting -----";
            message+="\n Type        Set                                       Set ";
			message+="\n States                                                AncState=1";
			message+="\n AllChar     Y/N                                       Yes";
			message+="\n\nCurrent settings have ancestral state 1, has one model on branches with one state label\nand a different model on branches with another state label.\nThe same model is applied to all chars.";
            PrintMessage();
        }
		
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Discrete command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}


void BROWNIE::FindFixedDiscreteModel() {
	
	
}




/**Make combined VCV
*/
void BROWNIE::MakeCombinedVCV(gsl_matrix *VCVcombined, gsl_matrix *VCVtoadd, int ntaxprocessed)
{
    int ntaxtoadd=VCVtoadd->size1;
    gsl_matrix_view newview= gsl_matrix_submatrix(VCVcombined,ntaxprocessed,ntaxprocessed,ntaxtoadd,ntaxtoadd);
    gsl_matrix_memcpy(&newview.matrix,VCVtoadd);
}

//Simulate tip values
gsl_vector* BROWNIE::SimulateTips(gsl_matrix * VCV, double rate, gsl_vector *MeanValues)
{
    //Code inspired by John Burkardt, also based on code from Handbook of Simulation: Principles, Methodology, Advances, Applications, and Practice, Jerry Banks, ed. 1998.
    int ntax=VCV->size1;
    gsl_vector *newtips=gsl_vector_calloc(ntax);
    gsl_matrix *A=gsl_matrix_calloc(ntax,ntax);
    int CopyResult= gsl_matrix_memcpy(A, VCV);
	gsl_matrix_scale (A,rate);
    int CholResult= gsl_linalg_cholesky_decomp( A); //A has LU in upper and lower diagonals
	gsl_vector *randomvect=gsl_vector_calloc(ntax);
	for (int i=0;i<ntax;i++) {
		gsl_vector_set(randomvect,i,gsl_ran_ugaussian(r));
	}
	gsl_blas_dgemv (CblasNoTrans,1, A, randomvect,0, newtips);
	gsl_matrix_free(A);
	gsl_vector_free(randomvect);
    return newtips;

}


//Convert VCV matrix (n x n) to VCV vector (length n^2) and append it to a pre-existing vector
// if you don't have a pre-existing vector, create one: gsl_vector OrigVector(0,0);
//gsl_vector* BROWNIE::ConvertVCVMatrixToVector(gsl_matrix *VCV, gsl_vector *OrigVector)
//{
//    int OrigVectorLength=OrigVector->size;
//   int VCVMatrixArea=(VCV->size1)*(VCV->size2);
//    gsl_matrix *VCVcomb=gsl_matrix_alloc(0,0);
//    gsl_matrix *OrigMatrix=gsl_matrix_alloc(0,0);
//    gsl_vector *OutVector;
//   OutVector=gsl_vector_calloc(0);
//    OrigMatrix=ConvertVCVVectorToMatrix(OrigVector);
//    VCVcomb=MakeCombinedVCV(OrigMatrix,VCV);
//    int ntax=VCVcomb->size1;
//    for (int r=0;r<ntax;r++) {
//        for (int c=0;c<ntax;c++) {
//            OutVector.push_back(VCVcomb[r][c]);
//        }
//    }
//    return OutVector;
//}

//Convert VCV vector (length n^2) to VCV matrix (n x n)
//gsl_matrix* BROWNIE::ConvertVCVVectorToMatrix(gsl_vector *VCVvector)
//{
//    double ntaxsquared=VCVvector->size;
//    int ntax=int(sqrt(ntaxsquared));
//    gsl_matrix *VCV=gsl_matrix_calloc (ntax,ntax);
//    for (int r=0;r<ntax;r++) {
//        for (int c=0;c<ntax;c++) {
//            VCV[r][c]=VCVvector[(r*ntax)+c];
//        }
//    }
//    return VCV;
//}

    //Extract one VCV matrix from a VCV combined matrix that is represented as a vector
    //gsl_matrix * BROWNIE::ExtractMatrixFromVector(gsl_vector *VCVvector, int ntaxprocessed, int currentntax) {
    //    gsl_matrix *VCVentiremat=gsl_matrix_calloc(0,0);
    //    gsl_matrix *VCVextracted=gsl_matrix_calloc(currentntax,currentntax);
    //    VCVentiremat=ConvertVCVVectorToMatrix(VCVvector);
    //    for (int r=0;r<currentntax;r++) {
    //        for (int c=0;c<currentntax;c++) {
    //            VCVextracted[r][c]=VCVentiremat[ntaxprocessed+r][ntaxprocessed+c];
    //        }
    //    }
    //    return VCVextracted;
    //}


    /**Make combined tips
    */
    //void BROWNIE::MakeCombinedTips(gsl_vector *tipscombined, gsl_vector *tipstoadd, int ntaxprocessed)
    void BROWNIE::MakeCombinedTips(gsl_vector *tipscombined, gsl_vector *tipstoadd, int ntaxprocessed)
{
        int ntaxtoadd=tipstoadd->size;
        gsl_vector_view newview=gsl_vector_subvector(tipscombined,ntaxprocessed,ntaxtoadd);
        gsl_vector_memcpy(&newview.vector,tipstoadd);
}

//Returns log likelihood. If the VCV matrix includes other components (like tip variance), deal with these AND THE RATE first, and just pass a rate of 1 to this function.
double BROWNIE::GetLScore(gsl_matrix *VCV,gsl_vector *tipresid,double rate){
    int ntax=VCV->size1;
    double lscore;
    gsl_matrix *RateTimesVCV=gsl_matrix_calloc(ntax,ntax);
    //gsl_matrix RateTimesVCV(ntax,ntax,0);
    gsl_matrix *InverseRateTimesVCVstart=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix *InverseRateTimesVCV=gsl_matrix_calloc(ntax,ntax);
    //gsl_matrix InverseRateTimesVCV(ntax,ntax,0);
    //long double DetRateTimesVCV;
    gsl_matrix_memcpy (RateTimesVCV, VCV);

    if(debugmode) {
        message="debugging line 957\n";
        for (int currentrow=0;currentrow<VCV->size1;currentrow++) {
            for (int currentcol=0;currentcol<VCV->size2;currentcol++) {
                message+=gsl_matrix_get(VCV,currentrow,currentcol);
                message+="\t";
            }
            message+="\n";
        }
        message+="\n\n RateTimesVCV\n";
        for (int currentrow=0;currentrow<RateTimesVCV->size1;currentrow++) {
            for (int currentcol=0;currentcol<RateTimesVCV->size2;currentcol++) {
                message+=gsl_matrix_get(RateTimesVCV,currentrow,currentcol);
                message+="\t";
            }
            message+="\n";
        }

        PrintMessage();
    }


    gsl_matrix_scale (RateTimesVCV,rate);
    if(debugmode) {
        message="debugging line 980\n";
        message+="\nrate is ";
        message+=rate;
        for (int currentrow=0;currentrow<RateTimesVCV->size1;currentrow++) {
            for (int currentcol=0;currentcol<RateTimesVCV->size2;currentcol++) {
                message+=gsl_matrix_get(RateTimesVCV,currentrow,currentcol);
                message+="\t";
            }
            message+="\n";

        }

        PrintMessage();
    }


    //RateTimesVCV=rate*VCV;
    //matrixsingular=TestSingularity(RateTimesVCV);
    //if (matrixsingular) {
    //    errormsg="Singular matrix (RateTimesVCV) during GetLScore";
    //    throw XNexus(errormsg);
    //}
    if (debugmode) {
        message="\nNow calculating inverse VCV for estimating lnL.\n";
        PrintMessage();
    }
//USE THESE

gsl_permutation * p = gsl_permutation_alloc (ntax);
int signum;
gsl_matrix_memcpy (InverseRateTimesVCVstart, RateTimesVCV);
if(debugmode) {
    message="debugging\n";
    for (int currentrow=0;currentrow<RateTimesVCV->size1;currentrow++) {
        for (int currentcol=0;currentcol<RateTimesVCV->size2;currentcol++) {
            message+=gsl_matrix_get(RateTimesVCV,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }
    PrintMessage();
}
if(debugmode) {
    message="ntax is ";
    message+=ntax;
    message+=" and size of InverseRateTimesVCVstart is ";
    int irtvssize=InverseRateTimesVCVstart->size1;
    message+=irtvssize;
    int permsize=p->size;
    message+=" size of perm is ";
    message+=permsize;
    PrintMessage();
}
gsl_linalg_LU_decomp (InverseRateTimesVCVstart,p, &signum);
gsl_linalg_LU_invert (InverseRateTimesVCVstart,p, InverseRateTimesVCV);
//gsl_permutation_free (p);
//InverseRateTimesVCV=Inverse(RateTimesVCV);
if (debugmode) {
    message="\nFinished calculating inverse VCV for estimating rate.\n";
    PrintMessage();
    cout<<"RateTimesVCV\n";
    cout<<RateTimesVCV;
    cout<<"\n";
}

if(debugmode) {
    message="ntax is ";
    message+=ntax;
    message+=" and size of InverseRateTimesVCVstart is ";
    int irtvssize=InverseRateTimesVCVstart->size1;
    message+=irtvssize;
    int permsize=p->size;
    message+=" size of perm is ";
    message+=permsize;
    PrintMessage();
}

gsl_matrix *tipsasmatrixprime=gsl_matrix_calloc(1,ntax);
//gsl_matrix tipsasmatrixprime(1,ntax,0);
for (int i=0; i<ntax; i++) {
    gsl_matrix_set(tipsasmatrixprime, 0, i, gsl_vector_get(tipresid,i));
    // tipsasmatrixprime[0][i]=tipresid[i];
}
gsl_vector *step1vect;
step1vect=gsl_vector_calloc(ntax);
//gsl_vector step1vect(ntax,0);
gsl_blas_dgemv (CblasNoTrans,1, InverseRateTimesVCV, tipresid,0, step1vect); ///TEST THIS
                                                                             //step1vect=MatrixTimesVector(InverseRateTimesVCV,tipresid);
gsl_vector *step2vect;
step2vect=gsl_vector_calloc(1);
//gsl_vector step2vect(1,0);
gsl_blas_dgemv (CblasNoTrans,1, tipsasmatrixprime, step1vect,0, step2vect); ///TEST THIS
                                                                            //step2vect=MatrixTimesVector(tipsasmatrixprime,step1vect);
                                                                            //DetRateTimesVCV=Determinant(RateTimesVCV);
                                                                            //lscore=-1*log((exp(-0.5*step2vect[0]))/(sqrt(DetRateTimesVCV*(pow((2*PI),ntax)))));
gsl_matrix *RateTimesVCVLU=gsl_matrix_calloc(ntax,ntax);
gsl_matrix_memcpy (RateTimesVCVLU, RateTimesVCV);
if(debugmode) {
    message="debugging\n";
    for (int currentrow=0;currentrow<RateTimesVCV->size1;currentrow++) {
        for (int currentcol=0;currentcol<RateTimesVCV->size2;currentcol++) {
            message+=gsl_matrix_get(RateTimesVCV,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }
    PrintMessage();
}
gsl_linalg_LU_decomp (RateTimesVCVLU,p, &signum);
lscore=0.5*gsl_vector_get(step2vect,0)+0.5*(gsl_linalg_LU_lndet(RateTimesVCVLU))+0.5*ntax*log(2*PI);
if (debugmode) {
    message="Likelihood score is ";
    message+=lscore;
    PrintMessage();
}
gsl_matrix_free(RateTimesVCV);
gsl_matrix_free(InverseRateTimesVCVstart);
gsl_matrix_free(InverseRateTimesVCV);
gsl_matrix_free(tipsasmatrixprime);
gsl_vector_free(step1vect);
gsl_vector_free(step2vect);
gsl_matrix_free(RateTimesVCVLU);
gsl_permutation_free (p);
return lscore;
}




/** @method HandleRateTest
*
*
*/
void BROWNIE::HandleRateTest( NexusToken& token )
{
	citationarray[0]=true;
    bool noerror=true;
    //int seed=time(NULL);
    //int *seedptr=&seed;
    //srand(time(0));
    typedef std::set<nxsstring, std::less < nxsstring> > nxsstring_set;
    nxsstring tmessage;
    ofstream tablef;
    nxsstring tablefname;
    bool tablef_open=false;
    bool name_provided=false;
    nxsstring_set chosentaxSETS;
    //gsl_vector *chosentaxsetNTAXvector=gsl_vector_calloc(0);
    //gsl_vector chosentaxsetNTAXvector(0,0);
    set<int> chosentaxNTAX;
    map<nxsstring,int> chosentaxsetntaxmap;
    //set<gsl_vector> chosentaxCholVCVvectors;
    // gsl_vector vectorofCholVCVvectors;
    // map<int,gsl_vector> mapofCholVCVvectors;
    //gsl_vector *CholVCVvectorvector;
    // CholVCVvectorvector=gsl_vector_calloc(0);
    //gsl_vector CholVCVvectorvector(0,0);
    //gsl_vector *VCVvectorvector;
    //VCVvectorvector=gsl_vector_calloc(0);
    //gsl_vector VCVvectorvector(0,0);

    nxsstring chosentaxset;
    int repsnumber=0;
    int listedtaxsets=0;
    int ntaxcomb=0;
    bool treeloop=false;
    bool charloop=false;
    bool adequateinput=false;
    bool notquietmode=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (adequateinput==false) {
                message="Insufficient input: type \"ratetest ?\" for help";
                PrintMessage();
            }
            break;
        }
        else if( token.Equals("?") ) {
            adequateinput=true;
            message="Usage: RateTest taxset=<taxset 1> [taxset=<taxset 2>...] [options...]\n\n";
            message+="Performs a censored rate test for two or more subtrees.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\nTaxset       <taxset name>                        *None";
            message+="\nReps         <integer>                            *0";
            message+="\nTreeloop     No|Yes                               *No";
            message+="\nCharloop     No|Yes                               *No";
            message+="\nQuiet        No|Yes                               *No";
            message+="\nFile         <file name>                          *None";
            message+="\n                                                 *Option is nonpersistent\n\n";
            message+="This will compare the likelihood under a single rate Brownian motion\nmodel with the likelihood under a multiple parameter Brownian motion\nmodel using the 'censored' test using the current tree (use 'choose' to\nchange the current tree). If just one taxset is given, it will assign\none rate to the pruned tree containing just those taxa and another rate\nto the pruned tree containing all the other taxa. If multiple taxsets\nare given, it will assign one rate to each pruned tree and test this\nagainst the simple model where there is one rate for all these pruned\ntrees. \n\nNote that this method currently does not do much error checking: it is\nup to the user to check that the pruned trees represented by the taxsets\nare clades or paraphyletic with respect to the other pruned trees and\nthat no taxon is in multiple taxsets. If you are just entering one\ntaxset, and its pruned tree is a clade, you should have nothing to worry\nabout, even if its complement is paraphyletic.\n\nThe 'reps' statement is also optional. This tells the program how many\ntimes it should simulate data under the null model when doing parametric\nbootstrapping to test for significance using a likelihood ratio test.\nThe default is 0. When set to zero, it avoids parametric\nbootstrapping altogether.\n\nTreeloop and charloop tell the program whether to loop across all\ntrees and/or characters or just use the currently-selected ones.\nThe default for both is no.\nIf quiet=yes, only the summary is printed.\n\n\nYou may want to have logging activated to record the output from ratetest.";
            PrintMessage();
        }
        else if( token.Abbreviation("File") ) {
            tablefname = GetFileName(token);
            name_provided = true;
            bool exists = FileExists( tablefname.c_str() );
            bool userok = true;
            if( exists && !UserSaysOk( "Ok to replace?", "Ratetest output file specified already exists" ) )
                userok = false;
            if( userok ) {
                tablef_open = true;
                tablef.open( tablefname.c_str() );
            }

            if( exists && userok ) {
                message = "\nReplacing ratetest output file ";
                message += tablefname;
            }
            else if( userok ) {
                message = "\nRatetest output file ";
                message += tablefname;
                message += " opened";
            }
            else {
                errormsg = "Aborting the ratetest so as not to overwrite the file.\n";
                throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );

            }
            PrintMessage();

        }

        else if (token.Abbreviation("Taxset") ) {
            adequateinput=true;
            int ntax=0;
            if (trees->GetNumTrees()<1) {
                errormsg = "Error: No valid trees are loaded.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            chosentaxset=GetFileName(token);
            chosentaxSETS.insert(chosentaxset);

            listedtaxsets++;
            IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
            if (taxonlist.empty()) {
                errormsg= "Error: Taxset ";
                errormsg+=chosentaxset.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntax++;
            }
            chosentaxsetntaxmap[chosentaxset]=ntax;
            ntaxcomb+=ntax;
        }
        else if (token.Abbreviation("Reps")  ) {
            nxsstring repschar;
            repschar=GetFileName(token);
            repsnumber=atoi(repschar.c_str());
        }
        else if (token.Abbreviation("TReeloop") ) {
            nxsstring yesnotreeloop=GetFileName(token);
            if (yesnotreeloop[0] == 'n') {
                treeloop=false;
            }
            else {
                treeloop=true;
            }
        }
        else if (token.Abbreviation("Quiet") ) {
            nxsstring yesnoquiet=GetFileName(token);
            if (yesnoquiet[0] == 'n') {
                notquietmode=true;
            }
            else {
                notquietmode=false;
            }
        }
        else if (token.Abbreviation("CHarloop") ) {
            nxsstring yesnocharloop=GetFileName(token);
            if (yesnocharloop[0] == 'n') {
                charloop=false;
            }
            else {
                charloop=true;
            }
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading RateTest command. Type \"RateTest ?\" for help.";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }


    }
    if (listedtaxsets>0) {
        if (listedtaxsets==1) {
            listedtaxsets++;
            nxsstring chosentaxset2="NOT";
            chosentaxset2+=chosentaxset; //so if one taxset is listed, its complement is automatically used.
            IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset2 );
            if (taxonlist.empty()) {
                errormsg= "Error: Taxset ";
                errormsg+=chosentaxset2.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            chosentaxSETS.insert(chosentaxset2);
            int ntax=0;
IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntax++;
            }
            chosentaxsetntaxmap[chosentaxset2]=ntax;
            ntaxcomb+=ntax;
        }



        if (tablef_open) {
            tmessage="Output in tab-delimited table:\nModel A is constrained to have one rate for all groups, Model B has one rate for each group.\n";
            tablef<<tmessage;
        }
        nxsstring_set chosentaxSETS;
        map<nxsstring, int>::const_iterator iter;
        int n=1;
        for (iter=chosentaxsetntaxmap.begin();iter!=chosentaxsetntaxmap.end();++iter) {
            //nxsstring_set::const_iterator cti;
            //for( cti = chosentaxSETS.begin(); cti != chosentaxSETS.end(); cti++ ) {
            tmessage="Taxset ";
            tmessage+=n;
            n++;
            tmessage+="=";
            tmessage+=iter->first;
            //message+=chosentaxSETS(*cti);
            tmessage+="\n";
            if (tablef_open) {
                tablef<<tmessage;
            }
            }
        tmessage="Tree\tTree weight\tTree name\tChar\t";
        for (int n=0; n<listedtaxsets; n++) {
            tmessage+="anc_";
            tmessage+=n+1;
            tmessage+="\trate_";
            tmessage+=n+1;
            tmessage+="\t lnL_";
            tmessage+=n+1;
            tmessage+="\t";
        }
        tmessage+="rate_A\tparam_A\tparam_B\tAIC_A\tAIC_B\tAICc_A\tAICc_B\t lnL_A\t lnL_B\tAIC dif\tAICc diff\tchi p\tparam p\tchosen model under AIC, AICc, chi, param\n";
        if (tablef_open) {
            tablef<<tmessage;
        }
        int starttree, stoptree, startchar, stopchar;
        int originalchosentree=chosentree;
        int originalchosenchar=chosenchar;
        if (treeloop) {
            starttree=1;
            stoptree=trees->GetNumTrees();
        }
        else {
            starttree=chosentree;
            stoptree=chosentree;
        }
        if (charloop) {
            startchar=1;
            stopchar=characters->GetNChar();
        }
        else {
            startchar=chosenchar;
            stopchar=chosenchar;
        }
        //Loop across trees and chars
        nxsstring summaryofresults="Summary of results\n(B/b = strong/weak support for multiple rate model)\nTree\tChar\tAIC\tAICc\tChi\tParam\n";
        if (repsnumber>0) {
            message="\nExpect a delay while doing parametric simulation.\nIf this is taking too long, use the option reps=0 when using the\nratetest command in the future. Note that parametric simulation\nis very important if you are using a p-value approach and do not\nhave many taxa (when the chi-square test is non-conservative).\n";
            PrintMessage();
        }
        double weighttotal=0;
        double weightedAIC1=0;
        double weightedAIC2=0;
        double weightedAICc1=0;
        double weightedAICc2=0;
        gsl_vector *weightedratevector;
        weightedratevector=gsl_vector_calloc(listedtaxsets+1);
        //gsl_vector weightedratevector(listedtaxsets+1,0); //So that the first rate is for the single rate model
        gsl_vector *weightedancstatevector;
        weightedancstatevector=gsl_vector_calloc(listedtaxsets+1);
        //gsl_vector weightedancstatevector(listedtaxsets+1,0); // Just for consistency with above
        double weightedchip=0;
        double weightedparamp=0;


        for (chosentree=starttree;chosentree<=stoptree;chosentree++) {
            int ntaxprocessedtreeloop=0;
            bool goodtree=true;
            //gsl_vector *CholVCVvectorvector;
            //CholVCVvectorvector=gsl_vector_calloc(0);
            //            gsl_vector CholVCVvectorvector(0,0);
            //gsl_vector *VCVvectorvector;
            //VCVvectorvector=gsl_vector_calloc(0);
            //            gsl_vector VCVvectorvector(0,0);
            gsl_matrix *VCVcomb=gsl_matrix_calloc(ntaxcomb,ntaxcomb);
            //gsl_matrix VCVcomb(0,0,0);
            double treeweight=trees->GetTreeWeight(chosentree-1);
            weighttotal+=treeweight;
            nxsstring treename=trees->GetTreeName(chosentree-1);
            map<nxsstring, int>::const_iterator iter;
            for (iter=chosentaxsetntaxmap.begin();iter!=chosentaxsetntaxmap.end();++iter) {
                nxsstring currenttaxset=iter->first;
                int currentntax=iter->second;
                gsl_matrix *VCV=gsl_matrix_calloc(currentntax,currentntax);
                //gsl_matrix VCV(currentntax,currentntax,0);
                VCV=DeleteStem(GetVCV(currenttaxset));
                //cout<<"taxset "<<currenttaxset<<"  VCV"<<endl<<VCV<<endl;
                // matrixsingular=TestSingularity(VCV);
                //if (matrixsingular) {
                //   goodtree=false;
                //    noerror=false;
                //}
                //else {
                //CholVCVvectorvector=ConvertVCVMatrixToVector(Chol(VCV),CholVCVvectorvector);
                //VCVvectorvector=ConvertVCVMatrixToVector(VCV,VCVvectorvector);
                MakeCombinedVCV(VCVcomb,VCV,ntaxprocessedtreeloop);
                ntaxprocessedtreeloop+=currentntax;


                //if (iter==chosentaxsetntaxmap.begin()) { // It's our first loop
                //    VCVcomb=VCV;
                //}
                //else {
                //    VCVcomb=MakeCombinedVCV(VCVcomb,VCV);
                //}
                //}
				gsl_matrix_free(VCV);
            }
if (goodtree) {
	
    //cout<<"VCV Comb="<<endl<<VCVcomb<<endl<<endl<<"VCVvectorvector (as matrix)"<<endl<<ConvertVCVVectorToMatrix(VCVvectorvector)<<endl<<endl;
    //cout<<"VCVvectorvector"<<endl;
    //cout<<"Length of VCVvectorvector "<<VCVvectorvector->size<<endl<<endl;
    //for (int m=0;m<VCVvectorvector->size;m++) {
    //    cout<<VCVvectorvector[m]<<endl;
    //}
	
for (chosenchar=startchar;chosenchar<=stopchar;chosenchar++) {
    int ntaxprocessedcharloop=0;
    gsl_vector *tipscombresid=gsl_vector_calloc(ntaxcomb);
    //tipscombresid=gsl_vector_calloc(0);
    //gsl_vector tipscombresid(0,0);
    message="Tree = ";
    message+=chosentree;
    message+=": ";
    message+=treename;
    message+=", character = ";
    message+=chosenchar;
    message+="\n";
    if (tablef_open) {
        tmessage="";
        tmessage+=chosentree;
        tmessage+="\t";
        tmessage+=treeweight;
        tmessage+="\t";
        tmessage+=treename;
        tmessage+="\t";
        tmessage+=chosenchar;
        tmessage+="\t";
        tablef<<tmessage;
    }
    int ntaxprocessed=0;
    //Do the test, spit out the values
    map<nxsstring, int>::const_iterator iter2;
    double likelihoodmultiparametermodel=0;
    int taxsetnumber=0;
    for (iter2=chosentaxsetntaxmap.begin();iter2!=chosentaxsetntaxmap.end();++iter2) {
        nxsstring currenttaxset=iter2->first;
        int currentntax=iter2->second;
        taxsetnumber++;
        message+="  Taxset = ";
        message+=currenttaxset;
        message+="\n";
        // gsl_vector *currentVCVvect;
        //currentVCVvect=gsl_vector_calloc(currentntax*currentntax);
        //                        gsl_vector currentVCVvect(currentntax*currentntax,0);
        gsl_matrix *currentVCVmat=gsl_matrix_calloc(currentntax,currentntax);
        currentVCVmat=DeleteStem(GetVCV(currenttaxset));
        //gsl_matrix currentVCVmat(currentntax,currentntax,0);
        gsl_vector *tips=gsl_vector_calloc(currentntax);
        //                        gsl_vector tips(currentntax,0);
        gsl_vector *tipsresid=gsl_vector_calloc(currentntax);
        //gsl_vector tipsresid(currentntax,0);
        double ancstate;
        double rate;
        double likelihood;
        //cout<<"Ntaxprocessed "<<ntaxprocessed<<endl;
        //for (int i=0;i<(currentntax*currentntax);i++) {
        //    currentVCVvect[i]=VCVvectorvector[(ntaxprocessed*ntaxprocessed)+i];
        //    cout<<"currentVCVvect[i] "<<currentVCVvect[i]<<endl;
        // }
        //cout<<"Length of currentVCVvect = "<<currentVCVvect->size<<endl;
        //currentVCVmat=ConvertVCVVectorToMatrix(currentVCVvect);
        //currentVCVmat=ExtractMatrixFromVector(VCVvectorvector,ntaxprocessed,currentntax);
        //matrixsingular=TestSingularity(currentVCVmat);
        //cout<<endl<<endl<<"Current taxset "<<currenttaxset<<endl<<"Current ntax "<<currentntax<<endl<<"CurrentVCVmat:\n"<<currentVCVmat<<endl<<endl;

        //if (matrixsingular) {
        //    cerr<<"Matrix currentVCVmat is singular, ratetest 1"<<endl;
        //    break;
        //}
        ntaxprocessed+=currentntax;
        tips=GetTipValues(currenttaxset,chosenchar);
        //cout<<"Tips of length "<<tips->size<<endl<<"Line 1297"<<endl;
        ancstate=GetAncestralState(currentVCVmat,tips);
        //cout<<endl<<endl<<endl<<"CurrentVCVmat:\n"<<currentVCVmat<<endl;
        //cout<<"Anc state = "<<ancstate<<endl;
        tipsresid=GetTipResiduals(tips,ancstate);
        //for (int debugtaxon=0;debugtaxon<currentntax;debugtaxon++) {
        //   cout<<"Taxon "<<debugtaxon+1<<": "<<tips[debugtaxon]<<"\t"<<tipsresid[debugtaxon]<<endl;
        //}

        rate=EstimateRate(currentVCVmat,tipsresid);
        //cout<<"Rate: "<<rate<<endl;
        if (rate<=0) {
            message="Tree ";
            message+=chosentree;
            message+=" excluded: one rate estimate (";
            message+=rate;
            message+=") was nonpositive\n";
            PrintMessage();
            summaryofresults+=chosentree;
            summaryofresults+="\t--Has a nonpositive rate estimate (";
            summaryofresults+=rate;
            summaryofresults+="). No output.--\n";
            if (tablef_open) {
                tmessage="";
                tmessage+=chosentree;
                tmessage+="\t--Has a nonpositive rate estimate. No output.--\n";
                tablef<<tmessage;
            }
            goodtree=false;
            noerror=false;
        }
        else {
            gsl_matrix *RateTimesVCVfortest=gsl_matrix_calloc(currentVCVmat->size1,currentVCVmat->size2);
            gsl_matrix_memcpy(RateTimesVCVfortest, currentVCVmat);
            gsl_matrix_scale(RateTimesVCVfortest,rate);
            //RateTimesVCVfortest=rate*currentVCVmat;
            //matrixsingular=TestSingularity(RateTimesVCVfortest);
            //if (matrixsingular==false) {
            likelihood=GetLScore(currentVCVmat,tipsresid,rate);
            likelihoodmultiparametermodel+=likelihood;
            //tipscombresid=MakeCombinedTips(tipscombresid,tipsresid);
            MakeCombinedTips(tipscombresid,tipsresid,ntaxprocessedcharloop);
            ntaxprocessedcharloop+=currentntax;
            message+="    Anc state = ";
            message+=ancstate;
            message+="\n    Rate = ";
            message+=rate;
            message+="\n    lnL = ";
            message+=likelihood;
            message+="\n";
            if (tablef_open) {
                tmessage="";
                tmessage+=ancstate;
                tmessage+="\t";
                tmessage+=rate;
                tmessage+="\t";
                tmessage+=likelihood;
                tmessage+="\t";
                tablef<<tmessage;
            }
            gsl_vector_set(weightedratevector,taxsetnumber,gsl_vector_get(weightedratevector,taxsetnumber)+(treeweight*rate));
            gsl_vector_set(weightedancstatevector,taxsetnumber,gsl_vector_get(weightedancstatevector,taxsetnumber)+(treeweight*ancstate));
            // }
            // else {
            //     goodtree=false;
            //     noerror=false;
            // }
			gsl_matrix_free(RateTimesVCVfortest);
        }
		gsl_matrix_free(currentVCVmat);
		gsl_vector_free(tips);
		gsl_vector_free(tipsresid);		
    }
    if (goodtree) {
        double ratecomb;
        double likelihoodsingleparametermodel;
        ratecomb=EstimateRate(VCVcomb,tipscombresid);
        likelihoodsingleparametermodel=GetLScore(VCVcomb,tipscombresid,ratecomb);
        double K1=listedtaxsets+1; //One ancestral state parameter for each taxset, plus one rate parameter
        double K2=2*listedtaxsets; //One ancestral state parameter and one rate parameter for each taxset.
        double model1aicc=(2*likelihoodsingleparametermodel)+2*K1+2*K1*(K1+1)/(ntaxcomb-K1-1); //AICc, n=1;
        double model2aicc=(2*likelihoodmultiparametermodel)+2*K2+2*K2*(K2+1)/(ntaxcomb-K2-1);
        double model1aic=(2*likelihoodsingleparametermodel)+2*K1;
        double model2aic=(2*likelihoodmultiparametermodel)+2*K2;
        if ((ntaxcomb-K1-1)==0) {
            message+="\nWARNING: In the single rate parameter model,\n   there are ";
            message+=K1;
            message+=" parameters to estimate and only\n   ";
            message+=ntaxcomb;
            message+=" datapoints. You need more taxa.\nAICc=Inf [division by zero].\n\n";
        }
        if ((ntaxcomb-K2-1)==0) {
            message+="\n\nWARNING: In the multiple rate parameter model,\n   there are ";
            message+=K2;
            message+=" parameters to estimate but only\n   ";
            message+=ntaxcomb;
            message+=" datapoints. You need more taxa.\nAICc=Inf [division by zero].\n\n";
        }
        tmessage="";
        tmessage+=ratecomb;
        tmessage+="\t";
        tmessage+=K1;
        tmessage+="\t";
        tmessage+=K2;
        tmessage+="\t";
        double aicdif=model1aic-model2aic;
        double absaicdif=fabs(aicdif);
        double aiccdif=model1aicc-model2aicc;
        double absaiccdif=fabs(aiccdif);
        double plrtchi=1-gsl_cdf_chisq_P(2*(likelihoodsingleparametermodel-likelihoodmultiparametermodel), ((2*listedtaxsets)-(listedtaxsets+1))); //TEST THIS
        char outputstring[14];
        // if (debugmode) {
        //     cout<<"One param aic="<<model1aic<<" 2 param aic="<<model2aic<<endl;
        //    cout<<"One param aicc="<<model1aicc<<" 2 param aicc="<<model2aicc<<endl;
        //     cout<<"listedtaxsets="<<listedtaxsets<<" K1="<<K1<<" K2="<<K2<<endl;
        //    cout<<"ntaxcomb="<<ntaxcomb<<endl;
        //cout<<"likelihoodsingleparametermodel="<<likelihoodsingleparametermodel<<endl;
        //cout<<"likelihoodmultiparametermodel="<<likelihoodmultiparametermodel<<endl;
        //    cout<<"(ntaxcomb-K1-1)="<<(ntaxcomb-K1-1)<<endl;
        //    cout<<"(ntaxcomb-K2-1)="<<(ntaxcomb-K2-1)<<endl;

        // }
        message+="\n  Single rate parameter model (model A):\n     lnL = ";
        sprintf(outputstring,"%14.6f",likelihoodsingleparametermodel);
        message+=outputstring;
        message+="\n    AIC =  ";
        sprintf(outputstring,"%14.6f",model1aic);
        message+=outputstring;
        message+="\n    AICc = ";
        sprintf(outputstring,"%14.6f",model1aicc);
        message+=outputstring;
        message+="\n    rate = ";
        message+=ratecomb;
        gsl_vector_set(weightedratevector,0,gsl_vector_get(weightedratevector,0)+(treeweight*ratecomb));
        //weightedratevector[0]+=ratecomb*treeweight;
        message+="\n\n  Multiple rate parameter model (model B):\n     lnL = ";
        sprintf(outputstring,"%14.6f",likelihoodmultiparametermodel);
        message+=outputstring;
        message+="\n    AIC =  ";
        sprintf(outputstring,"%14.6f",model2aic);
        message+=outputstring;
        message+="\n    AICc = ";
        sprintf(outputstring,"%14.6f",model2aicc);
        message+=outputstring;
        message+="\n\n";
        weightedAIC1+=treeweight*model1aic;
        weightedAIC2+=treeweight*model2aic;
        weightedAICc1+=treeweight*model1aicc;
        weightedAICc2+=treeweight*model2aicc;


        if (tablef_open) {
            sprintf(outputstring,"%14.6f",model1aic);
            tmessage+=outputstring;
            tmessage+="\t";
            sprintf(outputstring,"%14.6f",model2aic);
            tmessage+=outputstring;
            tmessage+="\t";
            sprintf(outputstring,"%14.6f",model1aicc);
            tmessage+=outputstring;
            tmessage+="\t";
            sprintf(outputstring,"%14.6f",model2aicc);
            tmessage+=outputstring;
            tmessage+="\t";
            sprintf(outputstring,"%14.6f",likelihoodsingleparametermodel);
            tmessage+=outputstring;
            tmessage+="\t";
            sprintf(outputstring,"%14.6f",likelihoodmultiparametermodel);
            tmessage+=outputstring;
            tmessage+="\t";
            tmessage+=aicdif;
            tmessage+="\t";
            tmessage+=aiccdif;
            tmessage+="\t";
            sprintf(outputstring,"%3g",plrtchi);
            tmessage+=outputstring;
            tmessage+="\t";
            tablef<<tmessage;
        }
        double parametricbootstrappingpvalue=0;
        if (repsnumber==0) {
            tmessage="NA\t";
        }
        else { //Do parametric simulation
               //int *seedptr=&seed;
			ProgressBar(repsnumber);
            for(int i=0; i<repsnumber; i++) {
				ProgressBar(0);
                int ntaxprocessed=0;
                gsl_vector *simtipscombined=gsl_vector_calloc(ntaxcomb);
                //simtipscombined=gsl_vector_calloc(0);
                //gsl_vector simtipscombined(0,0);
                double simlikelihoodmultiparametermodel=0;
                gsl_vector *simtipresidcomb=gsl_vector_calloc(ntaxcomb);
                //simtipresidcomb=gsl_vector_calloc(0);
                //gsl_vector simtipresidcomb(0,0);
                map<nxsstring, int>::const_iterator iter3;
                for (iter3=chosentaxsetntaxmap.begin();iter3!=chosentaxsetntaxmap.end();++iter3) {
                    nxsstring taxset=iter3->first;
                    int ntax=iter3->second;
                    double *simtipoutput;
                    //double cholvect[ntax*ntax];
                    //gsl_vector *cholvect;
                    // cholvect=gsl_vector_calloc(0);
                    // gsl_vector cholvect(0,0);
                    //gsl_vector currentVCVvect(ntax*ntax,0);
                    gsl_matrix *currentVCVmat=gsl_matrix_calloc(ntax,ntax);
                    currentVCVmat=DeleteStem(GetVCV(taxset)); //NOTE that this will take a lot of time (generate the VCV for each taxset for each rep). Better to store the VCVs once.
                                                              //currentVCVmat=ExtractMatrixFromVector(VCVvectorvector,ntaxprocessed,ntax);
                                                              //cholvect=ConvertVCVMatrixToVector(ExtractMatrixFromVector(CholVCVvectorvector,ntaxprocessed,ntax),cholvect);
                                                              //double cholvectdouble[ntax*ntax];
                                                              //for (int i=0;i<(ntax*ntax);i++) {
                                                              //  cholvectdouble[i]=cholvect[i];
                                                              //      currentVCVvect[i]=VCVvectorvector[(ntaxprocessed*ntaxprocessed)+i];
                                                              // }
                                                              //  currentVCVmat=ConvertVCVVectorToMatrix(currentVCVvect);

                    gsl_vector *nullmean=gsl_vector_calloc(ntax);
                    gsl_vector *simtipvector=gsl_vector_calloc(ntax);
                    simtipvector=SimulateTips(currentVCVmat, ratecomb, nullmean);
                    gsl_vector *simtipsresid=gsl_vector_calloc(ntax);
                    //gsl_vector simtipsresid(ntax,0);
                    double simancstate;
                    double simrate;
                    double simlikelihood;
                    simancstate=GetAncestralState(currentVCVmat,simtipvector);
                    simtipsresid=GetTipResiduals(simtipvector,simancstate);
                    simrate=EstimateRate(currentVCVmat,simtipsresid);
					//cout<<simrate<<"\t";
                    if (simrate==0) {
                        cerr<<"Warning: Rate in one parametric simulation was zero.\nChanging to a very tiny number";
                        simrate=1.0e-10;
                    }
                    simlikelihood=GetLScore(currentVCVmat,simtipsresid,simrate);
                    simlikelihoodmultiparametermodel+=simlikelihood;
                    MakeCombinedTips(simtipresidcomb,simtipsresid,ntaxprocessed);
                    ntaxprocessed+=ntax;
					gsl_matrix_free(currentVCVmat);
					gsl_vector_free(nullmean);
					gsl_vector_free(simtipvector);
					gsl_vector_free(simtipsresid);					
                }
                double simratecomb;
                double simlikelihoodsingleparametermodel;
                simratecomb=EstimateRate(VCVcomb,simtipresidcomb);
                simlikelihoodsingleparametermodel=GetLScore(VCVcomb,simtipresidcomb,simratecomb);
				//cout<<"-\tsim like (sim, then mult) "<<simlikelihoodsingleparametermodel<<"\t"<<simlikelihoodmultiparametermodel<<"\t-\treal like (single, multi) "<<likelihoodsingleparametermodel<<"\t"<<likelihoodmultiparametermodel<<"\t-\tsim rate: "<<simratecomb<<"\test rate: "<<ratecomb<<endl;
				//cout<<endl;
				//for (int i=0;i<simtipresidcomb->size;i++) {
				//	cout<<gsl_vector_get(simtipresidcomb,i)<<" ";
				//}
				//cout<<endl;
				//cout<<"\trates\t"<<simratecomb<<"\t"<<ratecomb<<endl;
                if ((simlikelihoodsingleparametermodel-simlikelihoodmultiparametermodel)>=(likelihoodsingleparametermodel-likelihoodmultiparametermodel)) {
                    parametricbootstrappingpvalue+=(1.0/repsnumber); //Need 1.0 rather than 1 or C++ assumes you want an int output
                }
				gsl_vector_free(simtipscombined);
				gsl_vector_free(simtipresidcomb);
            }
            tmessage="";
			sprintf(outputstring,"%3g",parametricbootstrappingpvalue);
            tmessage+=outputstring;
            tmessage+="\t";
            //PrintMessage();
            if (tablef_open) {
                tablef<<tmessage;
            }
        }

        /////////
        //print selected model under various methods
        summaryofresults+=chosentree;
        summaryofresults+="\t";
        summaryofresults+=chosenchar;
        summaryofresults+="\t";

        if (model1aic<model2aic) {
            if (absaicdif<10) {
                message+="a: AIC dif of ";
                message+=absaicdif;
                message+=" somewhat favors the single rate parameter model.\n";
                summaryofresults+="a\t";
            }
            else {
                message+="A: AIC dif of ";
                message+=absaicdif;
                message+=" strongly favors the single rate parameter model.\n";
                summaryofresults+="A\t";

            }
        }
        else {
            if (absaicdif<10) {
                message+="b: AIC dif of ";
                message+=absaicdif;
                message+=" somewhat favors the multiple rate parameter model.\n";
                summaryofresults+="b\t";
            }
            else {
                message+="B: AIC dif of ";
                message+=absaicdif;
                message+=" strongly favors the multiple rate parameter model.\n";
                summaryofresults+="B\t";
            }
        }

        if (model1aicc<model2aicc) {
            if (absaiccdif<10) {
                message+="a: AICc dif of ";
                message+=absaiccdif;
                message+=" somewhat favors the single rate parameter model.\n";
                summaryofresults+="a\t";
            }
            else {
                message+="A: AICc dif of ";
                message+=absaiccdif;
                message+=" strongly favors the single rate parameter model.\n";
                summaryofresults+="A\t";
            }
        }
        else {
            if (absaiccdif<10) {
                message+="b: AICc dif of ";
                message+=absaiccdif;
                message+=" somewhat favors the multiple rate parameter model.\n";
                summaryofresults+="b\t";
            }
            else {
                message+="B: AICc dif of ";
                message+=absaiccdif;
                message+=" strongly favors the multiple rate parameter model.\n";
                summaryofresults+="B\t";

            }
        }

        if (plrtchi>0.05) {
            message+="a: Chi-square p of ";
            message+=plrtchi;
            message+=" does not reject the single rate parameter model.\n";
            summaryofresults+="a\t";
        }
        else {
            if (plrtchi<0.01) {
                message+="B: Chi-square p of ";
                message+=plrtchi;
                message+=" rejects the single rate parameter model.\n";
                summaryofresults+="B\t";
            }
            else {
                message+="b: Chi-square p of ";
                message+=plrtchi;
                message+=" weakly rejects the single rate parameter model.\n";
                summaryofresults+="b\t";

            }
        }
        weightedchip+=plrtchi*treeweight;
        if (repsnumber==0) {
            message+="?: Parametric bootstrapping not done, despite bias in chi-square test.\n";
            summaryofresults+="?\t";
        }
        else {
            weightedparamp+=parametricbootstrappingpvalue*treeweight;
            if (parametricbootstrappingpvalue>0.05) {
                message+="a: Parametric bootstrap p of ";
                message+=parametricbootstrappingpvalue;
                message+=" does not reject the single rate parameter model.\n";
                summaryofresults+="a\t";
            }
            else {
                if (parametricbootstrappingpvalue==0) {
                    message+="B: Parametric bootstrap p of <";
                    message+=1.0/repsnumber;
                    message+=" rejects the single rate parameter model.\n";
                    summaryofresults+="B\t";
                }
                else if (parametricbootstrappingpvalue<0.01) {
                    message+="B: Parametric bootstrap p of ";
                    message+=parametricbootstrappingpvalue;
                    message+=" rejects the single rate parameter model.\n";
                    summaryofresults+="B\t";

                }
                else {
                    message+="b: Parametric bootstrap p of ";
                    message+=parametricbootstrappingpvalue;
                    message+=" weakly rejects the single rate parameter model.\n";
                    summaryofresults+="b\t";
                }
            }
        }


        //print selected model under various methods
        if (model1aic<model2aic) {
            if (absaicdif<10) {
                tmessage="a";
            }
            else {
                tmessage="A";
            }
        }
        else {
            if (absaicdif<10) {
                tmessage="b";
            }
            else {
                tmessage="B";
            }
        }

        if (model1aicc<model2aicc) {
            if (absaiccdif<10) {
                tmessage+="a";
            }
            else {
                tmessage+="A";
            }
        }
        else {
            if (absaiccdif<10) {
                tmessage+="b";
            }
            else {
                tmessage+="B";
            }
        }

        if (plrtchi>0.05) {
            tmessage+="a";
        }
        else {
            if (plrtchi<0.01) {
                tmessage+="B";
            }
            else {
                tmessage+="b";
            }
        }

        if (repsnumber==0) {
            tmessage+="?";
        }
        else {
            if (parametricbootstrappingpvalue>0.05) {
                tmessage+="a";
            }
            else {
                if (parametricbootstrappingpvalue<0.01) {
                    tmessage+="B";
                }
                else {
                    tmessage+="b";
                }
            }
        }
        if (notquietmode) {
            PrintMessage();
        }
        summaryofresults+="\n";
        if (tablef_open) {
            tmessage+="\n";
            tablef<<tmessage;
        }

    }
	gsl_vector_free(tipscombresid);	
}
}
else {
    weighttotal-=treeweight;
    message="Tree ";
    message+=chosentree;
    message+=" excluded: has at least one subtree matrix that is singular\n";
    PrintMessage();
    summaryofresults+=chosentree;
    summaryofresults+="\t--Has a singular subtree VCV. No output.--\n";
    if (tablef_open) {
        tmessage="";
        tmessage+=chosentree;
        tmessage+="\t--Has a singular subtree VCV. No output.--\n";
        tablef<<tmessage;
    }
}
gsl_matrix_free(VCVcomb);
        }
chosentree=originalchosentree; //restore initial values.
chosenchar=originalchosenchar;
message=summaryofresults;
PrintMessage();
if (charloop==false && treeloop==true) {
    message="\n\nAverage result (weighted by tree weights)\n";
    if (noerror==false) {
        message+="\tError: Some trees had errors\n";
        PrintMessage();
    }
    else {
        map<nxsstring, int>::const_iterator iter5;
        int taxsetnumber=0;
        for (iter5=chosentaxsetntaxmap.begin();iter5!=chosentaxsetntaxmap.end();++iter5) {
            nxsstring taxset=iter5->first;
            taxsetnumber++;
            message+="\tTaxset ";
            message+=taxset;
            message+=": rate = ";
            message+=gsl_vector_get(weightedratevector,taxsetnumber)/weighttotal;
            message+=", anc state = ";
            message+=gsl_vector_get(weightedancstatevector,taxsetnumber)/weighttotal;
            message+="\n";
        }
        message+="\tSingle rate parameter model rate = ";
        message+=gsl_vector_get(weightedratevector,0)/weighttotal;
        message+="\n\n        Single rate parameter           Multiple rate parameters\nAIC =     ";
        char outputstring[14];
        sprintf(outputstring,"%14.6f",weightedAIC1/weighttotal);
        message+=outputstring;
        message+="                    ";
        sprintf(outputstring,"%14.6f",weightedAIC2/weighttotal);
        message+=outputstring;
        message+="\nAICc =    ";
        sprintf(outputstring,"%14.6f",weightedAICc1/weighttotal);
        message+=outputstring;
        message+="                    ";
        sprintf(outputstring,"%14.6f",weightedAICc2/weighttotal);
        message+=outputstring;
        message+="\n\nChi-square p value = ";
        sprintf(outputstring,"%14.6f",weightedchip/weighttotal);
        message+=outputstring;
        message+="\nParametric bootstrap p value = ";
        if (repsnumber>0) {
            sprintf(outputstring,"%14.6f",weightedparamp/weighttotal);
            message+=outputstring;
        }
        else {
            message+="NA: Not done";
        }
        PrintMessage();
    }
}
gsl_vector_free(weightedratevector);
gsl_vector_free(weightedancstatevector);
        }
if (tablef_open) {
    tablef.close();
}
    }



/**@method NumOpt
*Allows setting of parameters for numerical optimization
*/
void BROWNIE::NumOpt( NexusToken& token)
{
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="Usage: NumOpt [options...]\n\n";
            message+="This sets options for numerical optimization\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\n Iter        <integer>                            ";
            message+=maxiterations;
            message+="\n Toler       <double>                             ";
            char outputstring[10];
            sprintf(outputstring,"%1.9f",stoppingprecision);
            message+=outputstring;
            message+="\n RandStart   <integer>                            ";
            message+=randomstarts;
            message+="\n Seed        <integer>                            ";
            message+=gslseedtoprint;
            message+="\n StepSize    <double>                             ";
            sprintf(outputstring,"%1.9f",stepsize);
            message+=outputstring;
            message+="\n Detail      Yes|No                               ";
            if (detailedoutput) {
                message+="Yes";
            }
            else {
                message+="No";
            }
            message+="\n\nIter sets the maximum number of iterations of the Nelder-Mead simplex algorithm.\n\nToler sets the precision of the stopping criterion: what amount\nof change in the likelihood is considered small enough to count as zero change.\n\nRandStart sets the number of random starts to use.\n\nStepSize sets the NM step size.\n\nDetail specifies whether or not to have detailed output from numerical optimization";
            PrintMessage();
        }
        else if( token.Abbreviation("Iter") ) {
            nxsstring numbernexus = GetNumber(token);
            maxiterations=atoi( numbernexus.c_str() ); //convert to int
        }
        else if( token.Abbreviation("TOler") ) {
            nxsstring numbernexus = GetNumber(token);
            stoppingprecision=atof( numbernexus.c_str() );
        }
        else if( token.Abbreviation("Randstart") ) {
            nxsstring numbernexus = GetNumber(token);
            randomstarts=atoi( numbernexus.c_str() ); //convert to int
        }
        else if( token.Abbreviation("SEed") ) {
            nxsstring numbernexus = GetNumber(token);
            gslseedtoprint=atoi( numbernexus.c_str() );
            gslseed=gslseedtoprint;
            gsl_rng_set(r,gslseed); //convert to int
            message="Sorry, setting the seed doesn't quite work yet\n";
            PrintMessage();
        }
        else if( token.Abbreviation("Detail") ) {
            nxsstring yesnodetail=GetFileName(token);
            if (yesnodetail[0] == 'n') {
                detailedoutput=false;
            }
            else {
                detailedoutput=true;
            }
        }
        else if( token.Abbreviation("STepsize") ) {
            nxsstring numbernexus = GetNumber(token);
            stepsize=atof( numbernexus.c_str() );
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading NUMOPT command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }

}


/**@method Garland
*Spits out values to use to compare Brownie's performance with other programs
*/
void BROWNIE::HandleGarland( NexusToken& token)
{
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="This will perform various tests. It assumes you have set a taxset named ALL\n\n";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading GARLAND command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
    //Do something
    int ntax=0;
    nxsstring chosentaxset="ALL";
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }
IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntax++;
    }
    gsl_matrix* VCV=gsl_matrix_calloc(ntax,ntax);
    VCV=DeleteStem(GetVCV(chosentaxset));
    gsl_permutation * p = gsl_permutation_alloc (ntax);
    int signum;
    gsl_matrix * VCVLU=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix_memcpy (VCVLU, VCV);
    gsl_linalg_LU_decomp (VCVLU,p, &signum);
    int sig=signum;
    double determinant =  gsl_linalg_LU_det (VCVLU,sig);
    message="Determinant is ";
    message+=determinant;
    PrintMessage();
    gsl_matrix * VCVrescaled=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix_memcpy (VCVrescaled, VCV);
    gsl_matrix_scale (VCVrescaled,1/determinant);
    ofstream vcvf;
    vcvf.open( "vcv.txt" );
    for (int currentrow=0;currentrow<ntax;currentrow++) {
        for (int currentcol=0;currentcol<ntax;currentcol++) {
            vcvf<<gsl_matrix_get(VCVrescaled,currentrow,currentcol);
            vcvf<<"\t";
        }
        vcvf<<"\n";
    }
    vcvf.close();
    message+="The rescaled VCV matrix has been saved to vcv.txt";
    PrintMessage();
    gsl_vector * tips=gsl_vector_calloc(ntax);
    tips=GetTipValues(chosentaxset,chosenchar);
    double meantips=0;
    for (int taxonnum=0;taxonnum<ntax;taxonnum++) {
        meantips+=gsl_vector_get(tips,taxonnum);
    }
    meantips=meantips/ntax;
    double ancestralstate=GetAncestralState(VCVrescaled,tips);
    message="Mean tip value is ";
    message+=meantips;
    message+="\nReconstructed ancestral state is ";
    message+=ancestralstate;
    PrintMessage();
    gsl_vector *variance =gsl_vector_calloc(ntax);
    int originalmodel=chosenmodel;
    chosenmodel=1;




    OptimizationFn my_fnA(VCVrescaled,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);


    gsl_vector *optimalrateA=gsl_vector_calloc(2);
				//gsl_vector_memcpy(optimalrate,my_fn.OptimizeRateWithGivenTipVariance());
    gsl_vector_memcpy(optimalrateA,my_fnA.GeneralOptimization(1));
    message="\nOptimal rate = ";
    message+=gsl_vector_get(optimalrateA,0);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateA,1);
    PrintMessage();
    OptimizationFn my_fnB(VCVrescaled,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);
    chosenmodel=3;
    message="OU model chosen";
    PrintMessage();
    gsl_vector *optimalrateB=gsl_vector_calloc(6);
    gsl_vector_memcpy(optimalrateB,my_fnB.GeneralOptimization(3));
    message="\nOptimal rate = ";
    message+=gsl_vector_get(optimalrateB,0);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateB,3);
    message+="\nAncestral state = ";
    message+=gsl_vector_get(optimalrateB,1);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateB,4);
    message+="\nd = ";
    message+=gsl_vector_get(optimalrateB,2);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateB,5);
    chosenmodel=4;
    PrintMessage();
    message="ACDC model chosen";
    PrintMessage();
    OptimizationFn my_fnC(VCVrescaled,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);

    gsl_vector *optimalrateC=gsl_vector_calloc(6);
    gsl_vector_memcpy(optimalrateC,my_fnC.GeneralOptimization(4));
    message="\nOptimal rate = ";
    message+=gsl_vector_get(optimalrateC,0);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateC,3);
    message+="\nAncestral state = ";
    message+=gsl_vector_get(optimalrateC,1);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateC,4);
    message+="\ng = ";
    message+=gsl_vector_get(optimalrateC,2);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateC,5);
    PrintMessage();



    chosenmodel=originalmodel;
	gsl_matrix_free(VCV);
	gsl_matrix_free(VCVLU);
	gsl_matrix_free(VCVrescaled);
	gsl_vector_free(tips);
	gsl_vector_free(variance);
	gsl_vector_free(optimalrateA);
	gsl_vector_free(optimalrateB);
	gsl_vector_free(optimalrateC);
	gsl_permutation_free (p);
}


/**
* @method HandleTipVariance [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the TipVariance command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleTipVariance( NexusToken& token )
{
    bool starting = false;
    bool stopping = false;
    bool appending = false;
    bool replacing = false;
    bool name_provided = false;
    bool help_provided = false;

    // Retrieve all tokens for this command, stopping only in the event
    // of a semicolon or an unrecognized keyword
    //
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            help_provided = true;
            message="Usage: TipVariance [options...]\n\n";
            message+="Observations of taxon means have uncertainty. By default, the program assumes that these values are known exactly.\nIf you have measured these variances for each taxon (so that, for each measured character,\nthe data matrix has the mean followed by the variance), choose Given.\nIf you want estimate one tip variance across all the taxa, choose Same.\nIf you want to assume no tip variance (one fewer parameter), choose None\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\n Type        Given|Same|None                      ";
            if (tipvariancetype==0) {
                message+="None";
            }
            else if(tipvariancetype==1) {
                message+="Same";
            }
            else if(tipvariancetype==2) {
                message+="Given";
            }
            else {
                message+=tipvariancetype;
            }
            message+="\n\nNote that if the tip variance type is set to given, values in the matrix are\nassumed to consist of a taxon average followed by the\ncorresponding variance, so the matrix might be\n";
            message+=" taxon1   char1_mean  char1_variance  char2_mean  char2_variance  char3_mean  char3_variance...\n";
            message+=" taxon2   char1_mean  char1_variance  char2_mean  char2_variance  char3_mean  char3_variance...\n";
            message+=" taxon3   char1_mean  char1_variance  char2_mean  char2_variance  char3_mean  char3_variance...\n";
            message+="\nThis differs somewhat from the NEXUS specification (Maddison et al, 1997), which would use an\nITEMS=(AVERAGE VARIANCE) command and a somewhat different matrix.\n";
            PrintMessage();
        }
        else if( token.Abbreviation("Type") ) {
            nxsstring tipvarianceinput=GetFileName(token);
            if (tipvarianceinput[0] == 'N' || tipvarianceinput[0] == 'n') {
                tipvariancetype=0;
            }
            else if(tipvarianceinput[0] == 'S' || tipvarianceinput[0] == 's') {
                tipvariancetype=1;
            }
            else if(tipvarianceinput[0] == 'G' || tipvarianceinput[0] == 'g') {
                tipvariancetype=2;
            }
            else {
                errormsg = "Unexpected option (";
                errormsg += tipvarianceinput;
                errormsg += ") encountered reading TipVariance command";
                throw XNexus( errormsg);
            }
            message="\n You have chosen ";
            if (tipvariancetype==0) {
                message+="to assume that the tip variance is zero.";
            }
            else if(tipvariancetype==1) {
                message+="to estimate one tip variance parameter.";
            }
            else if(tipvariancetype==2) {
                message+="to use given variances from the data matrix.\nThis will assume that odd characters (char 1, 3, 5...) are taxon averages\nand even characters are the corresponding variances.\n";
            }
            else {
                message+=tipvariancetype;
            }
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading TipVariance command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}





/**
* @method HandleEcho [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the ECHO command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleEcho( NexusToken& token )
{
    bool starting = false;
    bool stopping = false;
    bool appending = false;
    bool replacing = false;
    bool name_provided = false;
    bool help_provided = false;

    // Retrieve all tokens for this command, stopping only in the event
    // of a semicolon or an unrecognized keyword
    //
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            help_provided = true;
            message="Usage: Echo [options...]\n\n";
            message+="Records everything you type into a Brownie-executable batch file,\nsuitable for pasting at the end of your data file\n(delete the #nexus if you do this)\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\nFile         <echo-file-name>                     ";
            if (echof_open) {
                message+=echofname;
            }
            else {
                message+="[no file open]";
            }
            message+="\nReplace      No|Yes                              *No";
            message+="\nAppend       No|Yes                              *No";
            message+="\n          Start|Stop";
            message+="\n                                                 *Option is nonpersistent\n\nNote that you should specify =yes or =no for Replace \nor Append, but not for Start or Stop.\nFor example,\n\echo start file=test.nex replace=yes;\n\n";
            message+="\n                                                 *Option is nonpersistent\n\n";
            PrintMessage();
        }
        else if( token.Abbreviation("STOp") ) {
            stopping=true;
        }
        else if( token.Abbreviation("STArt") ) {
            starting=true;
        }
        else if( token.Abbreviation("Replace") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n') {
                replacing=false;
            }
            else {
                replacing=true;
            }
        }
        else if( token.Abbreviation("Append") ) {
            nxsstring yesnoappend=GetFileName(token);
            if (yesnoappend[0] == 'n') {
                appending=false;
            }
            else {
                appending=true;
            }
        }
        else if( token.Abbreviation("File") ) {
            echofname = GetFileName(token);
            name_provided = true;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading ECHO command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }

    // Check for incompatible combinations of keywords
    //
    if( stopping && ( starting || appending || replacing || name_provided ) ) {
        errormsg = "Cannot specify STOP with any of the following START, APPEND, REPLACE, FILE";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( appending && replacing ) {
        errormsg = "Cannot specify APPEND and REPLACE at the same time";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( echof_open && ( starting || name_provided || appending || replacing ) ) {
        errormsg = "Cannot start echo file since echo file is already open";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    // Is user closing an open echo file?
    //
    if( stopping ) {
        echof<<"end;"<<endl;
        echof.close();
        echof_open = false;

        message = "\nEcho file closed";
        PrintMessage();

        return;
    }

    // If this far, must be attempting to open a echo file
    //
    if( !name_provided && !help_provided) {
        errormsg = "Must provide a file name when opening a echo file\n";
        errormsg += "e.g., echo file=batch.nex start replace;";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( appending ) {
        echof_open = true;
        echof.open( echofname.c_str(), ios::out | ios::app );

        message = "\nAppending to echo file ";
        message += echofname;
        PrintMessage();
        echof<<"Begin Brownie;"<<endl;
    }
    else if( replacing ) {
        echof_open = true;
        echof.open( echofname.c_str() );

        message = "\nReplacing echo file ";
        message += echofname;
        PrintMessage();
        echof<<"#nexus"<<endl<<"Begin Brownie;"<<endl;
    }
    else  {
        bool exists = FileExists( echofname.c_str() );
        bool userok = true;
        if( exists && !UserSaysOk( "Ok to replace?", "Echo file specified already exists" ) )
            userok = false;
        if( userok ) {
            echof_open = true;
            echof.open( echofname.c_str() );
            echof<<"#nexus"<<endl<<"Begin Brownie;"<<endl;
        }

        if( exists && userok ) {
            message = "\nReplacing echo file ";
            message += echofname;
        }
        else if( userok ) {
            message = "\nEcho file ";
            message += echofname;
            message += " opened";
        }
        else {
            message = "\nEcho command aborted";
        }
        PrintMessage();
    }
}




/**
* @method HandleLog [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the LOG command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleLog( NexusToken& token )
{
    bool starting = false;
    bool stopping = false;
    bool appending = false;
    bool replacing = false;
    bool name_provided = false;
    bool help_provided = false;


    // Retrieve all tokens for this command, stopping only in the event
    // of a semicolon or an unrecognized keyword
    //
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            help_provided=true;
            message="Usage: Log [options...]\n\n";
            message+="Records output into a log file.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\nFile         <log-file-name>                      ";
            if (logf_open) {
                message+=logfname;
            }
            else {
                message+="[no file open]";
            }
            message+="\nReplace      No|Yes                              *No";
            message+="\nAppend       No|Yes                              *No";
            message+="\n          Start|Stop";
            message+="\n                                                 *Option is nonpersistent\n\nNote that you should specify =yes or =no for Replace \nor Append, but not for Start or Stop.\nFor example,\n\tlog start file=test.log replace=yes;\n\n";
            PrintMessage();
        }
        else if( token.Abbreviation("STOp") ) {
            stopping=true;
        }
        else if( token.Abbreviation("STArt") ) {
            starting=true;
        }
        else if( token.Abbreviation("Replace") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n') {
                replacing=false;
            }
            else {
                replacing=true;
            }
        }
        else if( token.Abbreviation("Append") ) {
            nxsstring yesnoappend=GetFileName(token);
            if (yesnoappend[0] == 'n') {
                appending=false;
            }
            else {
                appending=true;
            }
        }
        else if( token.Abbreviation("File") ) {
            logfname = GetFileName(token);
            name_provided = true;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading LOG command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }

    // Check for incompatible combinations of keywords
    //
    if( stopping && ( starting || appending || replacing || name_provided ) ) {
        errormsg = "Cannot specify STOP with any of the following START, APPEND, REPLACE, FILE";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( appending && replacing ) {
        errormsg = "Cannot specify APPEND and REPLACE at the same time";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( logf_open && ( starting || name_provided || appending || replacing ) ) {
        errormsg = "Cannot start log file since log file is already open";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    // Is user closing an open log file?
    //
    if( stopping ) {
        logf.close();
        logf_open = false;

        message = "\nLog file closed";
        PrintMessage();

        return;
    }

    // If this far, must be attempting to open a log file
    //
    if( !name_provided && !help_provided) {
        errormsg = "Must provide a file name when opening a log file\n";
        errormsg += "e.g., log file=doofus.txt start replace;";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( appending ) {
        logf_open = true;
        logf.open( logfname.c_str(), ios::out | ios::app );

        message = "\nAppending to log file ";
        message += logfname;
        PrintMessage();
    }
    else if( replacing ) {
        logf_open = true;
        logf.open( logfname.c_str() );

        message = "\nReplacing log file ";
        message += logfname;
        PrintMessage();
    }
    else  {
        bool exists = FileExists( logfname.c_str() );
        bool userok = true;
        if( exists && !UserSaysOk( "Ok to replace?", "Log file specified already exists" ) )
            userok = false;
        if( userok ) {
            logf_open = true;
            logf.open( logfname.c_str() );
        }

        if( exists && userok ) {
            message = "\nReplacing log file ";
            message += logfname;
        }
        else if( userok ) {
            message = "\nLog file ";
            message += logfname;
            message += " opened";
        }
        else {
            message = "\nLog command aborted";
        }
        PrintMessage();
    }
}

/**
* @method HandleNextCommand [void:public]
 *
 * Accepts a string in the form of a BROWNIE block containing one command
 * and processes it just like a real BROWNIE block in a NEXUS data file.
 */
void BROWNIE::HandleNextCommand()
{
std::istrstream cmdin( next_command );
    //std::istringstream cmdin( next_command );
    //cerr << "next_command is: " << next_command <<endl;
    //istream cmdin( next_command );

    NexusToken token(cmdin);
    try {
        Read( token );
    }
    catch( XNexus x )
    {
        NexusError( errormsg, x.pos, x.line, x.col );
        Reset();
    }
}

/**
* @method NexusError [virtual void:public]
 * @param msg [nxsstring&] the error message
 * @param pos [streampos] the point in the NEXUS file where the error occurred
 * @param line [long] the line in the NEXUS file where the error occurred
 * @param col [long] the column in the NEXUS file where the error occurred
 *
 * Called when an error is encountered in a NEXUS file. Allows program to
 * give user details of the error as well as the precise location of the
 * error. Virtual function that overrides the pure virtual function in the
 * base class Nexus.
 */
void BROWNIE::NexusError( nxsstring& msg, streampos /* pos */, long line, long col )
{
    message = "\n";
    message += msg;
    PrintMessage();

    if( inf_open )
    {
        message = "Line:   ";
        message += line;
        PrintMessage();

        message = "Column: ";
        message += col;
        PrintMessage();
    }
}

/**
* @method PreprocessNextCommand [void:public]
 *
 * Begins with the command just entered by the user, which is stored in
 * the data member next_command, adds a semicolon (if the user failed
                                                   * to supply one), and then adds "end;" so the whole bundle looks
 * like a very short BROWNIE block.  This is then passed to HandleNextCommand,
 * which processes it just like a real BROWNIE block in a NEXUS data file.
 */
void BROWNIE::PreprocessNextCommand()
{
    // If user failed to add the terminating semicolon,
    // we'll do it now. We will also remove the line feed
    // at the end and add the command "end;" to the end
    // of the line (see explanation below).
    //
    int len = strlen(next_command);
    //assert( len > 0 );

    // Remove any whitespace characters from end of string entered by user
    //
    int i = len;
    while( i > 0 && next_command[i-1] == ' ' || next_command[i-1] == '\t' || next_command[i-1] == '\n' )
        i--;


    // If character at position i-1 is a semicolon, put '\0' terminator at position i;
    // otherwise, put a semicolon at position i and terminator at i+1
    //
    if( next_command[i-1] != ';' ) {
        next_command[i] = ';';
        i++;
    }
    assert( i <= COMMAND_MAXLEN );
    next_command[i] = '\0';

    if( echof_open ) {
        echof<<next_command<<endl;
    }

    // Now add a semicolon at the beginning and terminate with an "END;" command
    // so that we can pretend this is simply a very short private NEXUS block
    // containing only one command.  This allows us to simply use the Read
    // function we inherited from the base class BstBase to process the command.
    //
    len = strlen(next_command);
    assert( len < COMMAND_MAXLEN-2 );
    nxsstring tmp = ";";
    tmp += next_command;
    tmp += "end;";
    strcpy( next_command, tmp.c_str() );
}

/**
* @method PrintMessage [void:public]
 * @param linefeed [bool] if true, places newline character after message
 *
 * All output handled here.  Writes string currently stored in message
 * (a nxsstring data member) to the output file stream, if open, and also
 * to the console via cerr. Places newline after string if linefeed is true.
 */
void BROWNIE::PrintMessage( bool linefeed /* = true */ )
{
    cerr << message;
    if( linefeed )
        cerr << endl;

    if( logf_open ) {
        logf << message;
        if( linefeed )
            logf << endl;
    }
}

/**
* @method PurgeBlocks [void:protected]
 *
 * Detaches all blocks, deletes them, creates new blocks, and
 * finally adds the new blocks.
 */
void BROWNIE::PurgeBlocks()
{
    Detach( taxa );
    Detach( trees );
    Detach( assumptions );
    Detach( characters );

    delete characters;
    delete assumptions;
    delete trees;
    delete taxa;

    taxa = new TaxaBlock();
    trees = new TreesBlock(*taxa);
    assumptions = new AssumptionsBlock( *taxa );
    characters = new CharactersBlock( *taxa, *assumptions );

    Add( taxa );
    Add( trees );
    Add( assumptions );
    Add( characters );
}

void BROWNIE::Assign(NexusToken& token)
{
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            message="Current assignment vector is:\n";
            for (int i=0; i<convertsamplestospecies.size(); i++) {
                message+=convertsamplestospecies[i];
                message+=" ";
            }
            PrintMessage();
            break;
        }
		if( token.Abbreviation("?") ) {
            message="Usage: Assign [manual | onetoone | automatic]";
			PrintMessage();
			
        }
        else if (token.Abbreviation("Manual")) {
            int ntax=taxa->GetNumTaxonLabels();
            message="Number of taxa is ";
            message+=ntax;
            PrintMessage();
            convertsamplestospecies.clear();
            for (int taxon=0; taxon<ntax; taxon++) {
                int assignment;
                while ((cout<<"Taxon "<<taxa->GetTaxonLabel(taxon)<<" Assign to species# = ") && (!(cin >> assignment) || assignment < 1 || assignment > ntax)) {
                    cout << "Incorrect input: enter a number between 1 and "<<ntax<<": ";
                    cin.clear();
                    //cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }
                convertsamplestospecies.push_back(assignment);
            }
            bool changedmaxnum=false;
            for (int i=0;i<convertsamplestospecies.size();i++) {
                if (convertsamplestospecies[i]>maxnumspecies) {
                    maxnumspecies=convertsamplestospecies[i];
                    changedmaxnum=true;
                }
            }
            if (changedmaxnum) {
                maxnumspecies++;
                message="Increased the maximum number of species allowed to ";
                message+=maxnumspecies;
                PrintMessage();
            }
        }
        else if (token.Abbreviation("Onetoone")) {
            message="Now assigning each sample to its own species";
            PrintMessage();
            convertsamplestospecies.clear();
            int ntax=taxa->GetNumTaxonLabels();
            for (int i=0; i<ntax; i++) {
                convertsamplestospecies.push_back(i+1);
            }
            bool changedmaxnum=false;
            for (int i=0;i<convertsamplestospecies.size();i++) {
                if (convertsamplestospecies[i]>maxnumspecies) {
                    maxnumspecies=convertsamplestospecies[i];
                    changedmaxnum=true;
                }
            }
            if (changedmaxnum) {
                maxnumspecies++;
                message="Increased the maximum number of species allowed to ";
                message+=maxnumspecies;
                PrintMessage();
            }
        }
        else if (token.Abbreviation("Automatic")) {
            message="Now using the first two symbols in the taxon names to assign species";
            PrintMessage();
            map<nxsstring, int> assignmentmap;
            convertsamplestospecies.clear();
            int ntax=taxa->GetNumTaxonLabels();
            int nspecies=0;
            for (int taxon=0; taxon<ntax; taxon++) {
                nxsstring maplabel=taxa->GetTaxonLabel(taxon);
                message="Taxon ";
                message+=maplabel;
                message+=" abbreviated to ";
                nxsstring shortenedlabel;
                char* s =new char[3];
                strncpy( s, maplabel.c_str(), 2 );
                s[2]='\0';
                shortenedlabel = s;
                delete(s);
                message+=shortenedlabel;
                //message+=" or ";
                //message+=maplabel.ShortenTo(5);
                message+=" assigned to species: taxon";

                if (assignmentmap.count(shortenedlabel)==0) { //No match
                    nspecies++;
                    assignmentmap[shortenedlabel]=nspecies;
                    convertsamplestospecies.push_back(nspecies);
                    message+=nspecies;
                    PrintMessage();
                }
                else { //Already in there
                    typedef map<nxsstring, int>::const_iterator ITER;
                    ITER found=assignmentmap.find(shortenedlabel);
                    int assignedtospeciesnum=found->second;
                    convertsamplestospecies.push_back(assignedtospeciesnum);
                    message+=assignedtospeciesnum;
                    PrintMessage();
                }
            }
            bool changedmaxnum=false;
            for (int i=0;i<convertsamplestospecies.size();i++) {
                if (convertsamplestospecies[i]>maxnumspecies) {
                    maxnumspecies=convertsamplestospecies[i];
                    changedmaxnum=true;
                }
            }
            if (changedmaxnum) {
                maxnumspecies++;
                message="Increased the maximum number of species allowed to ";
                message+=maxnumspecies;
                PrintMessage();
            }

        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Assign command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }

}

char* BROWNIE::OutputForGTP(ContainingTree *SpeciesTreePtr)
{
    nxsstring TreeFile="#nexus\nbegin trees;";
    //first the species tree
    //SpeciesTree.RandomTree(taxa->GetNumTaxonLabels());
    TreeFile+="\ntree speciestree = ";
    TreeFile+=PipeSpeciesTree(SpeciesTreePtr);
    //now the gene trees
    int originalchosentree=chosentree;
    for (chosentree = 0; chosentree < trees->GetNumTrees(); chosentree++) {
		      Tree t = intrees.GetIthTree(chosentree);
		      TreeFile+="\ntree ";
                      //  if (t.GetName() != "")
                      //		  TreeFile+=NEXUSString (t.GetName());
                      //else
                      //{
                      TreeFile+="genetree_";
                          TreeFile+=chosentree+1;
                          //}
                          TreeFile+=" = ";
                          // Tree
                          TreeFile+=PipeGTP(t);
    }
    chosentree=originalchosentree;
    TreeFile+="\nend;";
    char* OutputForGTPChar=(char*)TreeFile.c_str();
    // if (debugmode) {
    // cout<<"OutputForGTP output =\n"<<OutputForGTPChar<<endl<<endl;
    // }
    return OutputForGTPChar;
}

//Take the current best trees vector, format each tree and copy it into a new vector
vector<ContainingTree> BROWNIE::MakePrettyForOutput()
{
    vector<ContainingTree> PrettyTrees;
    for (int i=0; i<RawBestTrees.size(); i++) {
        ContainingTree PrettyTree=RawBestTrees[i];
        PrettyTree.FindAndSetRoot();
        PrettyTree.Update();
        PrettyTree.GetNodeDepths();
        NodeIterator <Node> n (PrettyTree.GetRoot());
        cur = n.begin();
        while (cur) {
            if (cur->IsLeaf()) {
                nxsstring newlabel=PipeLeafFinalSpeciesTree();
                cur->SetLabel(newlabel);
            }
            else {
                cur->SetLabel("");
            }
            cur = n.next();
        }
        PrettyTrees.push_back(PrettyTree);
    }
    return PrettyTrees;
}

void BROWNIE::FormatAndStoreBestTree(ContainingTree *NewBestTree,vector<double> scorevector)
{
    ContainingTree FormattedNewBestTree=*NewBestTree;
    FormattedNewBestTree.FindAndSetRoot();
    FormattedNewBestTree.Update();
    FormattedNewBestTree.GetNodeDepths();
    NodeIterator <Node> n (FormattedNewBestTree.GetRoot());
    cur = n.begin();
    while (cur) {
        if (cur->IsLeaf()) {
            nxsstring newlabel=PipeLeafFinalSpeciesTree();
            cur->SetLabel(newlabel);
        }
        else {
            cur->SetLabel("");
        }
        cur = n.next();
    }
    bool newtree=true;
    for (int j=0; j<FormattedBestTrees.size(); j++) {
        ContainingTree t1=FormattedNewBestTree;
        ContainingTree t2=FormattedBestTrees[j];
        int ntaxt1=t1.GetNumLeaves();
        int ntaxt2=t2.GetNumLeaves();
        int ntaxincommon=PrepareTreesForTriplet(&t1,&t2);
        if (ntaxincommon==ntaxt1 && ntaxincommon==ntaxt2) {
            if (ntaxincommon>=3) { //so, all leaves are in common, but is the topology the same? If there are only two or one species, yes, so newtree=false. If there are at least 3 taxa, there are different possible topologies, so maybe, thus check below
                vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,ntaxincommon);
                int numberdisagree=tripletoverlapoutput[1];
                if (numberdisagree==0) {
                    newtree=false;
                }
            }
            else {
                newtree=false; //all leaves in common, must have same topology, so same tree
            }
        }
    }
    if (newtree) {
        FormattedBestTrees.push_back(FormattedNewBestTree);
		BestConversions.push_back(convertsamplestospecies);
		TotalScores.push_back(scorevector[0]);
		GTPScores.push_back(scorevector[1]);
		StructScores.push_back(scorevector[2]);
    }
    for (int k=0; k<FormattedBestTrees.size(); k++) {
        (FormattedBestTrees[k]).Update();
        (FormattedBestTrees[k]).GetNodeDepths();
        (FormattedBestTrees[k]).ClearInternalLabels();
    }
    if (newtree) {
		if (jackknifesearch) {
				jackknifetreestooutput="";
		}
        ofstream outtreef;
        outtreef.open(treefilename.c_str());
        outtreef<<"#nexus\nbegin trees;\n";
        //outtreef<<"[heuristic search, best results after replicate "<<replicate<<"\nSearch options: \n]\n";
		
        for (int i=0; i<FormattedBestTrees.size(); i++) {
			outtreef<<"[ assignment vector: (";
			int numspecies=0;
			for (int j=0;j<convertsamplestospecies.size();j++) {
				outtreef<<" "<<BestConversions[i][j];
				numspecies=GSL_MAX(numspecies,BestConversions[i][j]);
			}
			outtreef<<" ) ]\n";
			if (jackknifesearch) {
				jackknifetreestooutput+="tree jackrep";
				jackknifetreestooutput+=jackrep;
				jackknifetreestooutput+=".";
				jackknifetreestooutput+=i+1;
				jackknifetreestooutput+=" = [&W 1/";
				int ntrees=FormattedBestTrees.size();
				jackknifetreestooutput+=ntrees;
				jackknifetreestooutput+="] [&R] ";
				jackknifetreestooutput+=ReturnFinalSpeciesTree(FormattedBestTrees[i]);
				jackknifetreestooutput+="\n";
			}
            outtreef<<"tree sptre"<<i+1<<" = [ "<<"Number of species: "<<numspecies<<"; Score: "<<TotalScores[i]<<" = "<<(1-structwt)<<" x GTP ("<<GTPScores[i]<<") + "<<structwt<<" x Struct ("<<StructScores[i]<<") ] [&R] ";
            outtreef<<ReturnFinalSpeciesTree(FormattedBestTrees[i]);
            outtreef<<endl;
        }
        outtreef<<"end;";
        outtreef.close();
		
		//output for mesquite
		nxsstring mesquitefilename=treefilename;
		mesquitefilename+=".mesquite.nex";
		ofstream outmesqf;
		outmesqf.open(mesquitefilename.c_str());
		outmesqf<<"#nexus\n";
		
		outmesqf<<"\nbegin taxa;\ntitle SamplesBlock;\ndimensions ntax="<<taxa->GetNumTaxonLabels()<<";\ntaxlabels\n\t"; //only need one samples block and one samples tree block
		for (int k=0;k<taxa->GetNumTaxonLabels();k++) {
			outmesqf<<taxa->GetTaxonLabel(k)<<" ";
		}
		outmesqf<<"\n;\nend;\n";

		outmesqf<<"\nbegin trees;\nTitle 'Trees of individual loci';\nLINK Taxa=SamplesBlock;\n";
		//output gene trees
		for (int selectedtree = 0; selectedtree < intrees.GetNumTrees(); selectedtree++) {
			outmesqf<<"tree genetree"<<selectedtree+1<<" = ";
			Tree t=intrees.GetIthTree(selectedtree);
			t.Write(outmesqf);
			outmesqf<<endl;
		}
		//outmesqf<<"OUTPUT GENE TREES\n";
		outmesqf<<"end;\n";
		
		for (int i=0; i<FormattedBestTrees.size(); i++) {
			int numspecies=0;
			for (int j=0;j<convertsamplestospecies.size();j++) {
				numspecies=GSL_MAX(numspecies,BestConversions[i][j]);
			}
			outmesqf<<"\nbegin taxa;\ntitle SpeciesBlock"<<i<<";\ndimensions ntax="<<numspecies<<";\ntaxlabels\n\t";
			for (int k=1;k<=numspecies;k++) {
					outmesqf<<"taxon"<<k<<" ";
			}
			outmesqf<<"\n;\nend;\n";
			
						
			outmesqf<<"\nbegin trees;\ntitle 'Species tree block "<<i<<"';\nLINK Taxa=SpeciesBlock"<<i<<";\n";
			//output tree
			outmesqf<<"tree sptree"<<i+1<<" = ";
			RawBestTrees[i].Write(outmesqf);
			outmesqf<<"\n";
			outmesqf<<"end;\n";
			
			outmesqf<<"\nbegin taxaassociation;\ntitle Assignment"<<i<<";\nTaxa SpeciesBlock"<<i<<" , SamplesBlock;\nAssociates\n\n";
			for (int k=1;k<=numspecies;k++) {
				outmesqf<<"taxon"<<k<<" / ";
				for(int l=0;l<convertsamplestospecies.size();l++) {
					if (convertsamplestospecies[l]==k) {
						outmesqf<<GetTaxonLabel(l)<<" ";
					}
				}
				if (k<numspecies) {
					outmesqf<<",";
				}
				outmesqf<<"\n";
			}
			outmesqf<<";\nend;\n";
			
        }
		
    }
}



double BROWNIE::GetGTPScoreNew(ContainingTree *SpeciesTreePtr)
{
    //cout<<"SpeciesTree"<<endl;
    //(*SpeciesTreePtr).Write(cout);
    //cout<<endl;
    //uses algorithm from Zmasek & Eddy 2001, with modification by Sanderson for inferring only "strong" duplications (those that haven't happened after the last speciation event)
    SpeciesTreePtr->SetLeafNumbers();
    //cout<<"has set leafnumbers\n";
    PreorderIterator <Node> n (SpeciesTreePtr->GetRoot());
    //cout<<"set preorderiterator\n";
    //cout<<"Species Tree = "<<endl;
    // SpeciesTreePtr->ReportTreeHealth();
    //FormatAndStoreBestTree(SpeciesTreePtr);
    //SpeciesTreePtr->Update();
    //SpeciesTreePtr->Draw(cout);
    //cout<<ReturnFinalSpeciesTree(FormattedBestTrees.back());
    //FormattedBestTrees.pop_back();
    //cout<<endl;


    NodePtr currentnode = n.begin();
    //cout<<"    NodePtr currentnode = n.begin();\n";
    assert(currentnode!=NULL);
    //cout<<"beginning to update\n";
    SpeciesTreePtr->Update();
    //cout<<"finished updating\n";
    //cout<<endl;
    //SpeciesTreePtr->ReportTreeHealth();
    //cout<<endl<<"finished reporttreehealth"<<endl;
    int labelcount=1;
    while (currentnode)
    {
        //	cout<<"in while(currentnode)\n";
        //	cout<<"currentnode="<<currentnode<<endl;
        //	cout<<"currentnode anc = "<<currentnode->GetAnc()<<endl;
        //	cout<<"current node index = "<<currentnode->GetIndex()<<endl;
        currentnode->SetIndex(labelcount);
        //	cout<<"currentnode->SetIndex(labelcount)\n";
        labelcount++;
        //	cout<<"labelcount++\n";
        currentnode = n.next();
        //	cout<<"currentnode = n.next();\n";
    }
    //cout<<"out of while(currentnode)\n";
    int originalchosentree=chosentree;
    double weightednumDup=0;
    for (int selectedtree = 0; selectedtree < intrees.GetNumTrees(); selectedtree++) {
		if (gtptoohigh) {
			break;
		}
		bool usethistree=true;
		if (jackknifesearch) {
			if (jackknifevector[selectedtree]==0) {
				usethistree=false;
			}
		}
		if (usethistree) {
			int numDup=0;
			Tree t=intrees.GetIthTree(selectedtree);
			
        //t.Write(cout);
        //   cout<<endl;
        //cout<<"Gene Tree = "<<endl;
			t.Update();
        // t.Draw(cout);
        //  cout<<endl;
			NodeIterator <Node> m (t.GetRoot());
			NodePtr genetreenode=m.begin();
			NodePtr a;
			NodePtr b;
			NodePtr g1;
			NodePtr g2;
			NodePtr Mg1;
			NodePtr Mg2;
			map<NodePtr,NodePtr> GtoSmap;
			map<NodePtr,NodePtr>::iterator GtoSmapiter;
        //map<NodePtr,map<int, int> > setsM;
        //   map<int,int> g1map;
        //   map<int,int> g2map;
        //   map<int,int> a;
        // map<int,int> b;
        //map<NodePtr,map<int, int> >::iterator setsMiterator;
			while (genetreenode!=NULL) {
				if (gtptoohigh) {
					break;
				}				
				genetreenode->SetMarked(false); //in this context, whether it predates a speciation event and so could be a strong duplication
                                            //  cout<<"\ngenetreenode = "<<genetreenode<<" is leaf? "<<genetreenode->IsLeaf()<<" "<<genetreenode->GetLabel()<<endl;
				if (!(genetreenode->IsLeaf())) {
					g1=genetreenode->GetChild();
					g2=g1->GetSibling();
					if(g1->IsLeaf()) {
						int SampleNumber=taxa->FindTaxon(g1->GetLabel());
                    //cout<<"samplenum = "<<SampleNumber<<" convertsamplestospecies.size() = "<<convertsamplestospecies.size()<<endl;
                    //   tempmap.insert(make_pair(convertsamplestospecies[SampleNumber],1));
                    //  setsM.insert(make_pair(g1,tempmap));
						if (SampleNumber>=convertsamplestospecies.size()) {
							cout<<"SampleNumber "<<SampleNumber<<" but...";
							cout<<" convertsamplestospecies.size() = "<<convertsamplestospecies.size()<<endl;
							for(int i=0;i<convertsamplestospecies.size();i++) {
								cout<<convertsamplestospecies[i]<<" ";
							}
							cout<<endl;
						}
						assert(SampleNumber<convertsamplestospecies.size());
						GtoSmap.insert(make_pair(g1,SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber])));
                    //  cout<<"Found leaf "<<SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber]);
                    //  cout<<" "<<(SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber]))->GetLabel()<<endl;
					}
					if(g2->IsLeaf()) {
						int SampleNumber=taxa->FindTaxon(g2->GetLabel());
						GtoSmap.insert(make_pair(g2,SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber])));
                    //  cout<<"Found leaf "<<SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber]);
                    //  cout<<" "<<(SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber]))->GetLabel()<<endl;
						
					}
					if(g1->IsMarked() || g2->IsMarked()) {
						genetreenode->SetMarked(true); // it predates a speciation event, so can be a strong dup (idea from Sanderson)
					}
					GtoSmapiter=GtoSmap.find(g1);
					a=(*GtoSmapiter).second;
					Mg1=a;
					GtoSmapiter=GtoSmap.find(g2);
					b=(*GtoSmapiter).second;
					Mg2=b;
                //cout<<"a = "<<a<<" b = "<<b<<endl;
					while (a->GetIndex() != b->GetIndex()) {
						if (a->GetIndex() > b->GetIndex()) {
							a=a->GetAnc();
						}
						else {
							b=b->GetAnc();
						}
					}
					GtoSmap.insert(make_pair(genetreenode,a));
                //  cout<<"Mgenetreenode,Mg1,Mg2,genetreenode,g1,g2\n"<<a<<"\t"<<Mg1<<"\t"<<Mg2<<"\t"<<genetreenode<<"\t"<<g1<<"\t"<<g2<<endl;
					if((a==Mg1) || (a==Mg2)) {
                    //   cout<<"Dup";
						if(genetreenode->IsMarked()) {
                        //    cout<<"--Strong, will count"<<endl;
							numDup++;
                        //    cout<<"Match, numDup = "<<numDup<<endl;
						}
					}
					else {
						genetreenode->SetMarked(true); //since it's a speciation event
                                                   //  cout<<"SpeciationEvent"<<endl;
					}
				}
				genetreenode=m.next();
			}
        //     cout<<"numDup="<<numDup<<" weight = "<<trees->GetTreeWeight(chosentree)<<endl;
			weightednumDup+=(trees->GetTreeWeight(chosentree))*(1.0*numDup);
			if (weightednumDup*(1.0-structwt)>bestscorelocal) {
				weightednumDup=(0.0001+bestscorelocal)/(1.0-structwt);
				gtptoohigh=true;
				break;
			}
		}
    }
    chosentree=originalchosentree;
    //cout<<"weightednumDup="<<weightednumDup<<endl;
    return weightednumDup;
}

//took this function out as no longer depend on external gtp
//double BROWNIE::GetGTPScore(ContainingTree *SpeciesTreePtr)
//{
//    double totalscore=0;
//    badgtpcount=0;
//    if (SpeciesTreePtr->GetNumLeaves()>1) { //don't bother getting the GTP score if there's only one species
//        int originalchosentree=chosentree;
//        if (outputallatonce) {
//            totalscore=GSL_NEGINF;
//            badgtpcount=-1;
//            while (totalscore<0) {
//                totalscore=ReturnScore(OutputForGTP(SpeciesTreePtr),unrooted);
//                badgtpcount++;
//                if (badgtpcount==500) {
//                    cout<<"Potential error getting tree score, trying to recover...";
//                }
//                if (badgtpcount>10000) {
//                    cout<<"Has failed, aborting. Email bcomeara@ucdavis.edu with this message, also include the info below:\n\n";
//                    cout<<OutputForGTP(SpeciesTreePtr)<<endl;
//                    errormsg="Aborted";
//                    throw XNexus( errormsg);
//                }
//            }
//            if (badgtpcount>=500) {
//                cout<<"... fixed! Required "<<badgtpcount<<" tries\n";
//            }
//        }
//       else {
//            for (chosentree = 0; chosentree < trees->GetNumTrees(); chosentree++) {
//                double thisscore=GSL_NEGINF;
//                while (thisscore<0) {
//                    nxsstring TreeFile="#nexus\nbegin trees;";
//                    TreeFile+="\ntree speciestree = ";
//                    TreeFile+=PipeSpeciesTree(SpeciesTreePtr);
//                    Tree t = intrees.GetIthTree(chosentree);
//                    TreeFile+="\ntree ";
//                    TreeFile+="genetree_";
//                    TreeFile+=chosentree+1;
//                    TreeFile+=" = ";
//                    TreeFile+=PipeGTP(t);
//                    TreeFile+="\nend;";
//                    thisscore=ReturnScore((char*)TreeFile.c_str(),unrooted);
//                    if (thisscore<0) {
//                        badgtpcount++;
//                    }
//                 }
//                totalscore+=(trees->GetTreeWeight(chosentree))*thisscore; //Get properly weighted score
//            }
//        }
//        chosentree=originalchosentree;
//    }
//return totalscore;
//}


/**
* @method OutputComment [virtual void:public]
 * @param msg [nxsstring&] the output comment to be displayed
 *
 * This function is called whenever an output comment (i.e., a comment
                                                       * beginning with an exclamation point) is found in the data file.
 */
void BROWNIE::OutputComment( nxsstring& msg )
{
    message=msg;
    PrintMessage();
}


/**
* @method Read [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * This function provides the ability to read everything following
 * the block name (which is read by the Nexus object) to the end or
 * endblock statement. Characters are read from the input stream
 * in. Overrides the pure virtual function in the base class.
 */
void BROWNIE::Read( NexusToken& token )
{
    isEmpty = false;

    // this should be the semicolon after the block name
    //
    token.GetNextToken();

    if( !token.Equals(";") ) {
        errormsg = "Expecting ';' after ";
        errormsg += id;
        errormsg += " block name, but found ";
        errormsg += token.GetToken();
        errormsg += " instead";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    for(;;)
    {
        token.GetNextToken();

        if( token.Abbreviation("ENdblock") ) {
            HandleEndblock( token );
            break;
        }
        else if( token.Abbreviation("Help") ) {
            HandleHelp( token );
        }
        else if( token.Abbreviation("Log") ) {
            HandleLog( token );
        }
        else if( token.Abbreviation("ECho") ) {
            HandleEcho( token );
        }
        else if( token.Abbreviation("EXecute") ) {
            HandleExecute( token );
        }
        else if( token.Abbreviation("Gettrees") ) {
            HandleGettrees( token );
        }
        else if( token.Abbreviation("CHoose") ) {
            HandleChoose( token );
        }
        else if( token.Abbreviation("SEt") ) {
            HandleSet( token );
        }
        else if( token.Abbreviation("COmpare") ) {
            HandleCompareRandomTrees( token );
        }
        else if( token.Abbreviation("HSearch") ) {
            HandleHeuristicSearch( token );
        }
        else if( token.Abbreviation("HEUristicsearch") ) {
            HandleHeuristicSearch( token );
        }
		else if( token.Abbreviation("JAckknife") ) {
			jackknifesearch=true;
			HandleHeuristicSearch( token ) ;  //a jackknife search is just like a heuristic search, just with jackknife options
		}
		else if( token.Abbreviation("EXhaustivesearch") ) {
			DoExhaustiveSearch();
		}
		else if( token.Abbreviation("NAst") ) {
			HandleNast(token);
		}
        else if( token.Abbreviation("ASsign") ) {
            Assign( token );
        }
        else if( token.Abbreviation("CItation") ) {
            HandleCitation( token );
        }
		else if( token.Abbreviation("ACcuracy") ) {
            HandleAccuracy( token );
        }
        else if( token.Abbreviation("Showtree") ) {
            HandleShowtree( token );
        }
        else if( token.Abbreviation("Blocks") ) {
            HandleBlocks( token );
        }
        else if( token.Abbreviation("VCV") ) {
            HandleVCV( token );
        }
        else if( token.Abbreviation("TIPVALues") ) {
            HandleTipValues( token );
        }
        else if( token.Abbreviation("TIPVARiance") ) {
            HandleTipVariance( token );
        }
        else if( token.Abbreviation("MOdel") ) {
            HandleModel( token );
        }
		else if( token.Abbreviation("DIscrete") ) {
			HandleDiscrete( token );
		}
		else if( token.Abbreviation("LOss") ) {
			HandleLoss( token );
		}
        else if( token.Abbreviation("TIMeslice") ) {
            HandleTimeSlice( token );
        }
        else if( token.Abbreviation("Debug") ) {
            HandleDebug( token );
        }
        else if (token.Abbreviation("PREorder") ) {
            PreOrderTraversal(token);
        }
        else if (token.Abbreviation("MRCA") ) {
            HandleMRCA( token );
        }
        else if (token.Abbreviation("EXport") ) {
            HandleExport(token);
        }
        else if ( token.Abbreviation("OPtimization") ) {
            HandleDebugOptimization( token );
        }
        else if ( token.Abbreviation("NUmopt") ) {
            NumOpt( token );
        }
        else if (token.Abbreviation("GARLAND") ) {
            HandleGarland(token);
        }
        else if( token.Abbreviation("RateTest") ) {
            HandleRateTest( token );
        }
        else if(token.Abbreviation("Taxset") ) {
            //errormsg = "Sorry, enter taxsets in the ASSUMPTIONS block";
            //throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            assumptions->HandleTaxset( token );
        }
        else if( token.Abbreviation("Printedgelengths") ) {
            HandlePrintEdgeLengths( token);
        }
        else if( token.Abbreviation("Quit") ) {
            quit_now = true;
			PrintCitations();
            message = "\nPlease remember to send bug reports to bcomeara@ucdavis.edu.\n";
            PrintMessage();

            break;
        }
        else
        {
            SkippingCommand( token.GetToken() );
            do {
                token.GetNextToken();
            } while( !token.AtEOF() && !token.Equals(";") );

            if( token.AtEOF() ) {
                errormsg = "Unexpected end of file encountered";
                throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
    }
}

/**
* @method Reset [void:protected]
 *
 * Overrides the pure virtual function in the base class.
 */
void BROWNIE::Reset()
{
    isEmpty = true;

    inf_open = false;
    quit_now = false;
    message = "";
    next_command[0] = '\0';
}

/**
* @method Report [virtual void:public]
 * @param out [ostream&] the output stream to which to write the report
 *
 * This function outputs a brief report of the contents of this BROWNIE block.
 * Overrides the pure virtual function in the base class.
 */
void BROWNIE::Report( ostream& /* out */ )
{
    message = "";
    PrintMessage();

    message = id;
    message += " block contains...";
    PrintMessage();
}

/**
* @method RunCmdLine [void:public]
 *
 * Runs the command line interpreter, allowing BROWNIE to interact with user.
 */
void BROWNIE::RunCmdLine(bool inputfilegiven, nxsstring fn)
{
    taxa = new TaxaBlock();
    trees = new TreesBlock(*taxa);
    assumptions = new AssumptionsBlock( *taxa );
    characters = new CharactersBlock( *taxa, *assumptions );
    Add( taxa );
    Add( trees );
    Add( assumptions );
    Add( characters );
    Add( this );
    cout<<endl<<endl<<endl<<"                               Brownie V2.0b7"<<endl;
	cout<<"                          Character evolution models,"<<endl;
	cout<<"                             species delimitation,"<<endl; 
	cout<<"                               and tree search"<<endl<<endl;
    cout<<"                                Brian O'Meara"<<endl;
    cout<<"                     http://www.brianomeara.info/brownie"<<endl<<endl;
    cout<<"                     Uses Paul Lewis' Nexus Class Library"<<endl;
    cout<<"                (modified to deal with continuous characters)"<<endl;
    cout<<"                        Rod Page's TreeLib & supertree,"<<endl;
    cout<<"                             Mike Sanderson's GTP,"<<endl;
    cout<<"                        and the GNU Scientific Library (v1.9)"<<endl<<endl;
    cout<<"                  Type \"help\" [no quotes] for a list of commands"<<endl;
    cout<<"                  and \"commandname ?\" for help for each command."<<endl<<endl;
    cout<<"                       Compiled on "<<__DATE__<<" at "<<__TIME__<<endl;
	cout<<"                         Using GSL version "<<GSL_VERSION<<endl<<endl;
    quit_now = false;
    if (inputfilegiven) {
        HandleExecuteCmdLine(fn);
    }
    while( !quit_now )
    {
        cerr << endl;
        cerr << "brownie> ";
        //cin.getline( next_command, 256 );
        fgets(next_command, COMMAND_MAXLEN,stdin); //use fgets since GCC on OS 10.2 has trouble with cin.
        PreprocessNextCommand();
        HandleNextCommand();
    }
}

/**
* @method Run [void:public]
 *
 * Runs the command line interpreter, allowing BROWNIE to interact with user.
 * Typically, this is the only function called in main after a BROWNIE object
 * is created.
 */
//void BROWNIE::Run()
//{
//	taxa = new TaxaBlock();
//	trees = new TreesBlock(*taxa);
//	assumptions = new AssumptionsBlock( *taxa );
//	characters = new CharactersBlock( *taxa, *assumptions );

    //	Add( taxa );
    //	Add( trees );
    //	Add( assumptions );
    //	Add( characters );
    //	Add( this );
    //      cout<<endl<<endl<<endl<<"                               Brownie V2.0b2"<<endl;
    //        cout<<"                 Testing rates of continuous character evolution"<<endl<<endl;
    //        cout<<"                                Brian O'Meara"<<endl;
    //        cout<<"                     http://www.brianomeara.info/brownie"<<endl<<endl;
    //        cout<<"                     Uses Paul Lewis' Nexus Class Library"<<endl;
    //        cout<<"                (modified to deal with continuous characters)"<<endl;
    //        cout<<"                         Rod Page's TreeLib, and code"<<endl;
    //        cout<<"                      from Ligia Mateiu & John Burkardt."<<endl<<endl;
    //        cout<<"                  Type \"help\" [no quotes] for a list of commands"<<endl;
    //        cout<<"                  and \"commandname ?\" for help for each command."<<endl<<endl;
    //	quit_now = false;
    //	while( !quit_now )
    //	{
    //		cerr << endl;
    //		cerr << "brownie> ";
    //		//cin.getline( next_command, 256 );
    //              fgets(next_command, COMMAND_MAXLEN,stdin); //use fgets since GCC on OS 10.2 has trouble with cin.
    //		PreprocessNextCommand();
    //		HandleNextCommand();
    //	}
    //}

/**
    * @method SkippingBlock [virtual void:public]
 * @param blockName [nxsstring] the unrecognized block name
 *
 * Called when program does not recognize a block name encountered in a
 * NEXUS file.  Virtual function that overrides the pure virtual function
 * in the base class Nexus.
 */
    void BROWNIE::SkippingBlock( nxsstring blockName )
    {
        message = "Skipping unknown block (";
        message += blockName;
        message += ")";
        PrintMessage();
    }

    /**
    * @method SkippingCommand [virtual void:public]
     * @param commandName [nxsstring] the name of the command being skipped
     *
     * This function is called when an unknown command named commandName is
     * about to be skipped.  This version of the function (which is identical
                                                           * to the base class version) does nothing (i.e., no warning is issued
                                                                                                      * that a command was unrecognized).  Modify this virtual function to
     * provide such warnings to the user (or eliminate it altogether since
                                          * the base class version already does what this does).
     */
    void BROWNIE::SkippingCommand( nxsstring commandName )
    {
        message = "Skipping unknown command (";
        message += commandName;
        message += ")\nType help for a list of available commands.";
        PrintMessage();
    }

    /**
    * @method SkippingDisabledBlock [virtual void:public]
     * @param blockName [nxsstring] the name of the block just exited
     *
     * Called by the Nexus object when skipping a block named blockName
     * that has been disabled. Allows program to notify user of progress
     * in parsing the NEXUS file. Virtual function that overrides the
     * pure virtual function in the base class Nexus.
     */
    void BROWNIE::SkippingDisabledBlock( nxsstring /*blockName*/ )
    {
    }

    /**
    * @method TaxonLabelToNumber [int:protected]
     * @param s [nxsstring] the taxon label to be translated to a taxon number
     *
     * The code here is modified from charactersblock.cpp
     */
    int BROWNIE::TaxonLabelToNumber( nxsstring s )
    {
        int i;
        try {
            i = 1 + taxa->FindTaxon(s);
        }
        catch( TaxaBlock::nosuchtaxon ) {
            i = 0;
        }
        return i;
    }

    /**
    * Asks user if "something" is ok, where "something" is expressed
     * in the title and message displayed.  This is a virtual function
     * so it can be overridden in a derived class to use a different
     * (perhaps graphical) means for displaying the message.
     * Note: mb_message should terminate with a quesiton mark; none
     * will be provided by this function.
     */
    bool BROWNIE::UserSaysOk( nxsstring mb_message, nxsstring mb_title )
    {
        cerr << endl;
        cerr << mb_title << endl;
        cerr << "  " << mb_message;
        cerr << " (y/n) ";
        fgets(next_command, COMMAND_MAXLEN,stdin);
        //cin.getline( next_command, COMMAND_MAXLEN ); //GCC problem
        //bool yep  = ( next_command[0] == 'y' && next_command[1] == '\0' );
        //bool nope = ( next_command[0] == 'n' && next_command[1] == '\0' );
        bool yep  = ( next_command[0] == 'y');
        bool nope = ( next_command[0] == 'n');
        int patience=5;
        int loopnumber=0;

        while( !yep && !nope )
        {
            loopnumber++;
            if (loopnumber>patience) {
                errormsg = "ERROR: You're too indecisive.";
                throw XNexus( errormsg);
            }
            cerr << endl;
            cerr << "Must answer by typing either y or n and then pressing the Enter key" << endl;
            cerr << endl;
            cerr << mb_title << endl;
            cerr << "  " << mb_message;
            cerr << " (y/n) ";
            fgets(next_command, COMMAND_MAXLEN,stdin);
            //cin.getline( next_command, COMMAND_MAXLEN ); //OS 10.2 GCC problem
            //yep  = ( next_command[0] == 'y' && next_command[1] == '\0' );
            //nope = ( next_command[0] == 'n' && next_command[1] == '\0' );
            yep  = ( next_command[0] == 'y');
            nope = ( next_command[0] == 'n');
        }

        return yep;
    }

    void BROWNIE::HandleMRCA( NexusToken& token )
    {
        nxsstring mrca_name;
        IntSet taxonnumbers;
        token.GetNextToken();

        if( token.Abbreviation("?") ) {
            message="Usage: MRCA <name> = <taxon 1> <taxon 2> >...>\n\n";
            message+="Loads the name of the most recent common ancestor of a given set of taxa (using names or numbers, as in a taxset definition).\nNote that the MRCA name is case-sensitive.";
            PrintMessage();
            token.GetNextToken();

        }
        else {
            // Token now stored should be the name of a mrca
            mrca_name = token.GetToken();

            // Now grab the equals sign
            token.GetNextToken();
            if( !token.Equals("=") ) {
                errormsg = "Expecting '=' in MRCA definition but found ";
                errormsg += token.GetToken();
                errormsg += " instead";
                throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            int totalTaxa = taxa->GetNumTaxonLabels();
            SetReader( token, totalTaxa, taxonnumbers, *this, SetReader::taxset ).Run();
            MrcaMap[mrca_name]=taxonnumbers;
        }
    }


    void BROWNIE::HandleExport( NexusToken& token)
    {
        nxsstring taxatoexport="";
        nxsstring internalstoexport="";
        double rate=1;
        double trend=0;
        nxsstring exportfname;
        int format=0; //0 means pagel
        bool replacing = false;
        bool nameprovided=true;
        bool hadfirstinternalnode=false;
        int source=0; //0 means use default character matrix, other sources involve simulation
        map<string,double> simcharmatrix;
        for(;;)
        {
            token.GetNextToken();

            if( token.Equals(";") ) {
                break;
            }
            else if( token.Abbreviation("FIle") ) {
                exportfname = GetFileName(token);
                nameprovided = true;
                //cout<<"Read file name"<<endl;
            }
            else if( token.Abbreviation("Replace") ) {
                nxsstring yesnoreplace=GetFileName(token);
                if (yesnoreplace[0] == 'n' || yesnoreplace[0] == 'N') {
                    replacing=false;
                }
                else {
                    replacing=true;
                }
            }
            else if( token.Abbreviation("Format") ) {
                nxsstring sourcestring=GetFileName(token);
                if (sourcestring[0] == 'p' || sourcestring[0] == 'P') {
                    format=0;
                }
            }

            else if( token.Abbreviation("Source") ) {
                nxsstring sourcestring=GetFileName(token);
                if (sourcestring[0] == 'c' || sourcestring[0] == 'C') {
                    source=0;
                }
                else if (sourcestring[0] == 't' || sourcestring[0] == 'T') {
                    source=1;
                }
                else {
                    //Simulate
                }
                // cout<<"Read source string"<<endl;
            }
            else if( token.Abbreviation("Trend") ) {
                nxsstring numbernexus = GetNumber(token);
                trend=atof( numbernexus.c_str() );
            }
            else if( token.Abbreviation("?") ) {
                message="Usage: Export [options...]\n\n";
                message+="Saves tree and data in a new format.\n\n";
                message+="Available options:\n\n";
                message+="Keyword ---- Option type ------------------------ Current setting --";
                message+="\nFile         <export-file-name>                   [None]";
                message+="\nReplace      yes | no                             false";
                //message+="\nFormat       Pagel | [other]                      Pagel";
                //message+="\nSource       CharMatrix | Trends                  CharMatrix";
                //message+="\nTrend        <double>                             ";
                //message+=trend;
                message+="\n\nCurrently, this only outputs the character matrix in Pagel format to a specified file.\nMake sure the tree has no polytomies.";
                PrintMessage();

            }
            else {
                errormsg = "Unexpected keyword (";
                errormsg += token.GetToken();
                errormsg += ") encountered reading PrintEdgeLengths command";
                throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        if (!nameprovided) {
            errormsg = "Must provide an output file name\n";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
        else {
            exportf.open( exportfname.c_str(), ios::out | ios::app );
            Tree t = intrees.GetIthTree(chosentree-1);
            //First do uppass to get Pagel like node numbers
            PreorderIterator <Node> m (t.GetRoot());
            Node *p = m.begin();
            int count=1+(characters->GetNTax());
            while (p)
            {
                if(!(p->IsLeaf())) {
                    ostringstream NewLabelStream;
                    NewLabelStream<<count;
                    p->SetLabel(NewLabelStream.str());
                    count--;
                }
                p=m.next();
            }
            //Now do downpass to print out results
            // int count = 0;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            Node *anc = n.begin();
            // Node *a=n.begin();
            // Node *p=n.begin();
            //Node *b=n.begin();
            //Node *mrca=n.begin();
            //exportf<<"#PAG"<<endl;
            //cout<<"Successfully saved\n";
            if (source==0) {
                exportf<<characters->GetNTax()<<" "<< characters->GetNChar()<<"\n";
            }
            else if (source==1) {
                int ntax=characters->GetNTax();
                exportf<<ntax;
                exportf<<" 1\n";
                simcharmatrix=SimulateBrownian(trend,1,0);
            }
            while (q->GetAnc())
            {
                string qlabel;
                qlabel=q->GetLabel();
                if (debugmode) {
                    cout<<qlabel<<endl;
                }
                // exportf << setiosflags(ios::right) << setw(4) << count					// arbitrary counter
                //      << " " << setw(8) << q->GetEdgeLength() 							// edge length
                //      << " " << setw(32) << setiosflags(ios::left) << q->GetLabel()	// node label
                //      ;
                //cout <<setw(32) << setiosflags(ios::left) << q->GetLabel()<<", ";
                if (q->IsLeaf()) {
                    taxatoexport+=(q->GetLabel()).c_str();
                    taxatoexport+=", ";
                    //exportf  << setiosflags(ios::left) << q->GetLabel()<<", ";
                }
                else if (q->GetAnc()) {
                    internalstoexport+=(q->GetLabel()).c_str();
                    internalstoexport+=", ";
                }
                if (q->GetAnc()) {
                    anc=q->GetAnc();
                    // if (strlen((anc->GetLabel()).c_str())==0) {
                    //     count++;
                    //     ostringstream NewLabelStream;
                    //     NewLabelStream<<count;
                    //     anc->SetLabel(NewLabelStream.str());
                    // }
                    if (q->IsLeaf()) {
                        taxatoexport+=(anc->GetLabel()).c_str();
                        taxatoexport+=", ";
                    }
                    else {
                        internalstoexport+=(anc->GetLabel()).c_str();
                        internalstoexport+=", ";
                    }
                    //exportf<<anc->GetLabel()<<", ";
                }
                else {
                    if (q->IsLeaf()) {
                        taxatoexport+=", ";
                    }
                    else if (q->GetAnc()) {
                        internalstoexport+=", ";
                    }
                    //exportf<<count<<", ";
                }
                double outputvalue=q->GetEdgeLength();
                if (q->IsLeaf()) {
                    taxatoexport+=outputvalue;
                }
                else if (q->GetAnc()) {
                    internalstoexport+=outputvalue;
                }
                //exportf << q->GetEdgeLength();
                if (q->IsLeaf())
                {
                    if (source==0) {
                        nxsstring currenttaxonlabel=(q->GetLabel()).c_str();
                        int taxonnumber=TaxonLabelToNumber( currenttaxonlabel );
                        for (int charnumber=0; charnumber<(characters->GetNChar()); charnumber++) {
                            taxatoexport+=", ";
                            taxatoexport+=characters->GetValue( taxonnumber, charnumber, true);
                        }
                    }
                    else if (source==1) {
                        map<string,double>::iterator pos;
                        pos=simcharmatrix.find(qlabel);
                        taxatoexport+=", ";
                        taxatoexport+=pos->second;
                        //exportf<<", "<<pos->second;
                    }

                }
                if (q->IsLeaf()) {
                    taxatoexport+="\n";
                }
                else if (q->GetAnc()) {
                    internalstoexport+="\n";
                }

                // exportf << endl;
                q = n.next();
            }
            exportf<<taxatoexport;
            exportf<<internalstoexport;
            exportf.close();
        }
    }

    nxsstring BROWNIE::PipeLeafGTP ()
    {
        int ntax=taxa->GetNumTaxonLabels();
        nxsstring output="";
        nxsstring TaxonLabel="";
        TaxonLabel+=cur->GetLabel();
        int TaxonNumber=taxa->FindTaxon(TaxonLabel);
        nxsstring NewLabel="taxon";
        if(TaxonNumber>=convertsamplestospecies.size()) {
            cout<<"Current taxon label is "<<TaxonLabel<<" which has number "<<TaxonNumber<<" but the convertsamplestospecies vector is of length "<<convertsamplestospecies.size()<<":"<<endl;
            for (int i=0;i<convertsamplestospecies.size();i++) {
                cout<<convertsamplestospecies[i]<<" ";
            }
            cout<<endl;
        }
        assert(TaxonNumber<convertsamplestospecies.size());
        NewLabel+=convertsamplestospecies[TaxonNumber];
        output=NEXUSString (NewLabel);
        return output;
    }

    nxsstring BROWNIE::PipeLeafSpeciesTree ()
    {
        nxsstring output="";
        output+=(cur->GetLabel());
        return output;
    }

    nxsstring BROWNIE::PipeGTP (Tree t)
    {
        nxsstring TreeDescription="";
        cur = t.GetRoot();

        while (cur)
        {
            if (cur->GetChild())
            {
                TreeDescription+=PipeLeftParenthesis ();
                stk.push (cur);
                cur = cur->GetChild();
            }
            else
            {
                TreeDescription+=PipeLeafGTP ();
                while (!stk.empty() && (cur->GetSibling() == NULL))
                {
                    TreeDescription+=PipeRightParenthesis ();
                    cur = stk.top();
                    PipeInternal ();
                    stk.pop();
                }
                if (stk.empty())
                    cur = NULL;
                else
                {
                    TreeDescription+=PipeSiblingSymbol ();
                    cur = cur->GetSibling();
                }
            }
        }
        TreeDescription+=PipeEndOfTreeGTP ();
        return TreeDescription;
    }


    nxsstring BROWNIE::PipeSpeciesTree (ContainingTree *SpeciesTreePtr)
    {
        nxsstring TreeDescription="";
        cur = SpeciesTreePtr->GetRoot();

        while (cur)
        {
            if (cur->GetChild())
            {
                TreeDescription+=PipeLeftParenthesis ();
                stk.push (cur);
                cur = cur->GetChild();
            }
            else
            {
                TreeDescription+=PipeLeafSpeciesTree ();
                while (!stk.empty() && (cur->GetSibling() == NULL))
                {
                    TreeDescription+=PipeRightParenthesis ();
                    cur = stk.top();
                    TreeDescription+=PipeInternal ();
                    stk.pop();
                }
                if (stk.empty())
                    cur = NULL;
                else
                {
                    TreeDescription+=PipeSiblingSymbol ();
                    cur = cur->GetSibling();
                }
            }
        }
        TreeDescription+=PipeEndOfTreeSpeciesTree ();
        return TreeDescription;
    }



    nxsstring BROWNIE::PipeEndOfTreeGTP ()
    {
        nxsstring output="";
        if (outputallatonce) {
            output+=":";
            output+=trees->GetTreeWeight(chosentree);
        }
        else {
            output+=":1.00"; //Since we're piping one tree at a time, we don't need weights
        }

        output+=';';
        return output;
    }

    nxsstring BROWNIE::PipeEndOfTreeSpeciesTree ()
    {
        nxsstring output="";
        output+=';';
        return output;
    }


    //------------------------------------------------------------------------------
    nxsstring BROWNIE::PipeLeftParenthesis ()
    {
        nxsstring output="";
        output+='(';
        return output;
    }

    //------------------------------------------------------------------------------
    nxsstring BROWNIE::PipeRightParenthesis ()
    {
        nxsstring output="";
        output+=')';
        return output;
    }

    //------------------------------------------------------------------------------
    nxsstring BROWNIE::PipeSiblingSymbol ()
    {
        nxsstring output="";
        output+=',';
        return output;
    }


    //------------------------------------------------------------------------------
    nxsstring BROWNIE::PipeInternal ()
    {
        nxsstring output="";
        if (cur->GetLabel() != "")
            output+=cur->GetLabel();
        return output;
    }


    nxsstring BROWNIE::ReturnFinalSpeciesTree (Tree t)
    {
        nxsstring TreeDescription="";
        cur = t.GetRoot();

        while (cur)
        {
            if (cur->GetChild())
            {
                TreeDescription+=PipeLeftParenthesis ();
                stk.push (cur);
                cur = cur->GetChild();
            }
            else
            {
                TreeDescription+=PipeLeafSpeciesTree ();
                while (!stk.empty() && (cur->GetSibling() == NULL))
                {
                    TreeDescription+=PipeRightParenthesis ();
                    cur = stk.top();
                    PipeInternal ();
                    stk.pop();
                }
                if (stk.empty())
                    cur = NULL;
                else
                {
                    TreeDescription+=PipeSiblingSymbol ();
                    cur = cur->GetSibling();
                }
            }
        }
        TreeDescription+=PipeEndOfTreeSpeciesTree();
        return TreeDescription;
    }

    nxsstring BROWNIE::PipeLeafFinalSpeciesTree ()
    {
        int speciesnumber;
        nxsstring specieslabel=cur->GetLabel();
        string speciesstring=specieslabel.c_str();
        size_t index = speciesstring.find('t');
        speciesstring.erase(index,5); //erase "taxon"
                                      //cout<<"species string now "<<speciesstring<<endl;
        speciesnumber=atoi(speciesstring.c_str());
        nxsstring finallabel="(";
        int ntax=taxa->GetNumTaxonLabels();
        int timesprinted=0;
        for (int taxon=0; taxon<ntax; taxon++) {
            if (convertsamplestospecies[taxon]==speciesnumber) {
                if (timesprinted>0) {
                    finallabel+=",";
                }
                finallabel+=GetTaxonLabel(taxon);
                timesprinted++;
            }
        }
        finallabel+=")";
        // finallabel+=specieslabel;
        return finallabel;
    }


    map<string,double> BROWNIE::SimulateBrownian(double trend,double rate,double rootstate)
    {
        //variables
        map<string, double> newcharmatrix;
        message="Now simulating, trend = ";
        message+=trend;
        message+=" rate = ";
        message+=rate;
        message+=" state at root = ";
        message+=rootstate;
        PrintMessage();
        Tree t = intrees.GetIthTree(chosentree-1);
        PreorderIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        Node *ancestor=n.begin();
        while (q)
        {
            //cout <<setw(32) << setiosflags(ios::left) << q->GetLabel()<<", ";
            string qlabel;
            double newstate;
            qlabel=q->GetLabel();
            if (!q->GetAnc()) {
                ostringstream NewLabelStream;
                NewLabelStream<<rootstate;
                q->SetLabel(NewLabelStream.str());
            }
            else {
                //cout<<"Not doing root node\n";
                ancestor=q->GetAnc();
                double ancestralstate=atof((ancestor->GetLabel()).c_str());
                double brlen=q->GetEdgeLength();
                newstate=ancestralstate+(trend*brlen)+(gsl_ran_gaussian ( r, sqrt(rate*brlen)));
                ostringstream NewLabelStream;
                NewLabelStream<<newstate;
                q->SetLabel(NewLabelStream.str());
                if (debugmode) {
                    cout<<"Ancestral state = "<<ancestralstate<<", brlen = "<<brlen<<", newstate = "<<newstate<<endl;
                }
            }
            if (q->IsLeaf())
            {
                newcharmatrix[qlabel]=newstate;
            }
            //cout<<"Moving on to next node\n";
            q = n.next();
        }
        message="Done simulating";
        PrintMessage();
        return newcharmatrix;
    }

    void BROWNIE::PreOrderTraversal(NexusToken& token)
    {
        //for(;;)
        // {
        //     token.GetNextToken();
        //      if( token.Equals(";") ) {
        //          break;
        //       }
        //   }
        Tree t = intrees.GetIthTree(chosentree-1);
        //t.SetInternalLabels(true);
        PreorderIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        Node *ancestor=n.begin();
        bool notfirstleaf=false;
        int count=0;
        while (q)
        {
            count++;
            cout << setiosflags(ios::right) << setw(4) << count					// arbitrary counter
                << " " << setw(8) << q->GetEdgeLength() 							// edge length
                << " " << setw(32) << setiosflags(ios::left) << q->GetLabel()	// node label
                ;
            if (q->IsLeaf())
            {
                cout << " [LEAF]";
                cout<<" length from root="<<q->GetHeight();
                notfirstleaf=true;
                ancestor=q;
                ancestor=ancestor->GetAnc();
                cout<<" label of ancestor = "<<ancestor->GetLabel();
            }
            else
            {
                cout << " [INTERNAL]";
                cout<<" length from root="<<q->GetHeight();
                char * NewLabel;
                sprintf(NewLabel,"%u",count);
                cout<<" did sprintf ";
                q->SetLabel(NewLabel);
                cout<<" new label = "<<q->GetLabel();
                cout<<" twice new label is "<<2*atoi((q->GetLabel()).c_str());
            }
            cout << endl;
            q = n.next();
        }
    }


    void BROWNIE::HandlePrintEdgeLengths(NexusToken& token)
    {
        for(;;)
        {
            token.GetNextToken();

            if( token.Equals(";") ) {
                break;
            }
            else {
                errormsg = "Unexpected keyword (";
                errormsg += token.GetToken();
                errormsg += ") encountered reading PrintEdgeLengths command";
                throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        Tree t = intrees.GetIthTree(chosentree-1);
        cout << endl << "Postorder traversal of tree" << endl;
        cout << "Node Length     Label                           Type" << endl;
        cout << "----------------------------------------------------" << endl;
        int count = 0;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        Node *a=n.begin();
        Node *p=n.begin();
        Node *b=n.begin();
        Node *mrca=n.begin();

        bool notfirstleaf=false;
        while (q)
        {
            count++;
            cout << setiosflags(ios::right) << setw(4) << count					// arbitrary counter
                << " " << setw(8) << q->GetEdgeLength() 							// edge length
                << " " << setw(32) << setiosflags(ios::left) << q->GetLabel()	// node label
                ;
            if (q->IsLeaf())
            {
                cout << " [LEAF]";
                if (q->IsMarked()) {
                    cout <<" Node is marked  ModelCategorySize: ";
                }
                else {
                    cout <<" Node not marked  ModelCategorySize: ";
                }
                vector<double> modelcatoutput;
                modelcatoutput=q->GetModelCategory();
                //gsl_vector modelcatoutput(q->GetModelCategory());
                //cout<<modelcatoutput->size;
                cout<<" length from root="<<q->GetPathLength();
                float pathlength=0;
                a=q;
                if(notfirstleaf) {
                    //  cout<<"Get MRCA of "<< q->GetLabel()<<" and "<<p->GetLabel()<<endl;
                    bool mrcanotfound=true;
                    while (a->GetAnc() && mrcanotfound) {
                        a=a->GetAnc();
                        b=p;
                        while (b->GetAnc() && mrcanotfound) {
                            b=b->GetAnc();
                            if (a==b) {
                                mrca=a;
                                mrcanotfound=false;
                                while(mrca->GetAnc()) {
                                    pathlength+=mrca->GetEdgeLength();
                                    //         cout<<"Edge length: "<<mrca->GetEdgeLength()<<" Total: "<<pathlength<<endl;
                                    mrca=mrca->GetAnc();
                                }
                            }
                        }
                    }
                }
                notfirstleaf=true;
                p=q;
                //  cout<<"Root to MRCA length="<<pathlength<<endl;
                //    while (aanc != Root && mrcanotfound) {
                //        banc=b;
                //        aanc=aanc->GetAnc();
                //       while (banc != Root && mrcanotfound) {
                //           banc=banc->GetAnc();
                //          if (aanc == banc) {
                //              mrcaptr=aanc;
                //               mrcanotfound=false;
                //           }
                //        }
                //    }
                //    return mrcaptr;


                //a=q;
                //cout<<endl<<"interim path length 1: "<<pathlength<<endl;
                //while (a->GetAnc()) {
                //   r=a->GetAnc();
                //  pathlength+=r->GetEdgeLength();
                // cout<<"interim: "<<pathlength<<" and edge "<<r->GetEdgeLength()<<endl;
                // a=r;
                //}
                //cout<<endl<<"Total length="<<pathlength<<endl;
            }
            else
            {
                cout << " [INTERNAL]";
                if (q->IsMarked()) {
                    cout <<" Node is marked  ModelCategorySize: ";
                }
                else {
                    cout <<" Node not marked  ModelCategorySize: ";
                }
                //gsl_vector *modelcatoutput;
                //modelcatoutput=gsl_vector_calloc(q->GetModelCategory());
                //gsl_vector modelcatoutput(q->GetModelCategory());
                //cout<<modelcatoutput->size;
                cout<<" length from root="<<q->GetPathLength();
            }
            cout << endl;
            q = n.next();
        }
        cout << "----------------------------------------------------" << endl;
    }

/**
* Copied from charactersblock.cpp
 * @method GetTaxonLabel [nxsstring:public]
 * @param i [int] the taxon's position in the taxa block
 *
 * Returns label for taxon number i (i ranges from 0 to ntax-1).
 */
nxsstring BROWNIE::GetTaxonLabel( int i )
{
    nxsstring s = taxa->GetTaxonLabel(i);
    return s;
}

//From TreeLib
nxsstring& BROWNIE::blanks_to_underscores( nxsstring& s )
{
    int len = s.length();
    for( int k = 0; k < len; k++ ) {
        if( s[k] == ' ' )
            s[k] = '_';
    }
    return s;
}

nxsstring& BROWNIE::underscores_to_blanks( nxsstring& s )
{
    int len = s.length();
    for( int k = 0; k < len; k++ ) {
        if( s[k] == '_' )
            s[k] = ' ';
    }
    return s;
}



void BROWNIE::HandleVCV(NexusToken& token)
{
    bool donenothing=true;
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            if (donenothing) {
                message = "Usage: VCV taxset=<taxset-name>\n\n";
                PrintMessage();
            }
            break;
        }
        else if (token.Equals("?") ) {
            donenothing=false;
            message = "Usage: VCV taxset=<taxset-name>\n\n";
            PrintMessage();
        }
        else if (token.Abbreviation("Taxset") ) {
            donenothing=false;
            if (trees->GetNumTrees()<1) {
                errormsg = "Error: No valid trees are loaded yet\nYou can use Gettrees to load the trees.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            nxsstring currenttaxonlabel;
            nxsstring currenttaxonlabeltaxset;
std::string currenttaxonlabelnodeptr;
            nxsstring chosentaxset;
            int ntaxintaxset=0;
            chosentaxset=GetFileName(token);
            message="Taxset ";
            message+=chosentaxset.c_str();
            message+=": \n";
            IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
            if (taxonlist.empty()) {
                errormsg= "Error: Taxset ";
                errormsg+=chosentaxset.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntaxintaxset++;
                currenttaxonlabel=taxa->GetTaxonLabel(*xi);
                message+="\t";
                message+=(*xi + 1);
                message+="=";
                message+=currenttaxonlabel.c_str();
                message+=" \n";
            }
            Tree t = intrees.GetIthTree(chosentree-1);
            gsl_matrix* VCVmatrix=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
            VCVmatrix=GetVCV(chosentaxset);
            message+="With stem\n";
            for (int currentrow=0;currentrow<ntaxintaxset;currentrow++) {
                for (int currentcol=0;currentcol<ntaxintaxset;currentcol++) {
                    message+=gsl_matrix_get(VCVmatrix,currentrow,currentcol);
                    message+="\t";
                }
                message+="\n";
            }
            message+="\n";
            VCVmatrix=DeleteStem(GetVCV(chosentaxset));
            message+="Without stem\n";
            for (int currentrow=0;currentrow<ntaxintaxset;currentrow++) {
                for (int currentcol=0;currentcol<ntaxintaxset;currentcol++) {
                    message+=gsl_matrix_get(VCVmatrix,currentrow,currentcol);
                    message+="\t";
                }
                message+="\n";
            }
            PrintMessage();
			gsl_matrix_free(VCVmatrix);
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading VCV command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

void BROWNIE::HandleTipValues(NexusToken& token)
{
    nxsstring chosentaxset;
    bool taxsetloaded=false;
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") && taxsetloaded!=true) {
            message="Error: Insufficient arguments. At least include taxset=taxset-name.\n\n";
            message += "Usage: TipValues taxset=taxset-name\n";
            PrintMessage();

            break;
        }
        else if (token.Equals(";") && taxsetloaded) {
            //Do tipvalue thing once everything else has been loaded
            int ntaxintaxset=0;
            IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
            if (taxonlist.empty()) {
                errormsg= "Error: Taxset ";
                errormsg+=chosentaxset.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntaxintaxset++;
            }
            gsl_vector *tipvector;
            tipvector=gsl_vector_calloc(ntaxintaxset);
            //gsl_vector tipvector(ntaxintaxset,0);
            int rowcount=-1;
IntSet::const_iterator ri;
            for( ri = taxonlist.begin(); ri != taxonlist.end(); ri++ ) {
                rowcount++;
                gsl_vector_set(tipvector,rowcount,characters->GetValue( *ri, chosenchar-1 ,true));
                message = "  ";
                message+=gsl_vector_get(tipvector,rowcount);
                message+="\n";
                PrintMessage();
            }
			gsl_vector_free(tipvector);
            break;
        }
        else if (token.Equals("?") ) {
            message = "Usage: TipValues taxset=taxset-name\nNote that char defaults to 1.";
            PrintMessage();
        }
        else if (token.Abbreviation("Taxset") ) {
            nxsstring currenttaxonlabel;
            nxsstring currenttaxonlabeltaxset;
            chosentaxset=GetFileName(token);
            message="Taxset " ;
            message+=chosentaxset.c_str();
            message+=", character ";
            message+=chosenchar;
            message+=":\n";
            PrintMessage();
            taxsetloaded=true;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading TipValues command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

//Returns a vector containing tip values IN THE SAME ORDER AS THE TAXA IN THE TAXSET
gsl_vector* BROWNIE::GetTipValues(nxsstring chosentaxset, int chosenchar)
{
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }

IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;

    }
    gsl_vector *tipvalues;
    tipvalues=gsl_vector_calloc(ntaxintaxset);
    //gsl_vector tipvalues(ntaxintaxset,0);
    int rowcount=-1;
IntSet::const_iterator ri;
    for( ri = taxonlist.begin(); ri != taxonlist.end(); ri++ ) {
        rowcount++;
        gsl_vector_set(tipvalues,rowcount,characters->GetValue( *ri, chosenchar-1, true));
    }

    return tipvalues;
}

//Returns a VCV matrix with columns and rows IN THE SAME ORDER AS THE TAXA IN THE TAXSET
gsl_matrix* BROWNIE::GetVCV(nxsstring chosentaxset)
{
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
	std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
	IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    Tree t = intrees.GetIthTree(chosentree-1);
    gsl_matrix *VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    int rowcount=-1;
    nxsstring rtaxon;
    nxsstring ctaxon;
	IntSet::const_iterator ri;
	
    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=-1;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
		IntSet::const_iterator ci;
        for( ci = taxonlistc.begin(); ci != taxonlistc.end(); ci++ ) {
            colcount++;
            ctaxon=taxa->GetTaxonLabel(*ci);
            nxsstring ctaxonUnderscores=blanks_to_underscores(ctaxon);
            nxsstring ctaxonBlanks=underscores_to_blanks(ctaxon);
            Node *cleafptr;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            bool goodmatch=false;
            while (q)
            {
                if (q->IsLeaf())
                {
                    string qlabel;
                    qlabel=q->GetLabel();
					// nxsstring qlabel2<<qlabel;
     //qlabel2=blanks_to_underscores(qlabel2);
                    if (qlabel==ctaxon.c_str() || qlabel==ctaxonUnderscores.c_str() || qlabel==ctaxonBlanks.c_str()) {
                        cleafptr=q;
                        goodmatch=true;
                    }
                }
                q = n.next();
            }
            if (goodmatch==false) {
                errormsg= "Error: there was trouble identifying taxon ";
                errormsg+=ctaxon.c_str();
                errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
                throw XNexus (errormsg );
            }
            float pathlength=0;
            Node *a;
            a=rleafptr;
            if (debugmode) {
                cout<<"Taxa "<<rleafptr->GetLabel()<<" (rtaxon="<<rtaxon<<") and "<<cleafptr->GetLabel()<<" (ctaxon="<<ctaxon<<"): ";
            }
            if (cleafptr==rleafptr) { //We need the pendant edge length
                pathlength+=cleafptr->GetEdgeLength();
            }
            bool mrcanotfound=true;
            if (debugmode) {
                cout<<pathlength;
            }
            while (a->GetAnc() && mrcanotfound) { //General strategy here:
                a=a->GetAnc();
                Node *b;
                b=cleafptr;
                while (b->GetAnc() && mrcanotfound) {
                    b=b->GetAnc();
                    if (a==b) {
                        Node* mrca;
                        mrca=a;
                        mrcanotfound=false;
                        while(mrca->GetAnc()) {
                            pathlength+=mrca->GetEdgeLength();
                            if (debugmode) {
                                cout<<" + "<<mrca->GetEdgeLength();
                            }
                            mrca=mrca->GetAnc();
                        }
                    }
                }
            }
            gsl_matrix_set(VCV,rowcount,colcount,pathlength);
            if (debugmode) {
                cout<<" = "<<pathlength<<endl;
            }
        }
    }
    //matrixsingular=TestSingularity(VCV);
    return VCV;
}

//Returns a VCV matrix with columns and rows IN THE SAME ORDER AS THE TAXA IN THE TAXSET, with edge lengths all raised to kappa power
gsl_matrix* BROWNIE::GetVCVwithKappa(nxsstring chosentaxset,double kappa)
{
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
	std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
	IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    Tree t = intrees.GetIthTree(chosentree-1);
    gsl_matrix *VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    int rowcount=-1;
    nxsstring rtaxon;
    nxsstring ctaxon;
	IntSet::const_iterator ri;
	
    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=-1;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
		IntSet::const_iterator ci;
        for( ci = taxonlistc.begin(); ci != taxonlistc.end(); ci++ ) {
            colcount++;
            ctaxon=taxa->GetTaxonLabel(*ci);
            nxsstring ctaxonUnderscores=blanks_to_underscores(ctaxon);
            nxsstring ctaxonBlanks=underscores_to_blanks(ctaxon);
            Node *cleafptr;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            bool goodmatch=false;
            while (q)
            {
                if (q->IsLeaf())
                {
                    string qlabel;
                    qlabel=q->GetLabel();
					// nxsstring qlabel2<<qlabel;
     //qlabel2=blanks_to_underscores(qlabel2);
                    if (qlabel==ctaxon.c_str() || qlabel==ctaxonUnderscores.c_str() || qlabel==ctaxonBlanks.c_str()) {
                        cleafptr=q;
                        goodmatch=true;
                    }
                }
                q = n.next();
            }
            if (goodmatch==false) {
                errormsg= "Error: there was trouble identifying taxon ";
                errormsg+=ctaxon.c_str();
                errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
                throw XNexus (errormsg );
            }
            float pathlength=0;
            Node *a;
            a=rleafptr;
            if (debugmode) {
                cout<<"Taxa "<<rleafptr->GetLabel()<<" (rtaxon="<<rtaxon<<") and "<<cleafptr->GetLabel()<<" (ctaxon="<<ctaxon<<"): ";
            }
            if (cleafptr==rleafptr) { //We need the pendant edge length
                pathlength+=pow(cleafptr->GetEdgeLength(),kappa);
            }
            bool mrcanotfound=true;
            if (debugmode) {
                cout<<pathlength;
            }
            while (a->GetAnc() && mrcanotfound) { //General strategy here:
                a=a->GetAnc();
                Node *b;
                b=cleafptr;
                while (b->GetAnc() && mrcanotfound) {
                    b=b->GetAnc();
                    if (a==b) {
                        Node* mrca;
                        mrca=a;
                        mrcanotfound=false;
                        while(mrca->GetAnc()) {
                            pathlength+=pow(mrca->GetEdgeLength(),kappa);
                            if (debugmode) {
                                cout<<" + "<<pow(mrca->GetEdgeLength(),kappa);
                            }
                            mrca=mrca->GetAnc();
                        }
                    }
                }
            }
            gsl_matrix_set(VCV,rowcount,colcount,pathlength);
            if (debugmode) {
                cout<<" = "<<pathlength<<endl;
            }
        }
    }
    //matrixsingular=TestSingularity(VCV);
    return VCV;
}


//Returns a VCV matrix with columns and rows IN THE SAME ORDER AS THE TAXA IN THE TAXSET
gsl_matrix* BROWNIE::GetVCVwithTree(nxsstring chosentaxset,Tree t)
{
    //SimpleLCAQuery lca;
    //lca.SetTree (&t);
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
	std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
	IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    gsl_matrix *VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    int rowcount=-1;
    nxsstring rtaxon;
    nxsstring ctaxon;
	IntSet::const_iterator ri;
	
    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=-1;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
		IntSet::const_iterator ci;
        for( ci = taxonlistc.begin(); ci != taxonlistc.end(); ci++ ) {
            colcount++;
            ctaxon=taxa->GetTaxonLabel(*ci);
            nxsstring ctaxonUnderscores=blanks_to_underscores(ctaxon);
            nxsstring ctaxonBlanks=underscores_to_blanks(ctaxon);
            Node *cleafptr;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            bool goodmatch=false;
            while (q)
            {
                if (q->IsLeaf())
                {
                    string qlabel;
                    qlabel=q->GetLabel();
                    // nxsstring qlabel2<<qlabel;
                    //qlabel2=blanks_to_underscores(qlabel2);
                    if (qlabel==ctaxon.c_str() || qlabel==ctaxonUnderscores.c_str() || qlabel==ctaxonBlanks.c_str()) {
                        cleafptr=q;
                        goodmatch=true;
                    }
                }
                q = n.next();
            }
            if (goodmatch==false) {
                errormsg= "Error: there was trouble identifying taxon ";
                errormsg+=ctaxon.c_str();
                errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
                throw XNexus (errormsg );
            }
            float pathlength=0;
            Node *a;
            a=rleafptr;
            if (debugmode) {
                cout<<"Taxa "<<rleafptr->GetLabel()<<" (rtaxon="<<rtaxon<<") and "<<cleafptr->GetLabel()<<" (ctaxon="<<ctaxon<<"): ";
            }
            if (cleafptr==rleafptr) { //We need the pendant edge length
                pathlength+=cleafptr->GetEdgeLength();
            }
            bool mrcanotfound=true;
            if (debugmode) {
                cout<<pathlength;
            }
            while (a->GetAnc() && mrcanotfound) { //General strategy here:
                a=a->GetAnc();
                Node *b;
                b=cleafptr;
                while (b->GetAnc() && mrcanotfound) {
                    b=b->GetAnc();
                    if (a==b) {
                        Node* mrca;
                        mrca=a;
                        mrcanotfound=false;
                        while(mrca->GetAnc()) {
                            pathlength+=mrca->GetEdgeLength();
                            if (debugmode) {
                                cout<<" + "<<mrca->GetEdgeLength();
                            }
                            mrca=mrca->GetAnc();
                        }
                    }
                }
            }
            gsl_matrix_set(VCV,rowcount,colcount,pathlength);
            if (debugmode) {
                cout<<" = "<<pathlength<<endl;
            }
        }
    }
    //matrixsingular=TestSingularity(VCV);
    return VCV;
}


gsl_matrix* BROWNIE::GetStartStopTimesforOneState(nxsstring chosentaxset, int selectedstate) {
	
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
	std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
	IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    Tree t = intrees.GetIthTree(chosentree-1);
    gsl_matrix * StartStopTimes=gsl_matrix_calloc(ntaxintaxset,6*ntaxintaxset);
    gsl_matrix * VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    VCV=DeleteStem(GetVCV(chosentaxset));
    double TotalTime=gsl_matrix_max (VCV);
    int rowcount=-1;
    nxsstring rtaxon;
	IntSet::const_iterator ri;
	
    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=0;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
        double pathlength=0;
        Node *a;
        a=rleafptr;
        while (a->GetAnc()) {
            vector<int> stateordervector(a->GetStateOrder());
            vector<double>  statetimesvector(a->GetStateTimes());
            double edgelength=a->GetEdgeLength();
            pathlength+=edgelength;
            double elapsedlength=0;
            for (int element=0;element<stateordervector.size();element++) {
                double segmentlength=statetimesvector[element];
              //  for (int position=0; position<staterestrictionvector.size(); position++) {
                  //  if (staterestrictionvector[position]==selectedstate) { //position is the called state, staterestrictionvector.at(position) gives the rate category to which that will be assigned
				if (staterestrictionvector[(stateordervector[element])]==selectedstate) {
					gsl_matrix_set(StartStopTimes,rowcount,colcount,TotalTime-edgelength+elapsedlength);
					colcount++;
					gsl_matrix_set(StartStopTimes,rowcount,colcount,TotalTime-edgelength+elapsedlength+segmentlength);
					colcount++;
				}
				elapsedlength+=segmentlength;
                   // }
               // }
            }
            a=a->GetAnc();
        }
    }
	gsl_matrix_free(VCV);
    return StartStopTimes;
}


/*Returns a VCV matrix for the selected taxset for one model
* with columns and rows IN THE SAME ORDER AS THE TAXA
* For each edge, there is a ModelCategory vector.
* For each ModelCategory (such as a morphological state), the corresponding
* entry in the vector is the amount of time spent in that Category on that edge.
*/
gsl_matrix* BROWNIE::GetVCVforOneModel(nxsstring chosentaxset, int selectedmodel)
{
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    Tree t = intrees.GetIthTree(chosentree-1);
    gsl_matrix *VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    int rowcount=-1;
    nxsstring rtaxon;
    nxsstring ctaxon;
IntSet::const_iterator ri;

    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=-1;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
IntSet::const_iterator ci;
        for( ci = taxonlistc.begin(); ci != taxonlistc.end(); ci++ ) {
            colcount++;
            ctaxon=taxa->GetTaxonLabel(*ci);
            nxsstring ctaxonUnderscores=blanks_to_underscores(ctaxon);
            nxsstring ctaxonBlanks=underscores_to_blanks(ctaxon);
            Node *cleafptr;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            bool goodmatch=false;
            while (q)
            {
                if (q->IsLeaf())
                {
                    string qlabel;
                    qlabel=q->GetLabel();
                    // nxsstring qlabel2<<qlabel;
                    //qlabel2=blanks_to_underscores(qlabel2);
                    if (qlabel==ctaxon.c_str() || qlabel==ctaxonUnderscores.c_str() || qlabel==ctaxonBlanks.c_str()) {
                        cleafptr=q;
                        goodmatch=true;
                    }
                }
                q = n.next();
            }
            if (goodmatch==false) {
                errormsg= "Error: there was trouble identifying taxon ";
                errormsg+=ctaxon.c_str();
                errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
                throw XNexus (errormsg );
            }
            double pathlength=0;
            Node *a;
            a=rleafptr;
            if (debugmode) {
                cout<<"Taxa "<<rleafptr->GetLabel()<<" (rtaxon="<<rtaxon<<") and "<<cleafptr->GetLabel()<<" (ctaxon="<<ctaxon<<"): ";
            }
            if (cleafptr==rleafptr) { //We need the pendant edge length
                                      // float unweightedpathlength=cleafptr->GetEdgeLength();
                                      //gsl_vector *modelcategoryvector;
                                      //modelcategoryvector=gsl_vector_calloc(cleafptr->GetModelCategory());
                vector<double> modelcategoryvector(cleafptr->GetModelCategory());
                // if (modelcategoryvector->size==0) {
                //     pathlength+=0;
                //}
                //else if(modelcategoryvector->size<(selectedmodel-1)) {
                //			errormsg= "Error: ModelCategory vector too short for the number of submodels ";
                //			throw XNexus (errormsg );
                //		}
                //		else {
                for (int position=0;position<staterestrictionvector.size();position++) {
                    if (staterestrictionvector[position]==selectedmodel) { //position is the called state, staterestrictionvector.at(position) gives the rate category to which that state will be assigned
                        pathlength+=modelcategoryvector[selectedmodel];
                    }
                }
                }
            bool mrcanotfound=true;
            if (debugmode) {
                cout<<pathlength;
            }
            while (a->GetAnc() && mrcanotfound) { //General strategy here:
                a=a->GetAnc();
                Node *b;
                b=cleafptr;
                while (b->GetAnc() && mrcanotfound) {
                    b=b->GetAnc();
                    if (a==b) {
                        Node* mrca;
                        mrca=a;
                        mrcanotfound=false;
                        while(mrca->GetAnc()) {
                            vector<double> modelcategoryvector(mrca->GetModelCategory());
                            pathlength+=modelcategoryvector[selectedmodel];
                            if (debugmode) {
                                cout<<" + "<<modelcategoryvector[selectedmodel];
                            }
                            mrca=mrca->GetAnc();
                        }
                    }
                }
            }
            gsl_matrix_set(VCV,rowcount,colcount,pathlength);
            if (debugmode) {
                cout<<" = "<<pathlength<<endl;
            }
            }
        }
//matrixsingular=TestSingularity(VCV);
return VCV;
    }


/*Returns a VCV matrix for the selected taxset for one model
* with columns and rows IN THE SAME ORDER AS THE TAXA
* For each edge, there is a ModelCategory vector.
* For each ModelCategory (such as a morphological state), the corresponding
* entry in the vector is the amount of time spent in that Category on that edge.
* This returns a VCV consisting only of branches with changes on them (wantchangeedges=true)
* Or a VCV consisting only of branches with no changes on them (wantchangeedges=false)
*/
gsl_matrix* BROWNIE::GetVCVforChangeNoChange(nxsstring chosentaxset, bool wantchangeedges)
{
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    Tree t = intrees.GetIthTree(chosentree-1);
    gsl_matrix *VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    int rowcount=-1;
    nxsstring rtaxon;
    nxsstring ctaxon;
IntSet::const_iterator ri;

    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=-1;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
IntSet::const_iterator ci;
        for( ci = taxonlistc.begin(); ci != taxonlistc.end(); ci++ ) {
            colcount++;
            ctaxon=taxa->GetTaxonLabel(*ci);
            nxsstring ctaxonUnderscores=blanks_to_underscores(ctaxon);
            nxsstring ctaxonBlanks=underscores_to_blanks(ctaxon);
            Node *cleafptr;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            bool goodmatch=false;
            while (q)
            {
                if (q->IsLeaf())
                {
                    string qlabel;
                    qlabel=q->GetLabel();
                    // nxsstring qlabel2<<qlabel;
                    //qlabel2=blanks_to_underscores(qlabel2);
                    if (qlabel==ctaxon.c_str() || qlabel==ctaxonUnderscores.c_str() || qlabel==ctaxonBlanks.c_str()) {
                        cleafptr=q;
                        goodmatch=true;
                    }
                }
                q = n.next();
            }
            if (goodmatch==false) {
                errormsg= "Error: there was trouble identifying taxon ";
                errormsg+=ctaxon.c_str();
                errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
                throw XNexus (errormsg );
            }
            double pathlength=0;
            Node *a;
            a=rleafptr;
            if (debugmode) {
                cout<<"Taxa "<<rleafptr->GetLabel()<<" (rtaxon="<<rtaxon<<") and "<<cleafptr->GetLabel()<<" (ctaxon="<<ctaxon<<"): ";
            }
            if (cleafptr==rleafptr) { //We need the pendant edge length
                                      // float unweightedpathlength=cleafptr->GetEdgeLength();
                                      //gsl_vector *modelcategoryvector;
                                      //modelcategoryvector=gsl_vector_calloc(cleafptr->GetModelCategory());
                vector<double> modelcategoryvector(cleafptr->GetModelCategory());
                // if (modelcategoryvector->size==0) {
                //     pathlength+=0;
                //}
                //else if(modelcategoryvector->size<(selectedmodel-1)) {
                //			errormsg= "Error: ModelCategory vector too short for the number of submodels ";
                //			throw XNexus (errormsg );
                //		}
                //		else {
                int numberofnonzeroentries=0;
                double temppathlength=0;
                for (int modelcat=0;modelcat<modelcategoryvector.size();modelcat++) {
                    if (modelcategoryvector[modelcat]>0) {
                        temppathlength+=modelcategoryvector[modelcat];
                        numberofnonzeroentries++;
                    }
                }
                if (wantchangeedges && (numberofnonzeroentries>1)) {
                    pathlength+=temppathlength;
                }
                if (!wantchangeedges && (numberofnonzeroentries<2)) {
                    pathlength+=temppathlength;
                }

                //		}
            }
            bool mrcanotfound=true;
            if (debugmode) {
                cout<<pathlength;
            }
            while (a->GetAnc() && mrcanotfound) { //General strategy here:
                a=a->GetAnc();
                Node *b;
                b=cleafptr;
                while (b->GetAnc() && mrcanotfound) {
                    b=b->GetAnc();
                    if (a==b) {
                        Node* mrca;
                        mrca=a;
                        mrcanotfound=false;
                        while(mrca->GetAnc()) {
                            vector<double> modelcategoryvector(mrca->GetModelCategory());
                            int numberofnonzeroentries=0;
                            double temppathlength=0;
                            for (int modelcat=0;modelcat<modelcategoryvector.size();modelcat++) {
                                if (modelcategoryvector[modelcat]>0) {
                                    temppathlength+=modelcategoryvector[modelcat];
                                    numberofnonzeroentries++;
                                }
                            }
                            if (wantchangeedges && (numberofnonzeroentries>1)) {
                                pathlength+=temppathlength;
                                if (debugmode) {
                                    cout<<" + "<<temppathlength;
                                }
                            }
                            if (!wantchangeedges && (numberofnonzeroentries<2)) {
                                pathlength+=temppathlength;
                                if (debugmode) {
                                    cout<<" + "<<temppathlength;
                                }
                            }
                            mrca=mrca->GetAnc();
                        }
                    }
                }
            }
            gsl_matrix_set(VCV,rowcount,colcount,pathlength);
            if (debugmode) {
                cout<<" = "<<pathlength<<endl;
            }
        }
    }
//matrixsingular=TestSingularity(VCV);
return VCV;
}

void BROWNIE::PrintMatrix(gsl_matrix *VCV)
{
	int ntax=VCV->size1;
	message="";
    for (int r=0;r<ntax;r++) {
        for (int c=0;c<ntax;c++) {
            message+=gsl_matrix_get(VCV,r,c);
			message+="\t";
        }
		message+="\n";
    }
    PrintMessage();
	
}




//Use's Pagel Delta transform for continuous characters. First, run DeleteStem on input matrix
gsl_matrix* BROWNIE::ConvertVCVwithDelta(gsl_matrix *VCVorig,double delta)
{
    int ntax=VCVorig->size1;
    gsl_matrix* VCVfinal=gsl_matrix_calloc(ntax,ntax);
    for (int r=0;r<ntax;r++) {
        for (int c=0;c<ntax;c++) {
            gsl_matrix_set(VCVfinal,r,c,(pow(gsl_matrix_get(VCVorig,r,c),delta)));
        }
    }
	if (debugmode) {
		PrintMatrix(VCVorig);
		cout<<"\nAfter transformation with delta of "<<delta<<endl<<endl;
		PrintMatrix(VCVfinal);
	}
	//cout<<gsl_matrix_get(VCVorig,0,0)<<"^"<<delta<<" = (using pow) "<<pow(gsl_matrix_get(VCVorig,0,0),delta)<<" and from output "<<gsl_matrix_get(VCVfinal,0,0)<<endl;
    return VCVfinal;
}

//Use's Pagel Lambda (continuous char) transform for continuous characters. First, run DeleteStem on input matrix
gsl_matrix* BROWNIE::ConvertVCVwithLambda(gsl_matrix *VCVorig,double lambda)
{
    int ntax=VCVorig->size1;
    gsl_matrix* VCVfinal=gsl_matrix_calloc(ntax,ntax);
    for (int r=0;r<ntax;r++) {
        for (int c=0;c<ntax;c++) {
			if (r!=c) {
            gsl_matrix_set(VCVfinal,r,c,(lambda*(gsl_matrix_get(VCVorig,r,c))));
			}
        }
    }
    return VCVfinal;
}

gsl_matrix* BROWNIE::DeleteStem(gsl_matrix *VCVorig)
{
    int ntax=VCVorig->size1;
    double minvalue=gsl_matrix_min(VCVorig);
    gsl_matrix* VCVfinal=gsl_matrix_calloc(ntax,ntax);
    for (int r=0;r<ntax;r++) {
        for (int c=0;c<ntax;c++) {
            gsl_matrix_set(VCVfinal,r,c,((gsl_matrix_get(VCVorig,r,c))-minvalue));
        }
    }
    return VCVfinal;
}

//This takes a matrix (VCV times rate) and adds a vector of known tip variances to it (by, basically, multiplying the identity matrix by the tip variance vector and adding the resulting matrix element by element to the VCV times rate matrix)
gsl_matrix* BROWNIE::AddTipVarianceVectorToRateTimesVCV(gsl_matrix * VCVorig,gsl_vector * TipVariance)
{
    int ntax=VCVorig->size1;
    gsl_matrix* VCVfinal=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix_memcpy (VCVfinal, VCVorig);
    for (int r=0;r<ntax;r++) {
        gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(VCVfinal,r,r)+(gsl_vector_get(TipVariance,r))));
    }
    return VCVfinal;
}






void BROWNIE::HandleDebugOptimization( NexusToken& token )
{
	/////////////ADD STUFF FOR TREE LOOP, SAVING FILES, USING OU MODELS
    ofstream tablef;
    nxsstring tablefname;
    bool tablef_open=false;
    bool name_provided=false;
    nxsstring chosentaxset;
    bool treeloop=false;
    bool adequateinput=false;
    int ntax=0;
    int nbest=1;
    int cstart=0;
    int cend=-1;
    nxsstring treefilename="";
    bool definedfilename=false;
    nxsstring tmessage;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (adequateinput==false) {
                message="Insufficient input: type \"optimize ?\" for help";
                PrintMessage();
            }
            break;
        }
        else if (token.Abbreviation("TReeloop") ) {
            nxsstring yesnotreeloop=GetFileName(token);
            if (yesnotreeloop[0] == 'n') {
                treeloop=false;
            }
            else {
                treeloop=true;
            }
        }
        else if( token.Abbreviation("FIle") ) {
            tablefname = GetFileName(token);
            name_provided = true;
            bool exists = FileExists( tablefname.c_str() );
            bool userok = true;
            if( exists && !UserSaysOk( "Ok to replace?", "Optimize output file specified already exists" ) )
                userok = false;
            if( userok ) {
                tablef_open = true;
                tablef.open( tablefname.c_str() );
            }
			
            if( exists && userok ) {
                message = "\nReplacing optimize output file ";
                message += tablefname;
            }
            else if( userok ) {
                message = "\nOptimize output file ";
                message += tablefname;
                message += " opened";
            }
            else {
                errormsg = "Aborting the optimization so as not to overwrite the file.\n";
                throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
				
            }
            PrintMessage();
			
        }
        else if( token.Abbreviation("Nbest") ) {
            
            nxsstring numbernexus = GetNumber(token);
            nbest=atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen to save ";
            message+=nbest;
            message+=" trees";
            PrintMessage();
            if (nbest<1) {
                errormsg = "Error: must select a number greater than zero";
                nbest=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("CStart") ) {
			
            nxsstring numbernexus = GetNumber(token);
            cstart=atoi( numbernexus.c_str() ); //convert to int
        }
        else if( token.Abbreviation("CEnd") ) {
			
            nxsstring numbernexus = GetNumber(token);
            cend=atoi( numbernexus.c_str() ); //convert to int
        }
        
		
		
        else if( token.Equals("?") ) {
            message="Usage: Optimize taxset=<chosen taxset> treeloop=[yes|no]\n\n";
            message+="Returns the likelihood and AIC under the current model.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\nTaxset       <taxset name>                        *None";
            message+="\nTreeloop     No|Yes                               *No";
            message+="\nFile         <file name>                          *None";
            message+="\nNBest        <integer>                            *1";
            message+="\nCStart       <integer>                            *none";
            message+="\nCEnd         <integer>                            *none";
            message+="\n                                                 *Option is nonpersistent\n\n";
            PrintMessage();
        }
        else if( token.Abbreviation("TReefile") ) {
            treefilename=GetFileName(token);
            definedfilename=true;
            break;
        }
        else if (token.Abbreviation("TAxset") ) {
            adequateinput=true;
            if (trees->GetNumTrees()<1) {
                errormsg = "Error: No valid trees are loaded.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            chosentaxset=GetFileName(token);
            IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
            if (taxonlist.empty()) {
                errormsg= "Error: Taxset ";
                errormsg+=chosentaxset.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
			IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntax++;
            }
        }
    }
    if (adequateinput) {
        int starttree, stoptree, startchar, stopchar;
        int originalchosentree=chosentree;
        int originalchosenchar=chosenchar;
        if (treeloop) {
            starttree=1;
            stoptree=trees->GetNumTrees();
        }
        else {
            starttree=chosentree;
            stoptree=chosentree;
        }
        for (chosentree=starttree;chosentree<=stoptree;chosentree++) {
            if (tablef_open) {
                tmessage="\n";
                tmessage+=chosentree;
                tablef<<tmessage;
                message="Now working on tree number ";
                message+=chosentree;
                PrintMessage();
            }
            gsl_matrix * VCV=gsl_matrix_calloc(ntax,ntax);
            gsl_vector * tips=gsl_vector_calloc(ntax);
            gsl_vector * variance=gsl_vector_calloc(ntax);
            tips=GetTipValues(chosentaxset,chosenchar);
            if (tipvariancetype==2) {
                variance=GetTipValues(chosentaxset,chosenchar+1);
            }
            if (chosenmodel<5 || chosenmodel==21 || chosenmodel==22) {
                VCV=DeleteStem(GetVCV(chosentaxset));
                OptimizationFn my_fn(VCV,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);
				
                if (tipvariancetype==1 && chosenmodel==2) {
                    gsl_vector *optimalrate=gsl_vector_calloc(5);
					//gsl_vector_memcpy(optimalrate,my_fn.OptimizeRateWithOptimizedTipVariance());
                    gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(2));
                    message="\nOptimal rate = ";
                    message+=gsl_vector_get(optimalrate,0);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,2);
                    message+="\nTip variance = ";
                    message+=gsl_vector_get(optimalrate,1);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,3);
                    message+="\nlnL = ";
                    message+=gsl_vector_get(optimalrate,4);
                    if (tablef_open) {
                        tmessage="\t";
                        tmessage+=gsl_vector_get(optimalrate,0);
                        tmessage+="\t";
                        tmessage+=gsl_vector_get(optimalrate,2);
                        tmessage+="\t";
                        tmessage+=gsl_vector_get(optimalrate,1);
                        tmessage+="\t";
                        tmessage+=gsl_vector_get(optimalrate,3);
                        tmessage+="\t";
                        tmessage+=gsl_vector_get(optimalrate,4);
                        tablef<<tmessage;
                    }
					gsl_vector_free(optimalrate);
				}
                else if (chosenmodel==1) {
                    gsl_vector *optimalrate=gsl_vector_calloc(3);
					//gsl_vector_memcpy(optimalrate,my_fn.OptimizeRateWithGivenTipVariance());
                    gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(1));
                    message="\nOptimal rate = ";
                    message+=gsl_vector_get(optimalrate,0);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,1);
                    message+="\nlnL = ";
                    message+=gsl_vector_get(optimalrate,2);
					gsl_vector_free(optimalrate);
                }
                else if (chosenmodel==3) {
                    gsl_vector *optimalrate=gsl_vector_calloc(7);
                    gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(3));
                    message="\nOptimal rate = ";
                    message+=gsl_vector_get(optimalrate,0);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,3);
                    message+="\nAncestral state = ";
                    message+=gsl_vector_get(optimalrate,1);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,4);
                    message+="\nd = ";
                    message+=gsl_vector_get(optimalrate,2);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,5);
                    message+="\nlnL = ";
                    message+=gsl_vector_get(optimalrate,6);
					gsl_vector_free(optimalrate);
				}
                else if (chosenmodel==4) {
                    gsl_vector *optimalrate=gsl_vector_calloc(7);
                    gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(4));
                    message="\nOptimal rate = ";
                    message+=gsl_vector_get(optimalrate,0);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,3);
                    message+="\nAncestral state = ";
                    message+=gsl_vector_get(optimalrate,1);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,4);
                    message+="\ng = ";
                    message+=gsl_vector_get(optimalrate,2);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,5);
                    message+="\nlnL = ";
                    message+=gsl_vector_get(optimalrate,6);
					gsl_vector_free(optimalrate);
                }
				else if (chosenmodel==21) {
                    gsl_vector *optimalrate=gsl_vector_calloc(7);
                    gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(21));
                    message="\nOptimal rate = ";
                    message+=gsl_vector_get(optimalrate,0);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,3);
                    message+="\nAncestral state = ";
                    message+=gsl_vector_get(optimalrate,1);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,4);
                    message+="\ndelta = ";
                    message+=gsl_vector_get(optimalrate,2);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,5);
                    message+="\nlnL = ";
                    message+=gsl_vector_get(optimalrate,6);
					gsl_vector_free(optimalrate);
				}
				else if (chosenmodel==22) {
                    gsl_vector *optimalrate=gsl_vector_calloc(7);
                    gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(22));
                    message="\nOptimal rate = ";
                    message+=gsl_vector_get(optimalrate,0);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,3);
                    message+="\nAncestral state = ";
                    message+=gsl_vector_get(optimalrate,1);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,4);
                    message+="\nlambda = ";
                    message+=gsl_vector_get(optimalrate,2);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalrate,5);
                    message+="\nlnL = ";
                    message+=gsl_vector_get(optimalrate,6);
					gsl_vector_free(optimalrate);
				}
                message+="\n\nNote that +/- reflects imprecision due to numerical optimization";
                PrintMessage();
            }
            else { //we must have a model that deals with multiple VCVs
                gsl_matrix * VCV0=gsl_matrix_calloc(ntax,ntax);
                gsl_matrix * VCV1=gsl_matrix_calloc(ntax,ntax);
                gsl_matrix * VCV2=gsl_matrix_calloc(ntax,ntax);
                gsl_matrix * VCV3=gsl_matrix_calloc(ntax,ntax);
                gsl_matrix * VCV4=gsl_matrix_calloc(ntax,ntax);
                gsl_matrix * VCV5=gsl_matrix_calloc(ntax,ntax);
                gsl_matrix * VCV6=gsl_matrix_calloc(ntax,ntax);
                gsl_matrix * VCV7=gsl_matrix_calloc(ntax,ntax);
                gsl_matrix * VCV8=gsl_matrix_calloc(ntax,ntax);
                gsl_matrix * VCV9=gsl_matrix_calloc(ntax,ntax);
                if (chosenmodel==5) {
                    VCV0=GetVCVforOneModel(chosentaxset,0);
                    //cout<<"VCV0 "<<gsl_matrix_get(VCV0,0,0);
                    VCV1=GetVCVforOneModel(chosentaxset,1);
                   // cout<<"VCV1 "<<gsl_matrix_get(VCV1,0,0);
                    VCV2=GetVCVforOneModel(chosentaxset,2);
                   // cout<<"VCV2 "<<gsl_matrix_get(VCV2,0,0);
                    VCV3=GetVCVforOneModel(chosentaxset,3);
                    VCV4=GetVCVforOneModel(chosentaxset,4);
                    VCV5=GetVCVforOneModel(chosentaxset,5);
                    VCV6=GetVCVforOneModel(chosentaxset,6);
                    VCV7=GetVCVforOneModel(chosentaxset,7);
                    VCV8=GetVCVforOneModel(chosentaxset,8);
                    VCV9=GetVCVforOneModel(chosentaxset,9);
                }
                else if (chosenmodel==6) {
                    VCV0=GetVCVforChangeNoChange(chosentaxset,false); //branches with no changes
                    VCV1=GetVCVforChangeNoChange(chosentaxset,true); //branches with changes
                }
                if (debugmode) {
                    cout<<"\n\nVCV0\n";
                    for (int rt=0;rt<ntax;rt++) {
                        cout<<"\n";
                        for (int ct=0;ct<ntax;ct++) {
                            cout<<gsl_matrix_get(VCV0,rt,ct)<<"\t";
                        }
                    }
                    cout<<"\n\nVCV1\n";
                    for (int rt=0;rt<ntax;rt++) {
                        cout<<"\n";
                        for (int ct=0;ct<ntax;ct++) {
                            cout<<gsl_matrix_get(VCV1,rt,ct)<<"\t";
                        }
                    }
                    cout<<"\n\nVCV2\n";
                    for (int rt=0;rt<ntax;rt++) {
                        cout<<"\n";
                        for (int ct=0;ct<ntax;ct++) {
                            cout<<gsl_matrix_get(VCV2,rt,ct)<<"\t";
                        }
                    }
                    cout<<"\n\nVCV3\n";
                    for (int rt=0;rt<ntax;rt++) {
                        cout<<"\n";
                        for (int ct=0;ct<ntax;ct++) {
                            cout<<gsl_matrix_get(VCV3,rt,ct)<<"\t";
                        }
                    }
                    cout<<"\n\nVCV4\n";
                    for (int rt=0;rt<ntax;rt++) {
                        cout<<"\n";
                        for (int ct=0;ct<ntax;ct++) {
                            cout<<gsl_matrix_get(VCV4,rt,ct)<<"\t";
                        }
                    }
                    cout<<"\n\nVCV5\n";
                    for (int rt=0;rt<ntax;rt++) {
                        cout<<"\n";
                        for (int ct=0;ct<ntax;ct++) {
                            cout<<gsl_matrix_get(VCV5,rt,ct)<<"\t";
                        }
                    }
                    cout<<"\n\nVCV6\n";
                    for (int rt=0;rt<ntax;rt++) {
                        cout<<"\n";
                        for (int ct=0;ct<ntax;ct++) {
                            cout<<gsl_matrix_get(VCV6,rt,ct)<<"\t";
                        }
                    }
                    cout<<"\n\nVCV7\n";
                    for (int rt=0;rt<ntax;rt++) {
                        cout<<"\n";
                        for (int ct=0;ct<ntax;ct++) {
                            cout<<gsl_matrix_get(VCV7,rt,ct)<<"\t";
                        }
                    }
                    cout<<"\n\nVCV8\n";
                    for (int rt=0;rt<ntax;rt++) {
                        cout<<"\n";
                        for (int ct=0;ct<ntax;ct++) {
                            cout<<gsl_matrix_get(VCV8,rt,ct)<<"\t";
                        }
                    }
                    cout<<"\n\nVCV9\n";
                    for (int rt=0;rt<ntax;rt++) {
                        cout<<"\n";
                        for (int ct=0;ct<ntax;ct++) {
                            cout<<gsl_matrix_get(VCV9,rt,ct)<<"\t";
                        }
                    }
					
                }
                if (chosenmodel==5) {
                    OptimizationFnMultiModel my_fn(VCV0,VCV1,VCV2,VCV3,VCV4,VCV5,VCV6,VCV7,VCV8,VCV9,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);
                    gsl_vector *optimalvalues=gsl_vector_calloc(24);
                    gsl_vector_memcpy(optimalvalues,my_fn.GeneralOptimization(5));
					message="\nlnL = ";
					message+=gsl_vector_get(optimalvalues,0);
					message+="\nAIC = ";
					message+=2*(gsl_vector_get(optimalvalues,0) + gsl_vector_get(optimalvalues,1));
					message+="\nAICc = ";
					message+=(2*gsl_vector_get(optimalvalues,0))+2*gsl_vector_get(optimalvalues,1)+2*gsl_vector_get(optimalvalues,1)*(gsl_vector_get(optimalvalues,1)+1)/(ntax-gsl_vector_get(optimalvalues,1)-1);
                    message+="\nAncestral state = ";
                    message+=gsl_vector_get(optimalvalues,2);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalvalues,13);
                    int np=int(gsl_vector_get(optimalvalues,1));
                    for (int modelstate=0;modelstate<(np-1);modelstate++) {
                        message+="\nRate in state ";
                        message+=modelstate;
                        message+=" = ";
                        message+=gsl_vector_get(optimalvalues,3+modelstate);
                        message+=" +/- ";
                        message+=gsl_vector_get(optimalvalues,14+modelstate);
                    }
                    PrintMessage();
					gsl_vector_free(optimalvalues);
				}
                if (chosenmodel==12) {
                    VCV0=DeleteStem(GetVCV(chosentaxset));
                    gsl_matrix_free (VCV1);
                    VCV1=gsl_matrix_calloc(ntax,6*ntax); //assumes that states change fewer than three times per branch on pectinate tree
                    VCV1=GetStartStopTimesforOneState(chosentaxset,0);
                    gsl_matrix_free (VCV2);
                    VCV2=gsl_matrix_calloc(ntax,6*ntax);
                    VCV2=GetStartStopTimesforOneState(chosentaxset,1);
                    gsl_matrix_free (VCV3);
                    VCV3=gsl_matrix_calloc(ntax,6*ntax);
                    VCV3=GetStartStopTimesforOneState(chosentaxset,2);
                    gsl_matrix_free (VCV4);
                    VCV4=gsl_matrix_calloc(ntax,6*ntax);
                    VCV4=GetStartStopTimesforOneState(chosentaxset,3);
                    gsl_matrix_free (VCV5);
                    VCV5=gsl_matrix_calloc(ntax,6*ntax);
                    VCV5=GetStartStopTimesforOneState(chosentaxset,4);
                    gsl_matrix_free (VCV6);
                    VCV6=gsl_matrix_calloc(ntax,6*ntax);
                    VCV6=GetStartStopTimesforOneState(chosentaxset,5);
                    gsl_matrix_free (VCV7);
                    VCV7=gsl_matrix_calloc(ntax,6*ntax);
                    VCV7=GetStartStopTimesforOneState(chosentaxset,6);
                    gsl_matrix_free (VCV8);
                    VCV8=gsl_matrix_calloc(ntax,6*ntax);
                    VCV8=GetStartStopTimesforOneState(chosentaxset,7);
                    gsl_matrix_free (VCV9);
                    VCV9=gsl_matrix_calloc(ntax,6*ntax);
                    VCV9=GetStartStopTimesforOneState(chosentaxset,8);
                    OptimizationFnMultiModel my_fn(VCV0,VCV1,VCV2,VCV3,VCV4,VCV5,VCV6,VCV7,VCV8,VCV9,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);
                    gsl_vector *optimalvalues=gsl_vector_calloc(28);
                    gsl_vector_memcpy(optimalvalues,my_fn.GeneralOptimization(12));
                    message="\nlnL = ";
                    message+=gsl_vector_get(optimalvalues,0);
                    message+="\nBM rate = ";
                    message+=gsl_vector_get(optimalvalues,2);
                   // message+=" +/- ";
                    //message+=gsl_vector_get(optimalvalues,15);
                    message+="\nOU attraction = ";
                    message+=gsl_vector_get(optimalvalues,3);
					//message+=" +/- ";
					//message+=gsl_vector_get(optimalvalues,16);
                    message+="\nRoot state = ";
                    message+=gsl_vector_get(optimalvalues,4);
					//message+=" +/- ";
					//message+=gsl_vector_get(optimalvalues,15);
                    int np=int(gsl_vector_get(optimalvalues,1));
                    for (int modelstate=0;modelstate<(np-3);modelstate++) {
                        message+="\nMean value in state ";
                        message+=modelstate;
                        message+=" = ";
                        message+=gsl_vector_get(optimalvalues,5+modelstate);
                        //message+=" +/- ";
                        //message+=gsl_vector_get(optimalvalues,18+modelstate);
                    }
                    PrintMessage();
					gsl_matrix_free(VCV1);
					gsl_matrix_free(VCV2);
					gsl_matrix_free(VCV3);
					gsl_matrix_free(VCV4);
					gsl_matrix_free(VCV5);
					gsl_matrix_free(VCV6);
					gsl_matrix_free(VCV7);
					gsl_matrix_free(VCV8);
					gsl_matrix_free(VCV9);
					
					
					gsl_vector_free(optimalvalues);
				}
                if (chosenmodel==14) {
                    Tree t = intrees.GetIthTree(chosentree-1);
                    double MaxRootTipLength=t.GetMaxPathLength();
                    if (debugmode) {
                        cout<<"MaxRootTipLength = "<<MaxRootTipLength<<endl;
                    }
                }
                if (chosenmodel==16 || chosenmodel==17) {
                    if (chosenmodel==16) {
                        treefilename="OutputTreesBMAN.nex";
                    }
                    else {
                        treefilename="OutputTreesBMAO.nex";
                    }
                    ofstream outtreef;
                    outtreef.open(treefilename.c_str());
                    outtreef<<"#nexus\nbegin trees;";
                    
                    gsl_combination *c;
                    size_t combinationsize;
                    //int ntax=(intrees.GetIthTree(0)).GetNumLeaves();
                    int totalnumberofcombinations=0;
                    for (combinationsize=0;combinationsize<(2*ntax-2);combinationsize++) {
                        totalnumberofcombinations+=int(gsl_sf_choose ((2*ntax-2), combinationsize));
                    }
                    message="There are ";
                    message+=totalnumberofcombinations;
                    message+=" total combinations to try.";
                    PrintMessage();
                    ProgressBar(totalnumberofcombinations);
                    int combinationcount=0;
                    for (combinationsize=0;combinationsize<(2*ntax-2);combinationsize++) { //we want to try all possible assignments of zero, except assigning all branches to zero
                        c=gsl_combination_calloc((2*ntax-2),combinationsize);
                        do
                        {
                            combinationcount++;
                            vector<int> nodestomakezero;
                            for (int i=0;i<combinationsize;i++) {
                                nodestomakezero.push_back(int(gsl_combination_get(c,i)));
                            }
                            nodestomakezero.push_back(-1); //just to keep things from going awry when we're done with all the nodes
                            int vectorstep=0;
                            int nodecount=0;
                            Tree t=intrees.GetIthTree(chosentree-1);
                            t.SetEdgeLengths(true);
                            PreorderIterator <Node> n (t.GetRoot());
                            NodePtr RootNode=t.GetRoot();
                            NodePtr currentnode = n.begin();
                            while (currentnode) {
                                if (currentnode!=RootNode) {
                                    if (nodecount==nodestomakezero[vectorstep]) {
                                        currentnode->SetEdgeLength(0.0);
                                        vectorstep++;
                                    }
                                    else if (chosenmodel==17) {
                                        currentnode->SetEdgeLength(1.0);
                                    }
                                    nodecount++; //don't iterate on the root node, as it's not in the combination
                                }
                                currentnode=n.next();    
                            }
                            gsl_matrix *currentVCVmat=gsl_matrix_calloc(ntax,ntax);
                            currentVCVmat=DeleteStem(GetVCVwithTree(chosentaxset,t));
                            gsl_vector *tips=gsl_vector_calloc(ntax);
                            gsl_vector *tipsresid=gsl_vector_calloc(ntax);
                            double ancstate;
                            double rate;
                            double likelihood;
                            tips=GetTipValues(chosentaxset,chosenchar);
                            ancstate=GetAncestralState(currentVCVmat,tips);
                            tipsresid=GetTipResiduals(tips,ancstate);
                            rate=EstimateRate(currentVCVmat,tipsresid);
                            gsl_matrix *RateTimesVCVfortest=gsl_matrix_calloc(currentVCVmat->size1,currentVCVmat->size2);
                            gsl_matrix_memcpy(RateTimesVCVfortest, currentVCVmat);
                            gsl_matrix_scale(RateTimesVCVfortest,rate);
                            likelihood=GetLScore(currentVCVmat,tipsresid,rate);
                            nxsstring rootlabel="lnL_";
                            rootlabel+=likelihood;
                            RootNode->SetLabel(rootlabel);
                            t.SetInternalLabels(true);
                            t.Update();
                           // outtreef<<"gethasedgelengths = "<<t.GetHasEdgeLengths()<<endl;
                            char outputstring[14];
                            sprintf(outputstring,"%14.6f",likelihood);
                            //outtreef<<endl<<"gethasedgelengths = "<<t.GetHasEdgeLengths();
							
                            outtreef<<"\ntree tree"<<combinationcount<<" = [likelihood = "<<outputstring<<" ] ";
                            NewickTreeWriter tw (&t);
                            tw.SetStream (&outtreef);
                            tw.SetWriteEdgeLengths(true);
                            tw.Write();
                            ProgressBar(0);
							gsl_matrix_free(currentVCVmat);
							gsl_vector_free(tips);
							gsl_vector_free(tipsresid);
							gsl_matrix_free(RateTimesVCVfortest);
							
						}
                        while (gsl_combination_next(c) == GSL_SUCCESS);
                        gsl_combination_free (c);
                    }
                    outtreef<<"\nend;";
                    outtreef.close();
                }
				if (chosenmodel==18 || chosenmodel==19) {
                    
                     //vector<Tree> BestTrees;
                     //vector<double> ScoresOfBestTrees;
					double bestscore=GSL_POSINF;
                     //cout<<"cstart = "<<cstart<<" cend = "<<cend<<" gslposinf = "<<int(GSL_POSINF)<<endl;
                    if (chosenmodel==18 && !definedfilename) {
                        treefilename="BMPN";
                        if ((cstart>0) && (cend>0)) {
                            treefilename+="_";
                            treefilename+=cstart;
                            treefilename+="_";
                            treefilename+=cend;
                            treefilename+=".tre";
                        }
                        else {
                            treefilename+="_ALL.tre";
                        }
                    }
                    else if (chosenmodel==19 && !definedfilename) {
                        treefilename="BMPO";
                        if ((cstart>0) && (cend>0)) {
                            treefilename+="_";
                            treefilename+=cstart;
                            treefilename+="_";
                            treefilename+=cend;
                            treefilename+=".tre";
                        }
                        else {
                            treefilename+="_ALL.tre";
                        }
                        
                    }
                    ofstream outtreef;
                    outtreef.open(treefilename.c_str());
                   // if (outtreef.fail()) {
                   //     cout<<"Has trouble writing to file!!!!";
                 //   }
					
                    outtreef<<"#nexus\nbegin trees;";
                    outtreef.close();
                    gsl_combination *c;
                    size_t combinationsize;
                    //int ntax=(intrees.GetIthTree(0)).GetNumLeaves();
                    int totalnumberofcombinations=0;
                    for (combinationsize=0;combinationsize<=(ntax-2);combinationsize++) { //ntax-1 internal nodes, minus the root (reml problem)
                        totalnumberofcombinations+=int(gsl_sf_choose ((ntax-2), combinationsize));
                    }
                    message="There are ";
                    message+=totalnumberofcombinations;
                    message+=" total combinations to try, you are doing ";
                    if (cend<0) {
                        cend=totalnumberofcombinations+1;
                    }
                    //cout<<"cstart = "<<cstart<<" cend = "<<cend<<" gslposinf = "<<int(GSL_POSINF)<<endl;
                    int actualcombinations=GSL_MIN(totalnumberofcombinations,cend-cstart+1);
                    message+=actualcombinations;
                    PrintMessage();
                    int loopcounter=0;
                    ProgressBar(actualcombinations);
                    int combinationcount=0;
                    //int cutoff=int(floor(combinationcount/10.0));
                    int cutoff=500;
                    message="\nComb.\t-lnL";
                    PrintMessage();
                    for (combinationsize=0;combinationsize<=(ntax-2);combinationsize++) { //if a node is selected, its child has brlen of zero; otherwise, its child's sibling does
                        c=gsl_combination_calloc((ntax-2),combinationsize);
                        do
                        {
                            combinationcount++;
                            if (combinationcount>=cstart && combinationcount<=cend) {
                                loopcounter++;
                                ProgressBar(0);
                                vector<int> nodestomakezero;
                                for (int i=0;i<combinationsize;i++) {
                                    nodestomakezero.push_back(int(gsl_combination_get(c,i)));
                                }
                                nodestomakezero.push_back(-1); //just to keep things from going awry when we're done with all the nodes
                                int vectorstep=0;
                                int nodecount=0;
                                Tree t=intrees.GetIthTree(chosentree-1);
                                t.SetEdgeLengths(true);
                                t.SetInternalLabels(true);
                                PreorderIterator <Node> n (t.GetRoot());
                                NodePtr RootNode=t.GetRoot();
                                NodePtr currentnode = n.begin();
                                while (currentnode) {
                                    if (!(currentnode->IsLeaf()) && (currentnode!=RootNode)) {
                                        if (nodecount==nodestomakezero[vectorstep]) {
                                            (currentnode->GetChild())->SetEdgeLength(0.0);
                                            if (chosenmodel==19) {
                                                ((currentnode->GetChild())->GetSibling())->SetEdgeLength(1.0);
                                            }
                                            vectorstep++;
                                        }
                                        else {
                                            ((currentnode->GetChild())->GetSibling())->SetEdgeLength(0.0);
                                            if (chosenmodel==19) {
                                                (currentnode->GetChild())->SetEdgeLength(1.0);
                                            }
                                        }
                                        nodecount++;
                                    }
                                    else if((currentnode==RootNode) && (chosenmodel==19)) {
                                        (currentnode->GetChild())->SetEdgeLength(1.0);
                                        ((currentnode->GetChild())->GetSibling())->SetEdgeLength(1.0);
                                    }
                                    currentnode=n.next();
                                }
                                gsl_matrix *currentVCVmat=gsl_matrix_calloc(ntax,ntax);
                                currentVCVmat=DeleteStem(GetVCVwithTree(chosentaxset,t));
                                gsl_vector *tips=gsl_vector_calloc(ntax);
                                gsl_vector *tipsresid=gsl_vector_calloc(ntax);
                                double ancstate;
                                double rate;
                                double likelihood;
                                tips=GetTipValues(chosentaxset,chosenchar);
                                ancstate=GetAncestralState(currentVCVmat,tips);
                                tipsresid=GetTipResiduals(tips,ancstate);
                                rate=EstimateRate(currentVCVmat,tipsresid);
                                gsl_matrix *RateTimesVCVfortest=gsl_matrix_calloc(currentVCVmat->size1,currentVCVmat->size2);
                                gsl_matrix_memcpy(RateTimesVCVfortest, currentVCVmat);
                                gsl_matrix_scale(RateTimesVCVfortest,rate);
                                likelihood=GetLScore(currentVCVmat,tipsresid,rate);
                                //t.Draw(cout);
                                //cout<<"likelihood is "<<likelihood<<endl;
                                if (likelihood==likelihood) { //test for NaN
                                                              // if (likelihood<=bestscore) {
                                                              // if(BestTrees.size()<nbest) {
                                                              //       BestTrees.push_back(t);
                                                              //      ScoresOfBestTrees.push_back(likelihood);
                                                              //  }
                                                              //   else {
                                                              //       for (int i=0;i<nbest-1;i++) {
                                                              //           BestTrees[i]=BestTrees[i+1];
                                                              //           ScoresOfBestTrees[i]=ScoresOfBestTrees[i+1]; //move everything over one, except for the last cell, which we'll replace
                                                              //       }
                                                              //       BestTrees.pop_back();
                                                              //       BestTrees.push_back(t);
                                                              //       ScoresOfBestTrees.pop_back();
                                                              //       ScoresOfBestTrees.push_back(likelihood);
                                                              //   }
                                                              // if (likelihood<bestscore) { //overwrite
                                                              //      outtreef.open(treefilename.c_str());
                                                              //      outtreef<<"#nexus\nbegin trees;";
                                                              // }
                                                              // else {
                                                              //save all trees
                                    outtreef.open(treefilename.c_str(), ios::out | ios::app );
                                    // }
                                    if (likelihood<=bestscore) {
                                        bestscore=likelihood;
                                    }
                                    nxsstring rootlabel="lnL_";
                                    rootlabel+=likelihood;
                                    RootNode->SetLabel(rootlabel);
                                    outtreef<<"\ntree tree"<<combinationcount<<" = [&R] [likelihood = "<<likelihood<<" ancstate = "<<ancstate<<" rate = "<<rate<<" ] ";
                                    t.Write(outtreef);
                                    outtreef.close();
									
															  }
                                if (loopcounter==cutoff) {
                                    message="";
                                    message+=combinationcount;
                                    message+="\t";
                                    message+=bestscore;
                                    PrintMessage();
                                    loopcounter=0;
                                }
								gsl_matrix_free(currentVCVmat);
								gsl_vector_free(tips);
								gsl_vector_free(tipsresid);
								gsl_matrix_free(RateTimesVCVfortest);
								}
							}
                        while ((gsl_combination_next(c) == GSL_SUCCESS) && (combinationcount<=cend));
                        gsl_combination_free (c);
						}
                   // cout<<"Saving "<<BestTrees.size()<<" best trees"<<endl;
                    //for (int i=0;i<BestTrees.size();i++) {
                    //    Tree goodtree=BestTrees.back();
                     //   goodtree.Update();
                     //   goodtree.Draw(cout);
                     //   BestTrees.pop_back();
                   //     double goodlikelihood=ScoresOfBestTrees.back();
                   //     ScoresOfBestTrees.pop_back();
                   //     cout<<"lnL = "<<goodlikelihood<<endl;
                   // outtreef<<"\ntree tree"<<i+1<<" = [&R] [likelihood = "<<goodlikelihood<<" ] ";
                   //  goodtree.Write(outtreef);
                   // }
                    outtreef.open(treefilename.c_str(), ios::out | ios::app );
                    outtreef<<"\nend;";
                    outtreef.close();
					}
				
                if (chosenmodel==6) {
                    OptimizationFnMultiModel my_fn(VCV0,VCV1,VCV2,VCV3,VCV4,VCV5,VCV6,VCV7,VCV8,VCV9,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);
                    gsl_vector *optimalvalues=gsl_vector_calloc(7);
                    gsl_vector_memcpy(optimalvalues,my_fn.GeneralOptimization(6));
                   	message="\nlnL = ";
					message+=gsl_vector_get(optimalvalues,0);
					message+="\nAIC = ";
					message+=2*(gsl_vector_get(optimalvalues,0) + gsl_vector_get(optimalvalues,1));
					message+="\nAICc = ";
					message+=(2*gsl_vector_get(optimalvalues,0))+2*gsl_vector_get(optimalvalues,1)+2*gsl_vector_get(optimalvalues,1)*(gsl_vector_get(optimalvalues,1)+1)/(ntax-gsl_vector_get(optimalvalues,1)-1);
                    message+="\nAncestral state = ";
                    message+=gsl_vector_get(optimalvalues,1);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalvalues,4);
                    message+="\nRate on branches with no changes = ";
                    message+=gsl_vector_get(optimalvalues,2);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalvalues,5);
                    message+="\nRate on branches with changes = ";
                    message+=gsl_vector_get(optimalvalues,3);
                    message+=" +/- ";
                    message+=gsl_vector_get(optimalvalues,6);
                    message+="\nlnL = ";
                    message+=gsl_vector_get(optimalvalues,0);
                    PrintMessage();
                    if (tablef_open) {
                        tmessage="\t";
                        tmessage+=gsl_vector_get(optimalvalues,0);
                        tmessage+="\t";
                        tmessage+=gsl_vector_get(optimalvalues,1);
                        tmessage+="\t";
                        tmessage+=gsl_vector_get(optimalvalues,4);
                        tmessage+="\t";
                        tmessage+=gsl_vector_get(optimalvalues,2);
                        tmessage+="\t";
                        tmessage+=gsl_vector_get(optimalvalues,5);
                        tmessage+="\t";
                        tmessage+=gsl_vector_get(optimalvalues,3);
                        tmessage+="\t";
                        tmessage+=gsl_vector_get(optimalvalues,6);
                        tablef<<tmessage;
                    }
					gsl_vector_free(optimalvalues);
				}
				gsl_matrix_free(VCV0);
				gsl_matrix_free(VCV1);
				gsl_matrix_free(VCV2);
				gsl_matrix_free(VCV3);
				gsl_matrix_free(VCV4);
				gsl_matrix_free(VCV5);
				gsl_matrix_free(VCV6);
				gsl_matrix_free(VCV7);
				gsl_matrix_free(VCV8);
				gsl_matrix_free(VCV9);
				
				}
			gsl_matrix_free(VCV);
			gsl_vector_free(tips);
			gsl_vector_free(variance);
			}
        chosentree=originalchosentree; //restore initial values.
        chosenchar=originalchosenchar;
        if (tablef_open) {
            tablef.close();
        }
		}
	}


//Returns the ancestral state, using formula from Martins and Lamont 1998, Animal Behavior 55: 1685-1706, bottom right of page 1689.
double BROWNIE::GetAncestralState(gsl_matrix *VCV, gsl_vector *tips)
{
    if (debugmode) {
        cout<<"\nNow in GetAncestralState\n";
    }
    int ntax=tips->size;
    if (debugmode) {
        message="\nNtax = ";
        message+=ntax;
        message+="\n";
        PrintMessage();
    }
    //cout<<"ntax="<<ntax<<endl;
    double ancestralstate;
    //matrixsingular=TestSingularity(VCV);
    // if (matrixsingular) {
    //     errormsg="Error: VCV matrix is singular in GetAncestralState ";
    //    throw XNexus( errormsg);
    // }
    if (debugmode) {
        message="\nNow calculating inverse VCV for estimating ancestral state.\n";
        message+="\n\nntax is ";
        message+=ntax;
        PrintMessage();
    }
    if(debugmode) {
        message="debugging line 3494\n";
        for (int currentrow=0;currentrow<VCV->size1;currentrow++) {
            for (int currentcol=0;currentcol<VCV->size2;currentcol++) {
                message+=gsl_matrix_get(VCV,currentrow,currentcol);
                message+="\t";
            }
            message+="\n";
        }
        PrintMessage();
    }

    gsl_permutation * p = gsl_permutation_alloc (ntax);
    int signum;
    gsl_matrix * VCVinverse=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix * VCVinversestart=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix_memcpy (VCVinversestart, VCV);
    if(debugmode) {
        message="debugging line 3511\n";
        for (int currentrow=0;currentrow<VCV->size1;currentrow++) {
            for (int currentcol=0;currentcol<VCV->size2;currentcol++) {
                message+=gsl_matrix_get(VCV,currentrow,currentcol);
                message+="\t";
            }
            message+="\n";
        }
        PrintMessage();
    }
    gsl_linalg_LU_decomp (VCVinversestart,p, &signum);
    gsl_linalg_LU_invert (VCVinversestart,p, VCVinverse);
    //VCVinverse=Inverse(VCV);
    if (debugmode) {
        message="\nFinished calculating inverse VCV for estimating ancestral state.\n";
        PrintMessage();
    }
    gsl_matrix * tipsasmatrix=gsl_matrix_calloc(ntax,1);
    for (int i=0; i<ntax; i++) {
        gsl_matrix_set(tipsasmatrix,i,0,gsl_vector_get(tips,i));
    }
    gsl_matrix * stepA=gsl_matrix_calloc(ntax,1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, VCVinverse, tipsasmatrix,0, stepA);
    //stepA=VCVinverse*tipsasmatrix;
    double stepB=0.0;
    for (int i=0; i<ntax; i++) {
        stepB+=gsl_matrix_get(stepA,i,0);
        //cout<<stepA[i][0]<<" stepB "<<stepB<<endl;

    }
    double stepC=0.0;
    for (int i=0; i<ntax; i++) {
        for (int j=0; j<ntax; j++) {
            stepC+=gsl_matrix_get(VCVinverse,i,j);
            // cout<<"StepC  "<<stepC<<endl;
        }
    }
    if (stepC==0) {
        errormsg="Error: Division by zero in GetAncestralState routine";
        throw XNexus(errormsg);
    }
    else {
        ancestralstate=stepB/stepC;
    }

    //gsl_matrix J(ntax,1,1);
    //gsl_matrix Jprime(1,ntax,1);
    //gsl_matrix step1(ntax,1,0);
    //   step1=VCVinverse*J;
    //  gsl_vector step1vect(ntax,0);
    //  for (int i=0; i<ntax; i++) {
    //      step1vect[i]=step1[i][0];
    //  }
    //cout<<"VCVinverse*J"<<endl<<step1<<endl;
    //   gsl_vector step2vect(1,0);
    //   step2vect=MatrixTimesVector(Jprime,step1vect);
    //   gsl_matrix step3(1,1,0);
    //   gsl_matrix step2(1,1,0);
    //   step2[0][0]=step2vect[0];
    //cout<<step2;
    //   matrixsingular=TestSingularity(step2);
    //   if (matrixsingular) {
    //       errormsg="Singular matrix (step2) during get ancestral state";
    //      throw XNexus(errormsg);
    //  }
    //  step3=Inverse(step2);
    //   gsl_matrix tipsasmatrix(ntax,1,0);
    //   gsl_matrix tipsasmatrixprime(1,ntax,0);
    //   for (int i=0; i<ntax; i++) {
    //       tipsasmatrix[i][0]=tips[i];
    //       tipsasmatrixprime[0][i]=tips[i];
    //   }
    //   gsl_matrix step4(ntax,1,0);
    //   step4=VCVinverse*tipsasmatrix;
    //   gsl_matrix step5(1,1,0);
    //   for (int i=0; i<ntax; i++) {
    //       step5[0][0]+=step4[i][0];
    //   }
    //   gsl_matrix finalvalue(1,1,0);
    //   finalvalue=step3*step5;
    //   ancestralstate=finalvalue[0][0];

	gsl_matrix_free(VCVinverse);
	gsl_matrix_free(VCVinversestart);
	gsl_matrix_free(tipsasmatrix);
	gsl_matrix_free(stepA);
	gsl_permutation_free(p);
    return ancestralstate;	
}


gsl_vector * BROWNIE::GetTipResiduals(gsl_vector * tips, double ancestralstate)
{
    int ntax=tips->size;
    gsl_vector *tipresiduals;
    tipresiduals=gsl_vector_calloc(ntax);
    //gsl_vector tipresiduals(ntax,0);
    //cout<<"GetTipResiduals ntax="<<ntax<<endl;
    //cout<<"Tips / tip residuals"<<endl;
    for (int i=0; i<ntax; i++) {
        // cout<<tips[i];
        //cout<<" / ";
        gsl_vector_set(tipresiduals,i,((gsl_vector_get(tips,i))-ancestralstate));
        //cout<<tipresiduals[i]<<endl;
    }
    return tipresiduals;
}

//Estimate rate
//Make sure to use tip residuals (tips minus estimated ancestral state)
double BROWNIE::EstimateRate(gsl_matrix * VCV, gsl_vector * tipresiduals)
{
    //rate=tipresiduals'*(inv(currenttreematrix))*tipresiduals/ntax
    double rateparameter;
    int ntax=VCV->size1;
    //matrixsingular=TestSingularity(VCV);
    //if (matrixsingular) {
    //    errormsg="Singular matrix (input VCV) during estimate rate";
    //   throw XNexus(errormsg);
    //}
    if (debugmode) {
        message="\nNow calculating inverse VCV for estimating rate.\n";
        PrintMessage();
    }
gsl_permutation * p = gsl_permutation_alloc (ntax);
int signum;
gsl_matrix * VCVinverse=gsl_matrix_calloc(ntax,ntax);
gsl_matrix * VCVinversestart=gsl_matrix_calloc(ntax,ntax);
gsl_matrix_memcpy (VCVinversestart, VCV);
if(debugmode) {
    message="debugging\n";
    for (int currentrow=0;currentrow<VCV->size1;currentrow++) {
        for (int currentcol=0;currentcol<VCV->size2;currentcol++) {
            message+=gsl_matrix_get(VCV,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }
    message+="\n\nduplicatedmatrix\n";
    for (int currentrow=0;currentrow<VCVinversestart->size1;currentrow++) {
        for (int currentcol=0;currentcol<VCVinversestart->size2;currentcol++) {
            message+=gsl_matrix_get(VCVinversestart,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }

    PrintMessage();
}
gsl_linalg_LU_decomp (VCVinversestart,p, &signum);
if(debugmode) {
    message="\n\nLUmatrix\n";
    for (int currentrow=0;currentrow<VCVinversestart->size1;currentrow++) {
        for (int currentcol=0;currentcol<VCVinversestart->size2;currentcol++) {
            message+=gsl_matrix_get(VCVinversestart,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }

    PrintMessage();
}

gsl_linalg_LU_invert (VCVinversestart,p, VCVinverse);

if(debugmode) {
    message="\n\nInverseMatrix\n";
    for (int currentrow=0;currentrow<VCVinverse->size1;currentrow++) {
        for (int currentcol=0;currentcol<VCVinverse->size2;currentcol++) {
            message+=gsl_matrix_get(VCVinverse,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }

    PrintMessage();
}


//VCVinverse=Inverse(VCV);
if (debugmode) {
    message="\nFinished calculating inverse VCV for estimating rate.\n";
    PrintMessage();
}
gsl_vector *step1vect;
step1vect=gsl_vector_calloc(ntax);
//gsl_vector step1vect(ntax,0);
gsl_blas_dgemv (CblasNoTrans,1, VCVinverse, tipresiduals,0, step1vect); ///TEST THIS
                                                                        //step1vect=MatrixTimesVector(VCVinverse,tipresiduals);
gsl_vector *step2vect;
step2vect=gsl_vector_calloc(1);
//gsl_vector step2vect(1,0);
gsl_matrix *tipsasmatrixprime=gsl_matrix_calloc(1,ntax);
for (int i=0; i<ntax; i++) {
    gsl_matrix_set(tipsasmatrixprime,0,i,(gsl_vector_get(tipresiduals,i)));
}
gsl_blas_dgemv (CblasNoTrans,1, tipsasmatrixprime, step1vect,0, step2vect); ///TEST THIS
                                                                            //step2vect=MatrixTimesVector(tipsasmatrixprime,step1vect);
rateparameter=gsl_vector_get(step2vect,0)/ntax;
if (debugmode) {
    message="rate estimate is ";
    message+=rateparameter;
    PrintMessage();
}
gsl_matrix_free(VCVinverse);
gsl_matrix_free(VCVinversestart);
gsl_vector_free(step1vect);
gsl_vector_free(step2vect);
gsl_matrix_free(tipsasmatrixprime);
gsl_permutation_free (p);
return rateparameter;
}





////////////////////////////////////////////////
double BROWNIE::GetLikelihoodUnderLindy2_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((BROWNIE*)obj)->GetLikelihoodUnderLindy2(variables);
	return temp;
}

double BROWNIE::GetLikelihoodUnderLindy2(const gsl_vector * variables)
{
	double rateA=gsl_vector_get(variables,0);
	double rateB=gsl_vector_get(variables,1);
	double likelihood=(CalculateDiscreteLindy2(rateA,rateB));
	return likelihood;
}

double BROWNIE::GetLikelihoodUnderLindy1_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((BROWNIE*)obj)->GetLikelihoodUnderLindy1(variables);
	return temp;
}

double BROWNIE::GetLikelihoodUnderLindy1(const gsl_vector * variables)
{
	double rateA=gsl_vector_get(variables,0);
	double likelihood=(CalculateDiscreteLindy1(rateA));
	return likelihood;
}

gsl_vector * BROWNIE::LindyGeneralOptimization(int ChosenModel)
{
	//Model 1=One param
	//Model 2=Two param
	double bestlikelihood=GSL_POSINF;
	int hitlimitscount=0;
	size_t np;
	if (ChosenModel==1) {
		np = 1;
	}
	else if (ChosenModel==2) {
		np = 2;
	}
	gsl_vector * results=gsl_vector_calloc(np);
	double estimates[randomstarts][np];
	double startingvalues[randomstarts][np];
	double likelihoods[randomstarts][1];
	for (int startnum=0;startnum<randomstarts;startnum++) {
		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
		gsl_multimin_fminimizer *s = NULL;
		gsl_vector *ss, *x;
		size_t iter = 0, i;
		int status;
		double size;
		bool hitlimits=false;
		/* Initial vertex size vector */
		ss = gsl_vector_alloc (np);
		gsl_vector_set_all (ss, stepsize);
		
		/* Starting point */
		x = gsl_vector_calloc (np);
		if (ChosenModel==1) {
			gsl_vector_set (x,0,gsl_ran_exponential (r,1.0)); //starting rate
			startingvalues[startnum][0]=gsl_vector_get(x,0);
		}
		else if (ChosenModel==2) {
			gsl_vector_set(x,0,gsl_ran_exponential (r,1.0)); //starting rate
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			gsl_vector_set(x,1,gsl_ran_exponential (r,1.0)); 
			startingvalues[startnum][1]=gsl_vector_get(x,1);
		}
		BROWNIE *pt;
		pt=(this);
		double (*F)(const gsl_vector *, void *);
		if (ChosenModel==1) {
			F = &BROWNIE::GetLikelihoodUnderLindy1_gsl;
		}
		else if (ChosenModel==2) {
			F = &BROWNIE::GetLikelihoodUnderLindy2_gsl;
		}
		gsl_multimin_function minex_func;
		minex_func.f=*F;
		minex_func.params=pt;
		minex_func.n = np;
		s = gsl_multimin_fminimizer_alloc (T, np);
		gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
		do
		{
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);
			if (status!=0) { //0 Means it's a success in c++, but not in C
				printf ("error: %s\n", gsl_strerror (status));
				break;
			}
			size = gsl_multimin_fminimizer_size (s);
			status = gsl_multimin_test_size (size, stoppingprecision); //since we want more precision
			if (status == GSL_SUCCESS)
			{
				//printf ("converged to minimum at\n");
			}
		}
		while (status == GSL_CONTINUE && iter < maxiterations);
		if (s->fval<bestlikelihood) {
			gsl_vector_memcpy(results,s->x);
			bestlikelihood=s->fval;
		}
		for (int parameternumber=0; parameternumber<np; parameternumber++) {
			estimates[startnum][parameternumber]=gsl_vector_get(s->x,parameternumber);
		}
		if (iter==maxiterations) {
			hitlimits=true;
			hitlimitscount++;
		}
		message="Replicate ";
		message+=startnum+1;
		if (hitlimits) {
			message+=" **WARNING**";
		}
		message+="\n   NM iterations needed = ";
		int iterationsrequired=iter;
		message+=iterationsrequired;
		if (hitlimits) {
			message+=" **Max iterations hit; see WARNING below**";
		}
		message+="\n   LnL = ";
		char outputstring[60];
		sprintf(outputstring,"%60.45f",-1.0*(s->fval));
		message+=outputstring;
		// message+="\n   Rate = ";
		// message+=gsl_vector_get(s->x,0);
		if (detailedoutput) {
			PrintMessage();
		}
		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);
	//	if(detailedoutput==false) {
	//		ProgressBar(0);
	//	}
	}
	if (hitlimitscount>0) {
		message="\n----------------------------------------------------------------------------\n WARNING: Out of ";
		message+=randomstarts;
		message+=" optimization starts, ";
		message+=hitlimitscount;
		if (hitlimitscount==1) {
			message+=" was ";
		}
		else {
			message+=" were ";
		}
		message+="stopped by hitting\n  the maximum # of iterations. This means that those replicates\n  may not even have hit the local maximum.\n\n  You can increase the maximum number of iterations or decrease the\n  precision with the NumOpt command. You could also consider\n  increasing the number of random starts using that same command.\n\n  If this happened on a small proportion of replicates, though,\n  or if the precision (below) is good enough, don't worry about it.\n----------------------------------------------------------------------------";
		if (detailedoutput) {
			PrintMessage();
		}
		else if ((randomstarts-hitlimitscount)<10 && (hitlimitscount/randomstarts)>.1) {
			PrintMessage();
		}
	}
	gsl_vector * finalvector=gsl_vector_calloc((2*np)+1);
	for (int position=0; position<np; position++) {
		gsl_vector_set(finalvector,position,gsl_vector_get(results,position));
		double paramestimate[randomstarts];
		for (int startnumber=0;startnumber<randomstarts;startnumber++) {
			paramestimate[startnumber]=estimates[startnumber][position];
		}
		gsl_vector_set(finalvector,position+np,gsl_stats_sd(paramestimate,1,randomstarts));
	}
	gsl_vector_set(finalvector,(2*np),bestlikelihood);
	return finalvector;
}





////////////////////////////////////////////////


//Calculates the likelihood of discrete character chosenchar on tree chosentree
double BROWNIE::CalculateDiscreteCharLnL(gsl_matrix * RateMatrix, gsl_vector * ancestralstatevector, bool returnancstates)
{
	double lnL=0;
	map<Node*, vector<double> > stateprobatnodes;
	Tree T=intrees.GetIthTree(chosentree-1);
	NodeIterator <Node> n (T.GetRoot());
	NodePtr currentnode = n.begin();
	for (int i=0;i<ancestralstatevector->size;i++) {
		(stateprobatnodes[T.GetRoot()]).push_back(gsl_vector_get(ancestralstatevector,i));
	}
	while (currentnode)
	{
		if (currentnode!=T.GetRoot() && !currentnode->IsLeaf()) {
			cout<<"RateMatrix is\n"<<gsl_matrix_get(RateMatrix,0,0)<<"\t"<<gsl_matrix_get(RateMatrix,0,1)<<"\n"<<gsl_matrix_get(RateMatrix,1,0)<<"\t"<<gsl_matrix_get(RateMatrix,1,1)<<endl;
			gsl_matrix * Pmatrix=ComputeTransitionProb(RateMatrix,currentnode->GetEdgeLength());
			for(int j=0;j<ancestralstatevector->size;j++) {
				double probofstatej=0;
				for (int i=0;i<ancestralstatevector->size;i++) {
					probofstatej+=((stateprobatnodes[currentnode->GetAnc()])[i])*gsl_matrix_get(Pmatrix,i,j);
				}
				(stateprobatnodes[currentnode]).push_back(probofstatej);
			}
		}
		else if (currentnode!=T.GetRoot() && currentnode->IsLeaf()) {
			int statenumber=characters->GetInternalRepresentation(taxa->FindTaxon(currentnode->GetLabel()),chosenchar);
			cout<<"RateMatrix is\n"<<gsl_matrix_get(RateMatrix,0,0)<<"\t"<<gsl_matrix_get(RateMatrix,0,1)<<"\n"<<gsl_matrix_get(RateMatrix,1,0)<<"\t"<<gsl_matrix_get(RateMatrix,1,1)<<endl;
			gsl_matrix * Pmatrix=ComputeTransitionProb(RateMatrix,currentnode->GetEdgeLength());
			for(int j=0;j<ancestralstatevector->size;j++) {
				double probofstatej=0;
				if (j==statenumber) {
					for (int i=0;i<ancestralstatevector->size;i++) {
						probofstatej+=((stateprobatnodes[currentnode->GetAnc()])[i])*gsl_matrix_get(Pmatrix,i,j);
					}
				} 
				(stateprobatnodes[currentnode]).push_back(probofstatej);
				lnL+=log(probofstatej);
			}
		}
		currentnode = n.next();
	}
	return lnL;
}

double BROWNIE::CalculateDiscreteLindy1(double rateA)
{
	double neglnL=CalculateDiscreteLindy2(rateA, rateA);
	return neglnL;
}

//Calculates the likelihood of going from state 1 to state 0 on one tree across all chars
double BROWNIE::CalculateDiscreteLindy2(double rateA, double rateB)
{
//	gsl_matrix * testmat = gsl_matrix_calloc(2,2);
//	gsl_matrix_set(testmat,0,0,-0.5);
//	gsl_matrix_set(testmat,0,1,0.5);
//	gsl_matrix_set(testmat,1,0,0.7);
//	gsl_matrix_set(testmat,1,1,-0.7);
//	gsl_matrix * testmatout=ComputeTransitionProb(testmat,1.0);
//	cout<<"testmatnow with brlen 1\n"<<gsl_matrix_get(testmatout,0,0)<<"\t"<<gsl_matrix_get(testmatout,0,1)<<"\n"<<gsl_matrix_get(testmatout,1,0)<<"\t"<<gsl_matrix_get(testmatout,1,1)<<endl;
//	testmatout=ComputeTransitionProb(testmat,0.0);
//	cout<<"testmatnow with brlen 0\n"<<gsl_matrix_get(testmatout,0,0)<<"\t"<<gsl_matrix_get(testmatout,0,1)<<"\n"<<gsl_matrix_get(testmatout,1,0)<<"\t"<<gsl_matrix_get(testmatout,1,1)<<endl;

	
	double neglnL=0.0;
	if (rateA<0.0 || rateB<0.0) {
		neglnL=GSL_POSINF;
	}
	else {
//		cout<<"rateA is "<<rateA<<endl;
		gsl_vector* ancestralstatevector=gsl_vector_calloc(2);
		gsl_vector_set(ancestralstatevector,0,0.0); //state 1 at root
		gsl_vector_set(ancestralstatevector,1,1.0);
		gsl_matrix* ratematrixA=gsl_matrix_calloc(2,2);
		//cout<<"Qmatrix starts as\n"<<gsl_matrix_get(ratematrixA,0,0)<<"\t"<<gsl_matrix_get(ratematrixA,0,1)<<"\n"<<gsl_matrix_get(ratematrixA,1,0)<<"\t"<<gsl_matrix_get(ratematrixA,1,1)<<endl;
		gsl_matrix_set(ratematrixA,1,0,rateA);
		//cout<<"Qmatrix now is is\n"<<gsl_matrix_get(ratematrixA,0,0)<<"\t"<<gsl_matrix_get(ratematrixA,0,1)<<"\n"<<gsl_matrix_get(ratematrixA,1,0)<<"\t"<<gsl_matrix_get(ratematrixA,1,1)<<endl;
		gsl_matrix_set(ratematrixA,1,1,0.0-rateA);
//		cout<<"ratematrixA now is is\n"<<gsl_matrix_get(ratematrixA,0,0)<<"\t"<<gsl_matrix_get(ratematrixA,0,1)<<"\n"<<gsl_matrix_get(ratematrixA,1,0)<<"\t"<<gsl_matrix_get(ratematrixA,1,1)<<endl;
		gsl_matrix* ratematrixB=gsl_matrix_calloc(2,2);
		gsl_matrix_set(ratematrixB,1,0,rateB);
		gsl_matrix_set(ratematrixB,1,1,0.0-rateB);
//		cout<<"ratematrixB now is is\n"<<gsl_matrix_get(ratematrixB,0,0)<<"\t"<<gsl_matrix_get(ratematrixB,0,1)<<"\n"<<gsl_matrix_get(ratematrixB,1,0)<<"\t"<<gsl_matrix_get(ratematrixB,1,1)<<endl;
		int nchartotal=characters->GetNChar();
		for (int currentchar=1;currentchar<=nchartotal;currentchar++) {
			//cout<<"rateA is "<<rateA<<endl;
			map<Node*, vector<double> > stateprobatnodes;
			Tree T=intrees.GetIthTree(chosentree-1);
			PreorderIterator <Node> n (T.GetRoot());
			NodePtr currentnode = n.begin();
			vector<double> ancstatevec;
			for (int i=0;i<ancestralstatevector->size;i++) {
				ancstatevec.push_back(gsl_vector_get(ancestralstatevector,i));
//				cout<<"anc state "<<gsl_vector_get(ancestralstatevector,i)<<endl;
			}
			stateprobatnodes[T.GetRoot()]=ancstatevec;
			while (currentnode)
			{
				if (currentnode!=T.GetRoot() && !currentnode->IsLeaf()) {
					gsl_matrix * Pmatrix=gsl_matrix_calloc(2,2);
					if ((currentnode->GetStateOrder())[0]==0) {
//						cout<<"\n\nQmatrix is\n"<<gsl_matrix_get(ratematrixA,0,0)<<"\t"<<gsl_matrix_get(ratematrixA,0,1)<<"\n"<<gsl_matrix_get(ratematrixA,1,0)<<"\t"<<gsl_matrix_get(ratematrixA,1,1)<<endl;
						Pmatrix=ComputeTransitionProb(ratematrixA,currentnode->GetEdgeLength());
//						cout<<"Qmatrix is\n"<<gsl_matrix_get(ratematrixA,0,0)<<"\t"<<gsl_matrix_get(ratematrixA,0,1)<<"\n"<<gsl_matrix_get(ratematrixA,1,0)<<"\t"<<gsl_matrix_get(ratematrixA,1,1)<<endl;

					}
					else {
						Pmatrix=ComputeTransitionProb(ratematrixB,currentnode->GetEdgeLength());
					}
//					cout<<"edgelength was "<<currentnode->GetEdgeLength()<<endl;
//					cout<<"Pmatrix is\n"<<gsl_matrix_get(Pmatrix,0,0)<<"\t"<<gsl_matrix_get(Pmatrix,0,1)<<"\n"<<gsl_matrix_get(Pmatrix,1,0)<<"\t"<<gsl_matrix_get(Pmatrix,1,1)<<endl;
					vector<double> statesatthisnode;
					for(int j=0;j<ancestralstatevector->size;j++) {
						double probofstatej=0;
						vector<double> statesatpreviousnode=stateprobatnodes[currentnode->GetAnc()];
						for (int i=0;i<ancestralstatevector->size;i++) {
//							cout<<"for j="<<j<<" prob of state "<<i<<" at prev node = "<<statesatpreviousnode[i]<<" and Pij="<<gsl_matrix_get(Pmatrix,i,j)<<endl;
							probofstatej+=(statesatpreviousnode[i])*gsl_matrix_get(Pmatrix,i,j);
						}
						statesatthisnode.push_back(probofstatej);
//						cout<<"j="<<j<<" prob of state j="<<probofstatej<<endl;
					}
					stateprobatnodes[currentnode]=statesatthisnode;
				}
				else if (currentnode!=T.GetRoot() && currentnode->IsLeaf()) {
					int statenumber=characters->GetInternalRepresentation(taxa->FindTaxon(currentnode->GetLabel()),currentchar-1);
					gsl_matrix * Pmatrix=gsl_matrix_calloc(2,2);
					if ((currentnode->GetStateOrder())[0]==0) {
						Pmatrix=ComputeTransitionProb(ratematrixA,currentnode->GetEdgeLength());
					}
					else {
						Pmatrix=ComputeTransitionProb(ratematrixB,currentnode->GetEdgeLength());
					}
					for(int j=0;j<ancestralstatevector->size;j++) {
						double probofstatej=0;
						if (j==statenumber) {
							vector<double> statesatpreviousnode=stateprobatnodes[currentnode->GetAnc()];
							for (int i=0;i<ancestralstatevector->size;i++) {
								probofstatej+=(statesatpreviousnode[i])*gsl_matrix_get(Pmatrix,i,j);
							}
							neglnL+=(-1.0)*log(probofstatej);
							//cout<<"For tip with state "<<statenumber<<" probofstate="<<probofstatej<<" prob at previousnode="<<statesatpreviousnode[statenumber]<<" -lnL total so far = "<<neglnL<<endl;
						} 
					}
				}
				currentnode = n.next();
			}
			//PreorderIterator <Node> m (T.GetRoot());
			//NodePtr currentnode2 = m.begin();
			//message="\nAncestralStates\n0\t1\n";
			//PrintMessage();
			//while (currentnode2) {
			//	message="";
			//	message+=(stateprobatnodes[currentnode])[0];
			//	message+="\t";
			//	message+=(stateprobatnodes[currentnode])[1];
			//	PrintMessage();
			//	currentnode2=m.next();
			//}
			

		}
	}
	return neglnL;
}

//uses pagel 1994's formula, P(t)=exp(Qt)+c=C*exp(Dt)*C^-1, where C is eigenvectors of Q and D is eigenvalues (in diagonal matrix)
gsl_matrix * BROWNIE::ComputeTransitionProb(gsl_matrix *RateMatrix, double brlen) {
	int dimension=RateMatrix->size1;
	gsl_matrix *ratematrixTMP=gsl_matrix_calloc(dimension,dimension);
	gsl_matrix_memcpy (ratematrixTMP, RateMatrix);
//	cout<<"ratematrixTMP is\n"<<gsl_matrix_get(ratematrixTMP,0,0)<<"\t"<<gsl_matrix_get(ratematrixTMP,0,1)<<"\n"<<gsl_matrix_get(ratematrixTMP,1,0)<<"\t"<<gsl_matrix_get(ratematrixTMP,1,1)<<endl;
	gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (dimension);
	gsl_vector_complex *eigenvalues=gsl_vector_complex_calloc(dimension);
	gsl_matrix_complex *eigenvectors=gsl_matrix_complex_calloc(dimension,dimension); //Matrix C as in Pagel 1994
	gsl_matrix *transitionmatrix=gsl_matrix_calloc(dimension,dimension);
	gsl_matrix_complex *transitionmatrixcomplex=gsl_matrix_complex_calloc(dimension,dimension);

	gsl_eigen_nonsymmv (ratematrixTMP, eigenvalues, eigenvectors, w);
//	cout<<"eigenvalues are "<<GSL_REAL(gsl_vector_complex_get(eigenvalues,0))<<" "<<GSL_REAL(gsl_vector_complex_get(eigenvalues,1))<<endl;
	gsl_matrix_complex *eigenvaluematrix=gsl_matrix_complex_calloc(dimension,dimension);
	//gsl_matrix *eigenvectorsreal=gsl_matrix_calloc(dimension,dimension);
	for (int i=0;i<dimension;i++) {
		gsl_matrix_complex_set(eigenvaluematrix,i,i,gsl_complex_exp(gsl_complex_mul_real(gsl_vector_complex_get(eigenvalues,i),brlen)));
	//	for (int j=0;j<dimension;j++) {
		//	gsl_matrix_set(eigenvectorsreal,i,j,GSL_REAL(gsl_matrix_complex_get(eigenvectors,i,j)));
	//	}
	}
	
	
	gsl_matrix_complex *inverseeigenvectorsstart=gsl_matrix_complex_calloc(dimension,dimension);
	gsl_matrix_complex *inverseeigenvectors=gsl_matrix_complex_calloc(dimension,dimension);
//	cout<<"eigenvaluematrix (after taking exp(t*lambda): \n"<<GSL_REAL(gsl_matrix_complex_get(eigenvaluematrix,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(eigenvaluematrix,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(eigenvaluematrix,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(eigenvaluematrix,1,1))<<endl;

//	cout<<"eigenvectors: \n"<<GSL_REAL(gsl_matrix_complex_get(eigenvectors,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(eigenvectors,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(eigenvectors,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(eigenvectors,1,1))<<endl;
	gsl_permutation * p = gsl_permutation_alloc (dimension);
	int signum;
	gsl_matrix_complex_memcpy (inverseeigenvectorsstart, eigenvectors);
	gsl_linalg_complex_LU_decomp (inverseeigenvectorsstart,p, &signum);
//	cout<<"inverseeigenvectorsstart: \n"<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectorsstart,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectorsstart,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectorsstart,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectorsstart,1,1))<<endl;
	gsl_linalg_complex_LU_invert (inverseeigenvectorsstart,p, inverseeigenvectors);
//	cout<<"inverseeigenvectors: \n"<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectors,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectors,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectors,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectors,1,1))<<endl;

	//cout<<"inverseeigenvectors: "<<gsl_matrix_get(inverseeigenvectors,0,0)<<" "<<gsl_matrix_get(inverseeigenvectors,0,1)<<endl<<gsl_matrix_get(inverseeigenvectors,1,0)<<"\t"<<gsl_matrix_get(inverseeigenvectors,1,1)<<endl;

	gsl_matrix_complex *stepA = gsl_matrix_complex_calloc(dimension,dimension);
	//gsl_complex one=gsl_complex_rect (1.0, 1.0);
	gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, eigenvaluematrix, inverseeigenvectors,GSL_COMPLEX_ZERO, stepA);
	//gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, one, eigenvaluematrix, inverseeigenvectors,one, stepA);
//	cout<<"stepA: \n"<<GSL_REAL(gsl_matrix_complex_get(stepA,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(stepA,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(stepA,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(stepA,1,1))<<endl;
	gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, eigenvectors, stepA,GSL_COMPLEX_ZERO, transitionmatrixcomplex);
	gsl_eigen_nonsymmv_free(w);
	gsl_vector_complex_free(eigenvalues);
	gsl_matrix_complex_free(eigenvectors);
	gsl_matrix_complex_free(inverseeigenvectorsstart);
	gsl_matrix_complex_free(inverseeigenvectors);
	gsl_matrix_complex_free(stepA);
	gsl_matrix_free(ratematrixTMP);
	gsl_permutation_free(p);
//	gsl_matrix *transitionmatrix=gsl_matrix_calloc(dimension,dimension);
	for (int i=0;i<dimension;i++) {
		for (int j=0;j<dimension;j++) {
			gsl_matrix_set(transitionmatrix,i,j,GSL_REAL(gsl_matrix_complex_get(transitionmatrixcomplex,i,j)));
		}
	}
	gsl_matrix_complex_free(transitionmatrixcomplex);
//	cout<<"transitionmatrixcomplex: \n"<<GSL_REAL(gsl_matrix_complex_get(transitionmatrixcomplex,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(transitionmatrixcomplex,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(transitionmatrixcomplex,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(transitionmatrixcomplex,1,1))<<endl;
//	cout<<"transitionmatrix: \n"<<gsl_matrix_get(transitionmatrix,0,0)<<" "<<gsl_matrix_get(transitionmatrix,0,1)<<endl<<gsl_matrix_get(transitionmatrix,1,0)<<"\t"<<gsl_matrix_get(transitionmatrix,1,1)<<endl;
	return transitionmatrix;
	
}


/** @method HandleTimeSlice
*
* Allows user to specify how to slice trees so that different time periods have different rates.
*
*/
void BROWNIE::HandleTimeSlice( NexusToken& token )
{
    bool noerror=true;
    bool adequateinput=false;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (adequateinput==false) {
                message="Insufficient input: type \"timeslice ?\" for help";
                PrintMessage();
            }
            break;
        }
        else if( token.Equals("?") ) {
            adequateinput=true;
            message="Usage: Timeslice splits=(<times>) models=(<numbers>) [options...]\n\n";
            message+="Assigns models to intervals.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\nSplits       (<numbers>)                          *None";
            message+="\nModels       (<numbers>)                          *None";
            // message+="\nTaxset       <taxset name>                        *All";
            // message+="\nFromRoot     No|Yes                               *No";
            message+="\n                                                 *Option is nonpersistent\n\n";
            message+="This will allow assignment of different models on the specified intervals across all trees.\nFor example, to specify one rate parameter for the interval from the present to 50 MYA,\nanother rate parameter from 50 MYA to 65 MYA, and the first rate parameter again from 65 MYA\nto the root of the phylogeny:\n\n  Time Model\n   0     1\n   5     1\n   .     .\n  45     1\n__50_____1_\n| 50     2 |\n| 55     2 |\n| 60     2 |\n|_65_____2_|\n  65     1\n  70     1\n   .     .\n\none would type\n\nTimeslice splits=(50 65) models=(1 2 1);\n\nSplits specifies where one model changes to the other;\nthere should be one fewer split than there are intervals.\nEdges spanning a split are properly divided\n(so a given path may have more than one model).\n\nBy default, times are measured from the most recent taxon. \nTo specify that they should be measured from the root instead, specify fromroot=yes;\nThis becomes most important when the root node may be of different depths in different trees.\nNote that you are currently limited to only 9 splits (if this becomes a problem, email bcomeara@ucdavis.edu).";
            PrintMessage();
        }
        else if( token.Abbreviation("Splits") ) {
            token.GetNextToken();
            if (!token.Equals("(")) {
                errormsg="Expecting next token to be a \'(\' but instead got ";
                errormsg+=token.GetToken();
                throw XNexus( errormsg);
            }
            token.GetNextToken(); //Should be a number
            int splitposcount=0;
            while (!token.Equals(")")) {
                nxsstring numbernexus;
                numbernexus=GetNumber(token);
                timeslicetimes[0]=atof( numbernexus.c_str() );
                splitposcount++;
                token.GetNextToken();
            }
        }
        else if( token.Abbreviation("Models") ) {
            token.GetNextToken();
            if (!token.Equals("(")) {
                errormsg="Expecting next token to be a \'(\' but instead got ";
                errormsg+=token.GetToken();
                throw XNexus( errormsg);
            }
            token.GetNextToken(); //Should be a number
            int modelposcount=0;
            while (!token.Equals(")")) {
                nxsstring numbernexus;
                numbernexus=GetNumber(token);
                timeslicemodels[0]=atoi( numbernexus.c_str() )-1;
                modelposcount++;
                token.GetNextToken();
            }
        }
        else if (token.Abbreviation("Taxset") ) {
        }
        else if (token.Abbreviation("FromRoot")  ) {
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading RateTest command. Type \"RateTest ?\" for help.";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
        //timeslicetimes[9]=GSL_POSINF;
        //timeslicemodels;

    }
    //if (__________EVERYTHING IS OKAY__________) {
    //Do the set model thing
    //}

}



int main(int argc, char *argv[])
{
    //Code on random number generation ripped from GSL documentation



    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc (T);
    BROWNIE brownie;
    bool inputfilegiven=false;
    if (argc>1) {
        for (int i = 1; i < argc; i++) {
            inputfilegiven=true;
            nxsstring fn=argv[i];
            brownie.RunCmdLine(inputfilegiven, fn );
        }
    }
    else {
        nxsstring fn="a";
        brownie.RunCmdLine(inputfilegiven, fn);
    }
    return 0;
}

