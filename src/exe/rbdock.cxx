/***********************************************************************
* The rDock program was developed from 1998 - 2006 by the software team 
* at RiboTargets (subsequently Vernalis (R&D) Ltd).
* In 2006, the software was licensed to the University of York for 
* maintenance and distribution.
* In 2012, Vernalis and the University of York agreed to release the 
* program as Open Source software.
* This version is licensed under GNU-LGPL version 3.0 with support from
* the University of Barcelona.
* http://rdock.sourceforge.net/
***********************************************************************/



//Main docking application
#include <iomanip>
using std::setw;

#ifdef _VISUAL_STUDIO
#else
#include <popt.h>		// for command-line parsing
#endif

#include <errno.h>

#if 1 // Xability (2017)
#include <mpi.h>
#endif // Xability (2017)

#include "RbtBiMolWorkSpace.h"
#include "RbtMdlFileSource.h"
#include "RbtMdlFileSink.h"
#include "RbtCrdFileSink.h"
#include "RbtParameterFileSource.h"
#include "RbtPRMFactory.h"
#include "RbtSFFactory.h"
#include "RbtTransformFactory.h"
#include "RbtRand.h"
#include "RbtModelError.h"
#include "RbtDockingError.h"
#include "RbtLigandError.h"
#include "RbtFilter.h"
#include "RbtSFRequest.h"
#include "RbtFileError.h"

#if 1 // AdvanceSoft (2018)
#include "RbtBaseTransform.h"
#include "RbtChromFactory.h"
#include <iostream>
#include <time.h>
#include "RbtGenome.h"
#endif // AdvanceSoft(2018)

const RbtString EXEVERSION = " ($Id: //depot/dev/client3/rdock/2013.1/src/exe/rbdock.cxx#4 $)";
//Section name in docking prm file containing scoring function definition
const RbtString _ROOT_SF = "SCORE";
const RbtString _RESTRAINT_SF = "RESTR";
const RbtString _ROOT_TRANSFORM = "DOCK";

void PrintUsage(void)
{
  cout << endl << "Usage:" << endl;
  cout << "rbdock -i <sdFile> -o <outputRoot> -r <recepPrmFile> -p <protoPrmFile> [-n <nRuns>] [-ap] [-an] [-allH]" << endl;
  cout << "       [-t <targetScore|targetFilterFile>] [-c] [-T <traceLevel>] [-s <rndSeed>]" << endl;
  cout << endl << "Options:\t-i <sdFile> - input ligand SD file" << endl;
  cout << "\t\t-o <outputRoot> - root name for output file(s)" << endl;
  cout << "\t\t-r <recepPrmFile> - receptor parameter file " << endl;
  cout << "\t\t-p <protoPrmFile> - docking protocol parameter file" << endl;
  cout << "\t\t-n <nRuns> - number of runs/ligand (default=1)" << endl;
  cout << "\t\t-ap - protonate all neutral amines, guanidines, imidazoles (default=disabled)" << endl;
  cout << "\t\t-an - deprotonate all carboxylic, sulphur and phosphorous acid groups (default=disabled)" << endl;
  cout << "\t\t-allH - read all hydrogens present (default=polar hydrogens only)" << endl;
  cout << "\t\t-t - score threshold OR filter file name" << endl;
  cout << "\t\t-c - continue if score threshold is met (use with -t <targetScore>, default=terminate ligand)" << endl;
  cout << "\t\t-T <traceLevel> - controls output level for debugging (0 = minimal, >0 = more verbose)" << endl;
  cout << "\t\t-s <rndSeed> - random number seed (default=from sys clock)" << endl;
#if 1 // AdvanceSoft (2018)
  cout << "\t\t-sideChain - put side chains of a receptor into chromosomes (default=disabled)" << endl;
  cout << "\t\t-l <iIsland> - number of islands for island-model GA (default=1)" << endl;
  cout << "\t\t-m - migration for island-model GA (default=disabled)" << endl;
  cout << "\t\t-g <mmgbpbsa_starting_point> - MM-GB/PBSA socre will be used for docking (default=disabled)" << endl;
  cout << "\t\t-M - parallel computation of MM-GB/PBSA rDock (default=disabled)" << endl;
  cout << "\t\t-outputRec - output a PDB file including the receptor structure(s) after docking" << endl;
#endif // AdvanceSoft (2018)
}

#if 1 // AdvanceSoft (2018)
void ConcatenateSDFiles(const RbtString& strSDFile, const RbtString& strOutput)
{
	RbtInt iMPIRank, nMPISize;
	MPI_Comm_rank(MPI_COMM_WORLD, &iMPIRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);

	RbtString strSD;

	int num_root_proc = m_iSdf;

	ifstream ifs(strSDFile.c_str());
	if (ifs) {
		strSD = RbtString(std::istreambuf_iterator<char>(ifs),
			std::istreambuf_iterator<char>());
	}
	ifs.close();

	vector<RbtInt> nSizeAll(num_root_proc, 0);
	RbtInt nSize = strSD.length();

	int iMPIRank_Ad = 0;

	if (iMPIRank != 0) {
		MPI_Send(&nSize, 1, MPI_INT, 0, iMPIRank, MPI_COMM_WORLD);
	}
	else {

		nSizeAll[0] = nSize;

		for (int i = 1; i < num_root_proc; ++i) {
			iMPIRank_Ad += m_iIsland;
			MPI_Status status;
			MPI_Recv(&nSizeAll[i], 1, MPI_INT, iMPIRank + iMPIRank_Ad, iMPIRank + iMPIRank_Ad,
				MPI_COMM_WORLD, &status);
		}
	}

	vector<RbtString> vstrSD(num_root_proc);
	iMPIRank_Ad = 0;
	for (int i = 1; i < num_root_proc; ++i) {
		iMPIRank_Ad += m_iIsland;
		if (iMPIRank == iMPIRank_Ad) {
			MPI_Send(&strSD[0], nSize, MPI_CHAR, 0, iMPIRank_Ad, MPI_COMM_WORLD);
		}
		if (iMPIRank == 0) {
			vstrSD[i].resize(nSizeAll[i], 0);
			MPI_Status status;
			MPI_Recv(&vstrSD[i][0], nSizeAll[i], MPI_CHAR, iMPIRank_Ad, iMPIRank_Ad,
				MPI_COMM_WORLD, &status);
		}
	}

	if (iMPIRank == 0) {
		for (int i = 1; i < num_root_proc; ++i) {
			strSD += vstrSD[i];
		}
		multimap<RbtInt, RbtString> recordmap;

		std::stringstream ss(strSD.c_str());
		while (true) {
			string strRec;
			RbtString line;
			RbtInt iRecord = 0;
			while (getline(ss, line, '\n')) {
				if (line == ">  <Record>") {
					getline(ss, line, '\n');
					iRecord = atoi(line.c_str());
					getline(ss, line, '\n');
					continue;
				}
				strRec += line + "\n";
				if (line == "$$$$") {
					break;
				}
			}
			if (0 < iRecord) {
				recordmap.insert(pair<RbtInt, RbtString>(iRecord, strRec));
			}
			else {
				break;
			}
		}
		strSD = "";
		multimap<RbtInt, RbtString>::iterator it = recordmap.begin();
		for (; it != recordmap.end(); ++it) {
			strSD += it->second;
		}

		ofstream ofs(strOutput.c_str());
		if (ofs) {
			ofs << strSD;
		}
		ofs.close();
	}
}
#endif // AdvanceSoft(2018)


#if 1 // AdvanceSoft (2018)
void ConcatenatePDBFiles(const RbtString& strPDBFile, const RbtString& strOutput)
{
	RbtInt iMPIRank, nMPISize;
	MPI_Comm_rank(MPI_COMM_WORLD, &iMPIRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);

	RbtString strSD;

	int num_root_proc = m_iSdf;

	ifstream ifs(strPDBFile.c_str());
	if (ifs) {
		strSD = RbtString(std::istreambuf_iterator<char>(ifs),
			std::istreambuf_iterator<char>());
	}
	ifs.close();

	vector<RbtInt> nSizeAll(num_root_proc, 0);
	RbtInt nSize = strSD.length();

	int iMPIRank_Ad = 0;

	if (iMPIRank != 0) {
		MPI_Send(&nSize, 1, MPI_INT, 0, iMPIRank, MPI_COMM_WORLD);
	}
	else {

		nSizeAll[0] = nSize;

		for (int i = 1; i < num_root_proc; ++i) {
			iMPIRank_Ad += m_iIsland;
			MPI_Status status;
			MPI_Recv(&nSizeAll[i], 1, MPI_INT, iMPIRank + iMPIRank_Ad, iMPIRank + iMPIRank_Ad,
				MPI_COMM_WORLD, &status);
		}
	}

	vector<RbtString> vstrSD(num_root_proc);
	iMPIRank_Ad = 0;
	for (int i = 1; i < num_root_proc; ++i) {
		iMPIRank_Ad += m_iIsland;
		if (iMPIRank == iMPIRank_Ad) {
			MPI_Send(&strSD[0], nSize, MPI_CHAR, 0, iMPIRank_Ad, MPI_COMM_WORLD);
		}
		if (iMPIRank == 0) {
			vstrSD[i].resize(nSizeAll[i], 0);
			MPI_Status status;
			MPI_Recv(&vstrSD[i][0], nSizeAll[i], MPI_CHAR, iMPIRank_Ad, iMPIRank_Ad,
				MPI_COMM_WORLD, &status);
		}
	}

	if (iMPIRank == 0) {
		for (int i = 1; i < num_root_proc; ++i) {
			strSD += vstrSD[i];
		}
		multimap<RbtInt, RbtString> recordmap;

		std::stringstream ss(strSD.c_str());
		while (true) {
			string strRec;
			RbtString line;
			RbtInt iRecord = 0;
			while (getline(ss, line, '\n')) {
				if (line == ">  <Record>") {
					getline(ss, line, '\n');
					iRecord = atoi(line.c_str());
					getline(ss, line, '\n');
					continue;
				}
				strRec += line + "\n";
				//if (line == "$$$$") {
				if (line == "END") {
					break;
				}
			}
			if (0 < iRecord) {
				recordmap.insert(pair<RbtInt, RbtString>(iRecord, strRec));
			}
			else {
				break;
			}
		}
		strSD = "";
		multimap<RbtInt, RbtString>::iterator it = recordmap.begin();
		for (; it != recordmap.end(); ++it) {
			strSD += it->second;
		}

		ofstream ofs(strOutput.c_str());
		if (ofs) {
			ofs << strSD;
		}
		ofs.close();
	}
}
#endif // AdvanceSoft(2018)



#if 1 // Xability (2017)
//int main(int argc, const char* argv[]) // original code
int main(int argc, char* argv[])
#endif // Xability (2017)
{

#if 1 // AdvanceSoft (2018)
	std::string buf_str;
	int counter = 0;
	//clock_t time_start_program;
	time_t time_start_program;
	//clock_t time_end_program;
	time_t time_end_program;
#endif // AdvanceSoft(2018)

#if 1 // Xability (2017)
	MPI_Init(&argc, &argv);
	int iMPIRank = 0;
	int nMPISize = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &iMPIRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
	
	ostrstream s;
	s << iMPIRank << ends;
	RbtString strMPIRank = s.str();
	RbtString strSDFile;
#endif // Xability (2017)

#if 1 // AdvanceSoft (2018)
	RbtString strPDBFile;
	bool b_mmgbpbsa(false);
	std::string str_mmgbpbsa_starting_point;
	mmgbpbsa_starting_point = 10;
#ifdef _VISUAL_STUDIO
	if (iMPIRank == 0) {
		//std::cout << "Push enter" << std::endl;
		//fflush(NULL);
		std::cout << "iMPIRank = " << iMPIRank << std::endl;
		std::cout << "nMPISize = " << nMPISize << std::endl;
		std::cout << "enter any charactoer" << std::endl;
		//fflush(NULL);
		std::cin >> buf_str;
	}

	MPI_Barrier(MPI_COMM_WORLD);
#else
#endif
	//time_start_program = clock();
	time_start_program = time(NULL);
#endif // AdvanceSoft(2018)
	
#ifdef _VISUAL_STUDIO 
#else
#if 1 // Xability (2017)
		std::streambuf* sbOut = cout.rdbuf();
		std::streambuf* sbErr = cerr.rdbuf();
		ofstream ofsOut(("proc" + strMPIRank + ".out").c_str());
		ofstream ofsErr(("proc" + strMPIRank + ".err").c_str());
		cout.rdbuf(ofsOut.rdbuf());
		cerr.rdbuf(ofsErr.rdbuf());
#endif // Xability (2017)
#endif



	cout.setf(ios_base::left,ios_base::adjustfield);

	//Strip off the path to the executable, leaving just the file name
	RbtString strExeName(argv[0]);
	RbtString::size_type i = strExeName.rfind("/");
	if (i != RbtString::npos)
		strExeName.erase(0,i+1);

	//Print a standard header
	Rbt::PrintStdHeader(cout,strExeName+EXEVERSION);

	//Command line arguments and default values
	RbtString	strLigandMdlFile;

	RbtBool		bOutput(false);
	RbtString	strRunName;
	RbtString	strReceptorPrmFile;		//Receptor param file
	RbtString	strParamFile;			//Docking run param file
	RbtString       strFilterFile; // Filter file
	
#ifdef _VISUAL_STUDIO 
	RbtInt		nDockingRuns(2);//Init to zero, so can detect later whether user explictly typed -n
#else
	RbtInt		nDockingRuns(0);//Init to zero, so can detect later whether user explictly typed -n
#endif

	//Params for target score
	RbtBool		bTarget(false);
	RbtBool		bStop(true);			//DM 25 May 2001 - if true, stop once target is met
	RbtBool         bDockingRuns(false); // is argument -n present?
	RbtDouble	dTargetScore(0.0);
	RbtBool         bFilter(false);

	RbtBool		bPosIonise(false);
	RbtBool		bNegIonise(false);
	
#ifdef _VISUAL_STUDIO
	RbtBool		bImplH(false);
#else
	RbtBool		bImplH(true);			//if true, read only polar hydrogens from SD file, else read all H's present
#endif

#if 1 // Xability (2017)
	RbtBool		bAttachH(false);		//if true, attach igonored H's present
#endif // Xability (2017)

	RbtBool		bSeed(false);			//Random number seed (default = from system clock)
	RbtInt		nSeed(0);
	RbtBool         bTrace(false);
	RbtInt		iTrace(0);//Trace level, for debugging

	// variables for popt command-line parsing
	char 			c;						// for argument parsing
	
#ifdef _VISUAL_STUDIO
#else
	poptContext		optCon;					// ditto
#endif

#ifdef _VISUAL_STUDIO 
#if 1
	//char 			*inputFile = "1sj0_ligand.sd";		// will be 'strLigandMdlFile'
	char 			*inputFile = "1sj0_ligand_2ligands.sd";		// will be 'strLigandMdlFile'
	//char 			*inputFile = "1sj0_ligand_4ligands.sd";		// will be 'strLigandMdlFile'
	//char 			*inputFile = "mcr_actives_final_10.sdf";		// will be 'strLigandMdlFile'
	char 			*outputFile = "out";		// will be 'strRunName'
	char 			*receptorFile = "CDK2_0.prm";		// will be 'strReceptorPrmFile'
	char 			*protocolFile = "dock.prm";		// will be 'strParamFile' 
#endif
#else
	char 			*inputFile = NULL;		// will be 'strLigandMdlFile'
	char 			*outputFile = NULL;		// will be 'strRunName'
	char 			*receptorFile = NULL;		// will be 'strReceptorPrmFile'
	char 			*protocolFile = NULL;		// will be 'strParamFile'
#endif
	char 			*strTargetScr = NULL;		// will be 'dTargetScore' 

#if 1 // AdvanceSoft (2018)
	int iIsland = 1; // number of processes for islands
#endif // AdvanceSoft(2018)

 
	
#ifdef _VISUAL_STUDIO 
#else
	struct poptOption optionsTable[] = {	// command line options
		{ "input",		'i',POPT_ARG_STRING | POPT_ARGFLAG_ONEDASH,&inputFile,   'i',"input file" },
		{ "output",		'o',POPT_ARG_STRING | POPT_ARGFLAG_ONEDASH,&outputFile,  'o',"output file" },
		{ "receptor",	'r',POPT_ARG_STRING | POPT_ARGFLAG_ONEDASH,&receptorFile,'r',"receptor file" },
		{ "protocol",	'p',POPT_ARG_STRING | POPT_ARGFLAG_ONEDASH,&protocolFile,'p',"protocol file" },
		{ "runs",		'n',POPT_ARG_INT | POPT_ARGFLAG_ONEDASH,&nDockingRuns,'n',"number of runs" },
		{ "trace",		'T',POPT_ARG_INT | POPT_ARGFLAG_ONEDASH,&iTrace,'T',"trace level for debugging" },
		{ "seed",		's',POPT_ARG_INT | POPT_ARGFLAG_ONEDASH,&nSeed,       's',"random seed" },
		{ "ap",			'P',POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH,0,            'P',"protonate groups" },
		{ "an",			'D',POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH,0,            'D',"DEprotonate groups" },
		{ "allH",        'H',POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH,0,            'H',"read all Hs" },
		{ "attachH",     'h',POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH,0,            'h',"ignore Hs, but atatch the Hs to result" },

#if 1 // AdvanceSoft (2018)
		{ "sideChain",  'S',POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH,0,            'S',"side-chains of a receptor" },
		{ "island",     'l',POPT_ARG_INT | POPT_ARGFLAG_ONEDASH,&iIsland,      'l',"number of islands in the island model" },
		{ "migration",  'm',POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH,0,            'm',"migration for the island-model GA" },
		{ "mmgbpbsa",	'g',POPT_ARG_STRING | POPT_ARGFLAG_ONEDASH,&str_mmgbpbsa_starting_point,'g',"str_mmgbpbsa_starting_point" },
		{ "parallel_mmgbpbsa",  'M',POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH,0,    'M',"parallel computation of MM-GB/PBSA rDock" },
		{ "outputRec",  'R',POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH,0,            'R',"output a PDB file including the receptor structure(s)" },
#endif // AdvanceSoft(2018)

		{ "target",      't',POPT_ARG_STRING | POPT_ARGFLAG_ONEDASH,&strTargetScr,'t',"target score" },
		{ "cont",        'C',POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH,0,            'C',"continue even if target met" },
		POPT_AUTOHELP
	{ NULL,0,0,NULL,0 }
	};
#endif

#ifdef _VISUAL_STUDIO 
#else
	optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
	poptSetOtherOptionHelp(optCon, "-r<receptor.prm> -p<protocol.prm> -i<infile> -o<outfile> [options]");
#endif
	
#ifdef _VISUAL_STUDIO 
#else
	//Brief help message
	if (argc < 2) {
		PrintUsage();
		return 1;
	}
	RbtDouble val(0.0);
#endif

#ifdef _VISUAL_STUDIO
#else
	while ((c = poptGetNextOpt(optCon)) >= 0) {
		switch (c) {
		case 'P':	// protonate
			bPosIonise = true;
			break;
		case 'D':	// deprotonate
			bNegIonise = true;
			break;
		case 'H':	// all-H model
			bImplH = false;
			break;

#if 1 // Xability (2017)
		case 'h':	// attach H
			bAttachH = true;
			break;
#endif // Xability (2017)

#if 1 // AdvanceSoft (2018)
		case 'S':	// side-chains of a receptor
			b_side_chains_of_receptor = true;
			break;
		case 'm':	// migration
			b_islandGA_w_migration = true;
			initial_b_islandGA_w_migration = true;
			break;
		case 'g':
			b_mmgbpbsa = true;
			break;
		case 'M':
			b_parallel_mmgbpbsa = true;
			break;
		case 'R':	// output a receptor PDB file
			b_receptor = true;
			break;
#endif // AdvanceSoft(2018)

		case 'C':
			bStop = false;
			break;
		case 't':
			// If str can be translated to an integer, I assume is a
			// threshold. Otherwise, I assume is the filter file name
			char *error;
			errno = 0;
			val = strtod(strTargetScr, &error);
			if (!errno && !*error)  // Is it a number?
			{
				dTargetScore = val;
				bTarget = true;
			}
			else // Assume is the filter file name
			{
				bFilter = true;
				strFilterFile = strTargetScr;
			}
			break;
		case 's':
			bSeed = true;
			break;
		case 'T':
			bTrace = true;
			break;
		default:
			break;
		}
	}
	cout << endl;
	poptFreeContext(optCon);
#endif
	
#if 1 // AdvanceSoft (2018)
	if (b_mmgbpbsa) {
		bImplH = false;
	}
#endif // AdvanceSoft (2018)

	// print out arguments
	// input ligand file, receptor and parameter is compulsory
	cout << endl << "Command line args:" << endl;
	if(!inputFile || !receptorFile || !protocolFile) { // if any of them is missing
		
#ifdef _VISUAL_STUDIO 
		std::cout << "the arguments given are not correct" << std::endl;
#else
		poptPrintUsage(optCon, stderr, 0);
#endif	
		
		exit(1);
	} else {
		strLigandMdlFile	= inputFile;
		strReceptorPrmFile	= receptorFile;
		strParamFile		= protocolFile;
		cout << " -i " << strLigandMdlFile		<< endl;
		cout << " -r " << strReceptorPrmFile	<< endl;
		cout << " -p " << strParamFile			<< endl;

#if 1 // AdvanceSoft (2018)
		if (b_side_chains_of_receptor) {
			std::cout << " -sideChain " << b_side_chains_of_receptor << std::endl;
		}
		if (iIsland) {
			std::cout << " -l " << iIsland << std::endl;
		}
		if (b_islandGA_w_migration) {
			std::cout << " -m " << b_islandGA_w_migration << std::endl;
		}
		if (b_mmgbpbsa) {
			std::cout << " -g " << b_mmgbpbsa << std::endl;
		}
		if (b_parallel_mmgbpbsa) {
			std::cout << " -M " << b_parallel_mmgbpbsa << std::endl;
		}
		if (b_receptor) {
			std::cout << " -outputRec " << b_receptor << std::endl;
		}
#endif // AdvanceSoft(2018)

	}
	// output is not that important but good to have
	if(outputFile) {
		strRunName	= outputFile;

#if 1 // Xability (2017)
		strSDFile = strRunName + "_proc" + strMPIRank + ".sd";
#endif // Xability (2017)

#if 1 // AdvanceSoft (2018)
		strPDBFile = strRunName + "_receptor_proc" + strMPIRank + ".pdb";
#endif // AdvanceSoft(2018)



		bOutput = true;
		cout << " -o " << strRunName << endl;
	} else {
		cout << "WARNING: output file name is missing."<< endl;
	}
	// docking runs
	if(nDockingRuns >= 1) {//User typed -n explicitly, so set bDockingRuns to true
	  bDockingRuns = true;
	  cout << " -n " << nDockingRuns << endl;
	} else {
	  nDockingRuns = 1;//User didn't type -n explicitly, so fall back to the default of n=1
	  cout << " -n " << nDockingRuns << " (default) " <<endl;
	}
	if(bSeed)		// random seed (if provided)
		cout << " -s " << nSeed << endl;
	if(bTrace)		// random seed (if provided)
		cout << " -T " << iTrace << endl;
	if(bPosIonise)	// protonate 
		cout << " -ap " << endl;
	if(bNegIonise)	// deprotonate
		cout << " -an " << endl;
	if(!bImplH)		// all-H
		cout << " -allH " << endl;

#if 1 // Xability (2017)
	if (bAttachH)		// all-H
		cout << " -attachH " << endl;
#endif // Xability (2017)

	if(!bStop)		// stop after target
		cout << " -cont " << endl;
	if(bTarget)
		cout << " -t " << dTargetScore << endl;

  //BGD 26 Feb 2003 - Create filters to simulate old rbdock
  //behaviour
  ostrstream strFilter;
  if (!bFilter)
  {
    if (bTarget) // -t<TS>
    {
      if (!bDockingRuns) // -t<TS> only
      {
        strFilter << "0 1 - SCORE.INTER " << dTargetScore << endl;
      }
      else  // -t<TS> -n<N> need to check if -cont present
            // for all other cases it doesn't matter
        if (!bStop) // -t<TS> -n<N> -cont
        {
          strFilter << "1 if - SCORE.NRUNS " << (nDockingRuns - 1)
                    << " 0.0 -1.0,\n1 - SCORE.INTER " << dTargetScore << endl;
        }
        else // -t<TS> -n<N>
        {
          strFilter << "1 if - " << dTargetScore << " SCORE.INTER 0.0 " 
                    << "if - SCORE.NRUNS " << (nDockingRuns - 1) 
                    << " 0.0 -1.0,\n1 - SCORE.INTER " << dTargetScore << endl;
        }
    } // no target score, no filter
    else if (bDockingRuns) // -n<N>
    {
      strFilter << "1 if - SCORE.NRUNS " << (nDockingRuns - 1) 
                << " 0.0 -1.0,\n0";
    }
    else // no -t no -n
    {
      strFilter << "0 0\n";
    }
  }

  //DM 20 Apr 1999 - set the auto-ionise flags
  if (bPosIonise)
    cout << "Automatically protonating positive ionisable groups (amines, imidazoles, guanidines)" << endl;
  if (bNegIonise)
    cout << "Automatically deprotonating negative ionisable groups (carboxylic acids, phosphates, sulphates, sulphonates)" << endl;
  if (bImplH)
    cout << "Reading polar hydrogens only from ligand SD file" << endl;
  else
    cout << "Reading all hydrogens from ligand SD file" << endl;

  if (bTarget) {
    cout << endl << "Lower target intermolecular score = " << dTargetScore << endl;
  }


#if 1 // AdvanceSoft (2018)
  int iSdf = 0; // number of processes for sdf files

#ifdef _VISUAL_STUDIO
  b_receptor = true;
  b_side_chains_of_receptor = 0;
  b_parallel_mmgbpbsa = true;
  buf_str = *argv[1];
  iIsland = std::atoi(buf_str.c_str());
  buf_str = *argv[2];
  b_islandGA_w_migration = std::atoi(buf_str.c_str());
  initial_b_islandGA_w_migration = b_islandGA_w_migration;
#else
#endif

  iSdf = nMPISize / iIsland;
  RbtBaseTransform::Set_m_iSdf(iSdf);
  RbtBaseTransform::Set_m_iIsland(iIsland);
  if (b_mmgbpbsa) {
	  mmgbpbsa_starting_point = std::atoi(str_mmgbpbsa_starting_point.c_str());
  }

  if (b_parallel_mmgbpbsa) {
	  iSdf = 1;
	  iIsland = 1;
	  RbtBaseTransform::Set_m_iSdf(iSdf);
	  RbtBaseTransform::Set_m_iIsland(iIsland);
  }
#endif // AdvanceSoft(2018)

#if 1 // debug write
  std::cout << std::endl;
  std::cout << "bPosIonise = " << bPosIonise << endl;
  std::cout << "bNegIonise = " << bNegIonise << endl;
  std::cout << "bImplH = " << bImplH << endl;
  std::cout << "bStop = " << bStop << endl;
  std::cout << "dTargetScore = " << dTargetScore << endl;
  std::cout << "bTarget = " << bTarget << endl;
  std::cout << "bFilter = " << bFilter << endl;
  std::cout << "strFilterFile = " << strFilterFile << endl;
  std::cout << "bSeed = " << bSeed << endl;
  std::cout << "bTrace = " << bTrace << endl;
  std::cout << "strLigandMdlFile = " << strLigandMdlFile << endl;
  std::cout << "strReceptorPrmFile = " << strReceptorPrmFile << endl;
  std::cout << "strParamFile = " << strParamFile << endl;
  std::cout << "strRunName = " << strRunName << endl;
  std::cout << "bOutput = " << bOutput << endl;
  std::cout << "nDockingRuns = " << nDockingRuns << endl;
  std::cout << "bDockingRuns = " << bDockingRuns << endl;
  std::cout << "b_side_chains_of_receptor = " << b_side_chains_of_receptor << std::endl;
  std::cout << "iIsland = " << iIsland << std::endl;
  std::cout << "b_islandGA_w_migration = " << b_islandGA_w_migration << std::endl;
  std::cout << "iSdf = " << iSdf << std::endl;
  std::cout << "b_receptor = " << b_receptor << std::endl;
  std::cout << "b_mmgbpbsa = " << b_mmgbpbsa << endl;
  std::cout << "mmgbpbsa_starting_point = " << mmgbpbsa_starting_point << endl;
  std::cout << "b_parallel_mmgbpbsa = " << b_parallel_mmgbpbsa << endl;
#endif // debug write



  try {
    //Create a bimolecular workspace
    RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
    
	//Set the workspace name to the root of the receptor .prm filename
    RbtStringList componentList = Rbt::ConvertDelimitedStringToList(strReceptorPrmFile,".");
    RbtString wsName = componentList.front();
    spWS->SetName(wsName);
    
    //Read the docking protocol parameter file
    RbtParameterFileSourcePtr spParamSource(new RbtParameterFileSource(Rbt::GetRbtFileName("data/scripts",strParamFile)));

	//Read the receptor parameter file
    RbtParameterFileSourcePtr spRecepPrmSource(new RbtParameterFileSource(Rbt::GetRbtFileName("data/receptors",strReceptorPrmFile)));
    cout << endl << "DOCKING PROTOCOL:" << endl << spParamSource->GetFileName() << endl << spParamSource->GetTitle() << endl;
    cout << endl << "RECEPTOR:" << endl << spRecepPrmSource->GetFileName() << endl << spRecepPrmSource->GetTitle() << endl;

    //Create the scoring function from the SCORE section of the docking protocol prm file
    //Format is:
    //SECTION SCORE
    //    INTER    RbtInterSF.prm
    //    INTRA RbtIntraSF.prm
    //END_SECTION
    //
    //Notes:
    //Section name must be SCORE. This is also the name of the root SF aggregate
    //An aggregate is created for each parameter in the section.
    //Parameter name becomes the name of the subaggregate (e.g. SCORE.INTER)
    //Parameter value is the file name for the subaggregate definition
    //Default directory is $RBT_ROOT/data/sf
    RbtSFFactoryPtr spSFFactory(new RbtSFFactory());//Factory class for scoring functions
    RbtSFAggPtr spSF(new RbtSFAgg(_ROOT_SF));//Root SF aggregate
    spParamSource->SetSection(_ROOT_SF);
    RbtStringList sfList(spParamSource->GetParameterList());
    //Loop over all parameters in the SCORE section
    for (RbtStringListConstIter sfIter = sfList.begin(); sfIter != sfList.end(); sfIter++) {
      //sfFile = file name for scoring function subaggregate
      RbtString sfFile(Rbt::GetRbtFileName("data/sf",spParamSource->GetParameterValueAsString(*sfIter)));
      RbtParameterFileSourcePtr spSFSource(new RbtParameterFileSource(sfFile));
      //Create and add the subaggregate
      spSF->Add(spSFFactory->CreateAggFromFile(spSFSource,*sfIter));
    }
    
    //Add the RESTRAINT subaggregate scoring function from any SF definitions in the receptor prm file
    spSF->Add(spSFFactory->CreateAggFromFile(spRecepPrmSource,_RESTRAINT_SF));
    
    //Create the docking transform aggregate from the transform definitions in the docking prm file
    RbtTransformFactoryPtr spTransformFactory(new RbtTransformFactory());
    spParamSource->SetSection();
    RbtTransformAggPtr spTransform(spTransformFactory->CreateAggFromFile(spParamSource,_ROOT_TRANSFORM));
    
    //Override the TRACE levels for the scoring function and transform
    //Dump details to cout
    //Register the scoring function and the transform with the workspace
    if (bTrace) {
      RbtRequestPtr spTraceReq(new RbtSFSetParamRequest("TRACE",iTrace));
      spSF->HandleRequest(spTraceReq);
      spTransform->HandleRequest(spTraceReq);
    }
    if (iTrace > 0) {
      cout << endl << "SCORING FUNCTION DETAILS:" << endl << *spSF << endl;
      cout << endl << "SEARCH DETAILS:" << endl << *spTransform << endl;
    }
    spWS->SetSF(spSF);
    spWS->SetTransform(spTransform);
    
    //DM 18 May 1999
    //Variants describing the library version, exe version, parameter file, and current directory
    //Will be stored in the ligand SD files
    RbtVariant vLib(Rbt::GetProduct()+" ("+Rbt::GetVersion()+", Build"+Rbt::GetBuild()+")");
    RbtVariant vExe(strExeName+EXEVERSION);
    RbtVariant vRecep(spRecepPrmSource->GetFileName());
    RbtVariant vPrm(spParamSource->GetFileName());
    RbtVariant vDir(Rbt::GetCurrentDirectory());
    
    spRecepPrmSource->SetSection();
   //Read docking site from file and register with workspace
    RbtString strASFile = spWS->GetName()+".as";
    RbtString strInputFile = Rbt::GetRbtFileName("data/grids",strASFile);
    //DM 26 Sep 2000 - ios_base::binary is invalid with IRIX CC
#if defined(__sgi) && !defined(__GNUC__)
    ifstream istr(strInputFile.c_str(),ios_base::in);
#else
    ifstream istr(strInputFile.c_str(),ios_base::in|ios_base::binary);
#endif
    //DM 14 June 2006 - bug fix to one of the longest standing rDock issues
    //(the cryptic "Error reading from input stream" message, if cavity file was missing) 
    if (!istr) {
      RbtString message = "Cavity file (" + strASFile + ") not found in current directory or $RBT_HOME";
      message += " - run rbcavity first"; 
      throw RbtFileReadError(_WHERE_,message);
    } 
    RbtDockingSitePtr spDS(new RbtDockingSite(istr));
    istr.close();
    spWS->SetDockingSite(spDS);
    cout << endl << "DOCKING SITE" << endl << (*spDS) << endl;

    //Prepare the SD file sink for saving the docked conformations for each ligand
    //DM 3 Dec 1999 - replaced ostrstream with RbtString in determining SD file name
    //SRC 2014 moved here this block to allow WRITE_ERROR TRUE
    if (bOutput) {
      
#if 1 // Xability (2017)
		//RbtMolecularFileSinkPtr spMdlFileSink(new RbtMdlFileSink(strRunName+".sd",RbtModelPtr()));
		RbtMolecularFileSinkPtr spMdlFileSink(new RbtMdlFileSink(strSDFile, RbtModelPtr()));
#endif // Xability (2017)

		spWS->SetSink(spMdlFileSink);

#if 1 // AdvanceSoft (2018)
		RbtMolecularFileSinkPtr spMdlFileSink_rec(new RbtMdlFileSink(strPDBFile, RbtModelPtr()));
		spWS->SetSink_rec(spMdlFileSink_rec);
#endif // AdvanceSoft(2018)

    }
    
    RbtPRMFactory prmFactory(spRecepPrmSource, spDS);
    prmFactory.SetTrace(iTrace);

	//Create the receptor model from the file names in the receptor parameter file
    RbtModelPtr spReceptor = prmFactory.CreateReceptor();
    spWS->SetReceptor(spReceptor);

    //Register any solvent
    RbtModelList solventList = prmFactory.CreateSolvent();
 	spWS->SetSolvent(solventList);
 	if (spWS->hasSolvent()) {
 		RbtInt nSolvent = spWS->GetSolvent().size();
 		cout << endl << nSolvent << " solvent molecules registered" << endl;
 	}
 	else {
 		cout << endl << "No solvent" << endl;
 	}
 	
    //SRC 2014 removed sector bOutput from here to some blocks above, for WRITEERRORS TRUE
    
    //Seed the random number generator
    RbtRand& theRand = Rbt::GetRbtRand();//ref to random number generator
    if (bSeed) {
      theRand.Seed(nSeed);
    }
    
    //Create the filter object for controlling early termination of protocol
    RbtFilterPtr spfilter;
    if (bFilter) {
      spfilter = new RbtFilter(strFilterFile);
      if (bDockingRuns) {
	spfilter->SetMaxNRuns(nDockingRuns);
      }
    }
    else {
      spfilter = new RbtFilter(strFilter.str(), true);
    }
    if (bTrace) {
      RbtRequestPtr spTraceReq(new RbtSFSetParamRequest("TRACE",iTrace));
      spfilter->HandleRequest(spTraceReq);
    }

    //Register the Filter with the workspace
    spWS->SetFilter(spfilter);

    //MAIN LOOP OVER LIGAND RECORDS
    //DM 20 Apr 1999 - add explicit bPosIonise and bNegIonise flags to MdlFileSource constructor
    RbtMolecularFileSourcePtr spMdlFileSource(new RbtMdlFileSource(strLigandMdlFile,bPosIonise,bNegIonise,bImplH));

#if 1 // Xability (2017)
	RbtMolecularFileSourcePtr spMdlFileSourceAtc;
	if (bImplH && bAttachH) {
		spMdlFileSourceAtc = new RbtMdlFileSource(strLigandMdlFile, bPosIonise, bNegIonise, false);
	}
#endif // Xability (2017)
	
	for (RbtInt nRec=1; spMdlFileSource->FileStatusOK(); spMdlFileSource->NextRecord(), nRec++) {
      
#if 1 // Xability (2017)
		if (!spMdlFileSourceAtc.Null()) {
			if (1 < nRec) {
				if (!spMdlFileSourceAtc->FileStatusOK()) {
					break;
				}
				spMdlFileSourceAtc->NextRecord();
			}
		}
#endif // Xability (2017)

#if 0 // Xability (2017)
		if ((nRec - 1) % nMPISize != iMPIRank) {
			continue;
		}
#endif // Xability (2017)

#if 1 // AdvanceSoft (2018)
		if (nMPISize == 1) {}
		else if (b_parallel_mmgbpbsa) {}
		else {
			bool b_branch(false);
			for (int i = 0; i < iIsland; ++i) {
				if (((nRec - 1) % iSdf) * iIsland + i == iMPIRank) {
					b_branch = true;
					break;
				}
			}
			if (!b_branch) {
				continue;
			}
		}
#endif // AdvanceSoft(2018)



		cout.setf(ios_base::left,ios_base::adjustfield);
      cout << endl
       << "**************************************************" << endl
       << "RECORD #" << nRec << endl;
      RbtError molStatus = spMdlFileSource->Status();
      if (!molStatus.isOK()) {
        cout << endl << molStatus << endl
             << "************************************************" << endl;
        continue;
      }
      
      //DM 26 Jul 1999 - only read the largest segment (guaranteed to be called H)
      //BGD 07 Oct 2002 - catching errors created by the ligands,
      //so rbdock continues with the next one, instead of
      //completely stopping
      try
      {
        spMdlFileSource->SetSegmentFilterMap
                                  (Rbt::ConvertStringToSegmentMap("H"));
      
        if (spMdlFileSource->isDataFieldPresent("Name"))
          cout << "NAME:   " << spMdlFileSource->GetDataValue("Name") << endl;
        if (spMdlFileSource->isDataFieldPresent("REG_Number"))
          cout << "REG_Num:" << spMdlFileSource->GetDataValue("REG_Number") 
               << endl;



#if 1 // AdvanceSoft (2018)
		if (m_iIsland > 1) {
			srand((unsigned int)time(NULL));
			theRand.Seed(rand() + iMPIRank * nMPISize);
		}
		else {}
#else
		theRand.Seed(1);
#endif // AdvanceSoft(2018)

        cout << setw(30) << "RANDOM_NUMBER_SEED:" << theRand.GetSeed() << endl;

        //Create and register the ligand model
        RbtModelPtr spLigand = prmFactory.CreateLigand(spMdlFileSource);
        
#if 1 // Xability (2017)
		RbtModelPtr spLigandAtc;
		RbtCoordList coordAtc;
		RbtIntList viAtcMap;
		if (!spMdlFileSourceAtc.Null()) {
			spLigandAtc = prmFactory.CreateLigand(spMdlFileSourceAtc);
			RbtAtomList atomListAtc = spLigandAtc->GetAtomList();
			RbtInt nAtomAtc = atomListAtc.size();
			coordAtc.resize(nAtomAtc);
			for (RbtInt i = 0; i < nAtomAtc; ++i) {
				coordAtc[i] = atomListAtc[i]->GetCoords();
			}
			viAtcMap.resize(nAtomAtc, -1);
			RbtInt i = 0;
			RbtAtomList atomList = spLigand->GetAtomList();
			for (RbtInt j = 0; j < atomList.size(); ++j) {
				RbtCoord vj = atomList[j]->GetCoords();
				for (; i < nAtomAtc; ++i) {
					RbtCoord vi = coordAtc[i];
					if (Rbt::Length2(vi, vj) < 1e-10) {
						viAtcMap[i] = j;
						break;
					}
				}
			}
		}
#endif // Xability (2017)

		RbtString strMolName = spLigand->GetName();
        spWS->SetLigand(spLigand);
        
#if 1 // Xability (2017)
		spWS->SetLigandAtc(spLigandAtc);
#endif // Xability (2017)
		
		//Update any model coords from embedded chromosomes in the ligand file
        spWS->UpdateModelCoordsFromChromRecords(spMdlFileSource, iTrace);
      
        //DM 18 May 1999 - store run info in model data
        //Clear any previous Rbt.* data fields
        spLigand->ClearAllDataFields("Rbt.");
        spLigand->SetDataValue("Rbt.Library",vLib);
        spLigand->SetDataValue("Rbt.Executable",vExe);
        spLigand->SetDataValue("Rbt.Receptor",vRecep);
        spLigand->SetDataValue("Rbt.Parameter_File",vPrm);
        spLigand->SetDataValue("Rbt.Current_Directory",vDir);
      
        //DM 10 Dec 1999 - if in target mode, loop until target score is reached
        RbtBool bTargetMet = false;
        
        ////////////////////////////////////////////////////
        //MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
        //Create a history file sink, just in case it's needed by any 
        //of the transforms
        RbtInt iRun = 1;
	// need to check this here. The termination 
	// filter is only run once at least
	// one docking run has been done.
        if (nDockingRuns < 1) 	
          bTargetMet = true;
        while (!bTargetMet) {
	  //Catching errors with this specific run
          try {
	    if (bOutput) {
              ostrstream histr;
              histr << strRunName << "_" << strMolName << nRec << "_his_" 
                    << iRun << ".sd" << ends;
              RbtMolecularFileSinkPtr spHistoryFileSink
		(new RbtMdlFileSink(histr.str(),spLigand));
              delete histr.str();
              spWS->SetHistorySink(spHistoryFileSink);
            }



#if 0 // debug write
		std::cout << std::endl;
		std::cout << "iRun = " << iRun << std::endl;
		std::cout << "*************** Chromosome elements ***************" << std::endl;
		std::cout << std::endl;
		for (int i = 0; i < spWS->GetNumModels(); i++) {

			std::cout << "swWS->GetModel(" << i << ")->GetChrom(): \n" 
				<< *(spWS->GetModel(i)->GetChrom()) << std::endl;
		}
#endif // debug write



		spWS->Run();//Dock!



	    RbtBool bterm = spfilter->Terminate();
	    RbtBool bwrite = spfilter->Write();
	    if (bterm)
	      bTargetMet = true;
	    if (bOutput && bwrite) {

#if 1 // Xability (2017)
			if (!spLigandAtc.Null()) {
				spLigandAtc->AlignCoord(*spLigand, viAtcMap, coordAtc);
			}
#endif // Xability (2017)

	      spWS->Save();
	    }
	    iRun++;
#if 1 // AdvanceSoft (2018)
		counter_GA = 0;
#endif // AdvanceSoft (2018)
          }
          catch (RbtDockingError& e) {
	    cout << e << endl;
	  }
        }

#if 1 // AdvanceSoft (2018)
		num_docking_Ad = 0;
#endif // AdvanceSoft(2018)

	//END OF MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
	////////////////////////////////////////////////////
      } 
      //END OF TRY
      catch (RbtLigandError& e) {
	cout << e << endl;
      }
    }	//	the end of the "for" sentence

    //END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////
    cout << endl << "END OF RUN" << endl;
//    if (bOutput && flexRec) {
//      RbtMolecularFileSinkPtr spRecepSink(new RbtCrdFileSink(strRunName+".crd",spReceptor));
//      spRecepSink->Render();
//    }
	
#if 1 // AdvanceSoft (2018)
	if (b_parallel_mmgbpbsa) {
		if (iMPIRank == 0) {
			ConcatenateSDFiles(strSDFile, strRunName + ".sd");
			if (b_receptor) {
				ConcatenatePDBFiles(strPDBFile, strRunName + "_receptor.pdb");
			}
		}
	}
	else if (iMPIRank % m_iIsland == 0) {
		ConcatenateSDFiles(strSDFile, strRunName + ".sd");
		if (b_receptor) {
			ConcatenatePDBFiles(strPDBFile, strRunName + "_receptor.pdb");
		}
	}
#endif // AdvanceSoft(2018)

  }	//	the end of the first "try"



  catch (RbtError& e) {
    cout << e << endl;
  }
  catch (...) {
    cout << "Unknown exception" << endl;
  }

#if 1 // debug write
  //time_end_program = clock();
  time_end_program = time(NULL);
  //std::cout << "duration of the program = " << (double)(time_end_program - time_start_program) / CLOCKS_PER_SEC << "sec.\n";
  std::cout << "duration of the program = " << time_end_program - time_start_program << "sec.\n";
  std::cout << "counter_iCycle = " << counter_iCycle << std::endl;
  std::cout << "counter_score_calculations = " << counter_score_calculations << std::endl;
  std::cout << "counter_mmgbpbsa_calculations = " << counter_mmgbpbsa_calculations << std::endl;
#endif // debug write

#ifdef _VISUAL_STUDIO 
#else
#if 1 // Xability (2017)
  ofsOut.close();
  ofsErr.close();
  cout.rdbuf(sbOut);
  cerr.rdbuf(sbErr);
#endif // Xability (2017)
#endif

#if 1 // AdvanceSoft (2018)
  if (b_mmgbpbsa) {
	  buf_str = "rm -r dir_proc_" + RbtWorkSpace::to_string(iMPIRank);
	  system(buf_str.c_str());
  }
#endif // AdvanceSoft (2018)

  MPI_Finalize();

  _RBTOBJECTCOUNTER_DUMP_(cout)

  return 0;
}
