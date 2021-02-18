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

#include "RbtGenome.h"
#include "RbtBaseSF.h"
#include "RbtChromElement.h"

#if 1 // AdvanceSoft (2018)
#include "RbtWorkSpace.h"
#include <time.h>
#include "RbtBaseTransform.h"
#include <mpi.h>
#include <iomanip>
int counter_score_calculations = 0;
int counter_mmgbpbsa_calculations = 0;
#endif // AdvanceSoft (2018)

#ifdef _VISUAL_STUDIO
#include <chrono>
#endif

RbtString RbtGenome::_CT("RbtGenome");

RbtGenome::RbtGenome(RbtChromElement* pChr)
        : m_chrom(pChr->clone()),
        m_score(0.0),
        m_RWFitness(0.0)
{
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtGenome::RbtGenome(const RbtGenome& g)
        : m_chrom( (g.m_chrom)->clone() ),
        m_score(g.m_score),
        m_RWFitness(g.m_RWFitness)
{
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtGenome& RbtGenome::operator=(const RbtGenome& g) {
  if (&g != this) {
    m_chrom = (g.m_chrom)->clone();
    m_score = g.m_score;
    m_RWFitness = g.m_RWFitness;
  }
  return *this;
}

RbtGenome::~RbtGenome() {
    delete m_chrom;
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

RbtGenome* RbtGenome::clone() const {
    return new RbtGenome(*this);
}

void RbtGenome::SetScore(RbtBaseSF* pSF) {
  if (pSF != NULL) {
    m_chrom->SyncToModel();

#if 1 // AdvanceSoft (2018)

	int iMPIRank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &iMPIRank);

	counter_score_calculations++;

	std::string buf_str;
	std::stringstream buf_stringstream;
	time_t time_start_mmgbpbsa;
	time_t time_end_mmgbpbsa;
	clock_t time_start_temp, time_end_temp;
	std::string energy;
	int net_charge = 0;
	int counter = 0;
	double temp_double = 0;

	if (counter_GA >= mmgbpbsa_starting_point) {

		counter_mmgbpbsa_calculations++;

		time_start_mmgbpbsa = time(NULL);

		static bool b_first(true);
		if (b_first) {
			
			buf_str = "mkdir dir_proc_" + RbtWorkSpace::to_string(iMPIRank);
			system(buf_str.c_str());
			//buf_str = "cp to_after_Rdock_argv_5_br171114.py ./dir_proc_" + RbtWorkSpace::to_string(iMPIRank);
			//system(buf_str.c_str());
			buf_str = "cp conv_Ad.py ./dir_proc_" + RbtWorkSpace::to_string(iMPIRank);
			system(buf_str.c_str());
			buf_str = "cp mmpbsa.in ./dir_proc_" + RbtWorkSpace::to_string(iMPIRank);
			system(buf_str.c_str());
			buf_str = "cp CDK2_0.pdb ./dir_proc_" + RbtWorkSpace::to_string(iMPIRank);
			system(buf_str.c_str());

						buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/shell_script_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".sh";
			std::ofstream ofs_shell(buf_str.c_str());

			buf_str = "cd ./dir_proc_" + RbtWorkSpace::to_string(iMPIRank);
			ofs_shell << buf_str << std::endl;

#if 0
			buf_str = "python to_after_Rdock_argv_5_br171114.py CDK2_0.pdb coordlist_inRbtGenome_receptor_proc_" +
				RbtWorkSpace::to_string(iMPIRank) + " atomlist_inRbtGenome_receptor_proc_" + RbtWorkSpace::to_string(iMPIRank) +
				" temp_receptor_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".pdb";
			ofs_shell << buf_str << std::endl;
#endif

			buf_str = "python conv_Ad.py temp_receptor_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".pdb " + RbtWorkSpace::to_string(iMPIRank);
			ofs_shell << buf_str << std::endl;

			buf_str = "antechamber -i 1_ligand_proc_" + RbtWorkSpace::to_string(iMPIRank) +
				".sd -fi sdf -o lig_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".mol2 -fo mol2 -c gas -nc " +
				RbtWorkSpace::to_string(net_charge) + " -at gaff2 -rn LIG -pf y";
			ofs_shell << buf_str << std::endl;

			buf_str = "parmchk -i lig_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".mol2 -f mol2 -o lig_proc_" +
				RbtWorkSpace::to_string(iMPIRank) + ".frcmod";
			ofs_shell << buf_str << std::endl;

			buf_str = "tleap_rDock -f tleap_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".in";
			ofs_shell << buf_str << std::endl;

			buf_str = "MMPBSA.py -O -i mmpbsa.in -cp com.leap_proc_" + RbtWorkSpace::to_string(iMPIRank)
				+ ".prm -rp protein_amb_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".leap.prm -lp lig.leap_proc_" + RbtWorkSpace::to_string(iMPIRank)
				+ ".prm -y com.leap_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".crd -o FINAL_RESULTS_MMPBSA_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".dat";
			ofs_shell << buf_str << std::endl;

			buf_str = "MMPBSA.py --clean";
			ofs_shell << buf_str << std::endl;

#if 1
			// delete files to prevent them from interfering the next calculation
			buf_str = "rm 1_ligand_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".sd";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm lig_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".mol2";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm lig_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".frcmod";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm protein_amb_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".pdb";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm com.leap_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".prm";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm protein_amb_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".leap.prm";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm lig.leap_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".prm";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm com.leap_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".crd";
			ofs_shell << buf_str << std::endl;

#if 0 // old
			buf_str = "rm atomlist_inRbtGenome_receptor_proc_" + RbtWorkSpace::to_string(iMPIRank);
			ofs_shell << buf_str << std::endl;

			buf_str = "rm atomlist_inRbtGenome_ligand_proc_" + RbtWorkSpace::to_string(iMPIRank);
			ofs_shell << buf_str << std::endl;

			buf_str = "rm coordlist_inRbtGenome_receptor_proc_" + RbtWorkSpace::to_string(iMPIRank);
			ofs_shell << buf_str << std::endl;

			buf_str = "rm coordlist_inRbtGenome_ligand_proc_" + RbtWorkSpace::to_string(iMPIRank);
			ofs_shell << buf_str << std::endl;

			buf_str = "rm bondlist_inRbtGenome_ligand_proc_" + RbtWorkSpace::to_string(iMPIRank);
			ofs_shell << buf_str << std::endl;
#endif // old

			buf_str = "rm com.leap_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".pdb";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm lig.leap_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".crd";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm protein_amb_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".leap.crd";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm temp_receptor_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".pdb";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm leap.log";
			ofs_shell << buf_str << std::endl;

			buf_str = "rm reference.frc";
			ofs_shell << buf_str << std::endl;
#endif

			ofs_shell.close();

			buf_str = "chmod +x ./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/shell_script_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".sh";
			system(buf_str.c_str());

			// create tleap.in
			buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/tleap_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".in";
			std::ofstream ofs_tleap(buf_str.c_str());
			ofs_tleap << "source leaprc.protein.ff14SB" << std::endl;
			ofs_tleap << "source leaprc.gaff2" << std::endl;
			ofs_tleap << "set default pbradii mbondi" << std::endl;
			ofs_tleap << "lig = loadmol2 lig_proc_" << iMPIRank << ".mol2" << std::endl;
			ofs_tleap << "frcmod = loadamberparams lig_proc_" << iMPIRank << ".frcmod" << std::endl;
			ofs_tleap << "saveamberparm lig lig.leap_proc_" << iMPIRank << ".prm lig.leap_proc_" << iMPIRank << ".crd" << std::endl;
			ofs_tleap << "prot = loadpdb protein_amb_proc_" << iMPIRank << ".pdb" << std::endl;
			ofs_tleap << "saveamberparm prot protein_amb_proc_" << iMPIRank << ".leap.prm protein_amb_proc_" << iMPIRank << ".leap.crd" << std::endl;
			ofs_tleap << "complex = combine {prot lig}" << std::endl;
			ofs_tleap << "saveamberparm complex com.leap_proc_" << iMPIRank << ".prm com.leap_proc_" << iMPIRank << ".crd" << std::endl;
			ofs_tleap << "savepdb complex com.leap_proc_" << iMPIRank << ".pdb" << std::endl;
			ofs_tleap << "quit" << std::endl;
			ofs_tleap.close();

			b_first = false;
		}



#if 0 // output information of the receptor and ligand.
		buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/atomlist_inRbtGenome_receptor_proc_" + RbtWorkSpace::to_string(iMPIRank);
		std::ofstream ofs_atomlist_inRbtGenome_receptor(buf_str.c_str());
		if (!ofs_atomlist_inRbtGenome_receptor) {
			exit(1);
		}
		for (int i = 0; i < pSF->GetWorkSpace()->GetModel(0)->GetAtomList().size(); ++i) {
			ofs_atomlist_inRbtGenome_receptor << *pSF->GetWorkSpace()->GetModel(0)->GetAtomList().at(i) << std::endl;
		}
		ofs_atomlist_inRbtGenome_receptor.close();

		buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/coordlist_inRbtGenome_receptor_proc_" + RbtWorkSpace::to_string(iMPIRank);
		std::ofstream ofs_coordlist_inRbtGenome_receptor(buf_str.c_str());
		for (int i = 0; i < pSF->GetWorkSpace()->GetModel(0)->GetAtomList().size(); i++) {
			ofs_coordlist_inRbtGenome_receptor << "ID" << i + 1 << ", " << pSF->GetWorkSpace()->GetModel(0)->GetAtomList().at(i)->GetX() << ", "
				<< pSF->GetWorkSpace()->GetModel(0)->GetAtomList().at(i)->GetY() << ", "
				<< pSF->GetWorkSpace()->GetModel(0)->GetAtomList().at(i)->GetZ() << std::endl;
		}
		ofs_coordlist_inRbtGenome_receptor.close();



		RbtInt ri = pSF->GetWorkSpace()->GetModel(0)->GetCurrentCoords();
		pSF->GetWorkSpace()->GetModel(1)->SetDataValue("RI", ri);

		buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/atomlist_inRbtGenome_ligand_proc_" + RbtWorkSpace::to_string(iMPIRank);
		std::ofstream ofs_atomlist_inRbtGenome_ligand(buf_str.c_str());
		for (int i = 0; i < pSF->GetWorkSpace()->GetModel(1)->GetAtomList().size(); ++i) {
			ofs_atomlist_inRbtGenome_ligand << *pSF->GetWorkSpace()->GetModel(1)->GetAtomList().at(i) << std::endl;
		}
		ofs_atomlist_inRbtGenome_ligand.close();

		buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/bondlist_inRbtGenome_ligand_proc_" + RbtWorkSpace::to_string(iMPIRank);
		std::ofstream ofs_bondlist_inRbtGenome_ligand(buf_str.c_str());
		for (int i = 0; i < pSF->GetWorkSpace()->GetModel(1)->GetBondList().size(); ++i) {
			ofs_bondlist_inRbtGenome_ligand << *pSF->GetWorkSpace()->GetModel(1)->GetBondList().at(i) << std::endl;
		}
		ofs_bondlist_inRbtGenome_ligand.close();

		buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/coordlist_inRbtGenome_ligand_proc_" + RbtWorkSpace::to_string(iMPIRank);
		std::ofstream ofs_coordlist_inRbtGenome_ligand(buf_str.c_str());
		for (int i = 0; i < pSF->GetWorkSpace()->GetModel(1)->GetAtomList().size(); i++) {
			ofs_coordlist_inRbtGenome_ligand << "ID" << i + 1 << ", " << pSF->GetWorkSpace()->GetModel(1)->GetAtomList().at(i)->GetX() << ", "
				<< pSF->GetWorkSpace()->GetModel(1)->GetAtomList().at(i)->GetY() << ", "
				<< pSF->GetWorkSpace()->GetModel(1)->GetAtomList().at(i)->GetZ() << std::endl;
		}
		ofs_coordlist_inRbtGenome_ligand.close();
#endif // output information of the receptor and ligand.



		// create a receptor PDB file
		int m_size_atomlist_receptor = 0;
		std::vector<std::string> m_v_sdf_receptor;
		int atom_id = 1;

		m_size_atomlist_receptor = pSF->GetWorkSpace()->GetModel(0)->GetAtomList().size();

		m_v_sdf_receptor.clear();
		
		for (int i = 0; i < m_size_atomlist_receptor; ++i) {
			if (pSF->GetWorkSpace()->GetModel(0)->GetAtomList()[i]->GetAtomName() == "H") {}
			else {

				buf_stringstream.str("");
				buf_stringstream.clear(std::stringstream::goodbit);
				buf_stringstream.setf(ios_base::left, ios_base::adjustfield);

				buf_stringstream << std::setw(6) << "ATOM";

				buf_stringstream.setf(ios_base::right, ios_base::adjustfield);
				buf_stringstream << std::setw(5) << atom_id;
				atom_id++;

				if (pSF->GetWorkSpace()->GetModel(0)->GetAtomList()[i]->GetAtomName().size() != 4) {
					buf_stringstream.setf(ios_base::left, ios_base::adjustfield);
					buf_stringstream << std::setw(2) << " "
						<< std::setw(3) << pSF->GetWorkSpace()->GetModel(0)->GetAtomList()[i]->GetAtomName();
				}
				else {
					buf_stringstream.setf(ios_base::left, ios_base::adjustfield);
					buf_stringstream << std::setw(1) << " "
						<< std::setw(4) << pSF->GetWorkSpace()->GetModel(0)->GetAtomList()[i]->GetAtomName();
				}

				buf_stringstream.setf(ios_base::right, ios_base::adjustfield);
				buf_stringstream << std::setw(4) << pSF->GetWorkSpace()->GetModel(0)->GetAtomList()[i]->GetSubunitName()
					<< std::setw(2) << " "
					<< std::setw(4) << pSF->GetWorkSpace()->GetModel(0)->GetAtomList()[i]->GetSubunitId()
					<< std::setw(4) << " ";

				buf_stringstream.precision(3);
				buf_stringstream.setf(ios_base::fixed, ios_base::floatfield);

				buf_stringstream << std::setw(8) << pSF->GetWorkSpace()->GetModel(0)->GetAtomList()[i]->GetX()
					<< std::setw(8) << pSF->GetWorkSpace()->GetModel(0)->GetAtomList()[i]->GetY()
					<< std::setw(8) << pSF->GetWorkSpace()->GetModel(0)->GetAtomList()[i]->GetZ();

				double buf_double_1 = 1;
				double buf_double_2 = 0;
				buf_stringstream.precision(2);
				buf_stringstream << std::setw(6) << buf_double_1
					<< std::setw(6) << buf_double_2
					<< std::setw(9) << " "
					//<< std::setw(3) << "" << ends;
					<< std::setw(3) << "";

				m_v_sdf_receptor.push_back(buf_stringstream.str());
			}
		}

		m_v_sdf_receptor.push_back("ENDMDL");

		buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/temp_receptor_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".pdb";
		std::ofstream ofs_PDB_file(buf_str.c_str());
		for (int i = 0; i < m_v_sdf_receptor.size(); ++i) {
			ofs_PDB_file << m_v_sdf_receptor.at(i) << std::endl;
		}
		ofs_PDB_file.close();



		RbtInt ri = pSF->GetWorkSpace()->GetModel(0)->GetCurrentCoords();
		pSF->GetWorkSpace()->GetModel(1)->SetDataValue("RI", ri);

		// create a ligand SDF file
		std::string atom_type_ligand;
		std::string bonded_atom_1;
		std::string bonded_atom_2;
		std::string bond_order;
		int m_size_atomlist_ligand = 0;
		int m_size_bondlist_ligand = 0;
		std::vector<std::string> m_v_sdf_ligand;

		m_size_atomlist_ligand = pSF->GetWorkSpace()->GetModel(1)->GetAtomList().size();
		m_size_bondlist_ligand = pSF->GetWorkSpace()->GetModel(1)->GetBondList().size();

		m_v_sdf_ligand.clear();
		m_v_sdf_ligand.resize(3 + 1 + m_size_atomlist_ligand + m_size_bondlist_ligand + 2);

		m_v_sdf_ligand.at(0) = "LIGAND";
		m_v_sdf_ligand.at(1) = "  extracted from rDock";
		m_v_sdf_ligand.at(2) = "";

		buf_stringstream.str("");
		buf_stringstream.clear(std::stringstream::goodbit);

		buf_stringstream.setf(ios_base::right, ios_base::adjustfield);

		buf_stringstream << std::setw(3) << m_size_atomlist_ligand << std::setw(3) << m_size_bondlist_ligand << "  0  0  0  0";
		m_v_sdf_ligand.at(3) = buf_stringstream.str();

		std::string x, y, z;

		for (int i = 0; i < m_size_atomlist_ligand; ++i) {

			int ii = i + 4;

			buf_stringstream.str("");
			buf_stringstream.clear(std::stringstream::goodbit);

			buf_stringstream.precision(4);
			buf_stringstream.setf(ios_base::fixed, ios_base::floatfield);

			buf_stringstream << pSF->GetWorkSpace()->GetModel(1)->GetAtomList()[i]->GetX();
			x = buf_stringstream.str();

			buf_stringstream.str("");
			buf_stringstream.clear(std::stringstream::goodbit);

			buf_stringstream.precision(4);
			buf_stringstream.setf(ios_base::fixed, ios_base::floatfield);

			buf_stringstream << pSF->GetWorkSpace()->GetModel(1)->GetAtomList()[i]->GetY();
			y = buf_stringstream.str();

			buf_stringstream.str("");
			buf_stringstream.clear(std::stringstream::goodbit);

			buf_stringstream.precision(4);
			buf_stringstream.setf(ios_base::fixed, ios_base::floatfield);

			buf_stringstream << pSF->GetWorkSpace()->GetModel(1)->GetAtomList()[i]->GetZ();
			z = buf_stringstream.str();

			atom_type_ligand = pSF->GetWorkSpace()->GetModel(1)->GetAtomList()[i]->GetAtomName();

			buf_str = atom_type_ligand[1];
			if (buf_str == "1" ||
				buf_str == "2" ||
				buf_str == "3" ||
				buf_str == "4" ||
				buf_str == "5" ||
				buf_str == "6" ||
				buf_str == "7" ||
				buf_str == "8" ||
				buf_str == "9") {
				atom_type_ligand = atom_type_ligand[0];
			}
			else {
				buf_stringstream.str("");
				buf_stringstream.clear(std::stringstream::goodbit);
				buf_stringstream << atom_type_ligand[0] << atom_type_ligand[1];
				atom_type_ligand = buf_stringstream.str();
			}

			buf_stringstream.str("");
			buf_stringstream.clear(std::stringstream::goodbit);

			buf_stringstream << std::right << std::setw(10) << x
				<< std::right << std::setw(10) << y
				<< std::right << std::setw(10) << z
				<< std::internal << std::setw(3) << atom_type_ligand
				<< std::right << std::setw(3) << "0"
				<< std::right << std::setw(3) << "0"
				<< std::right << std::setw(3) << "0";

			m_v_sdf_ligand.at(ii) = buf_stringstream.str();
		}

		for (int i = 0; i < m_size_bondlist_ligand; ++i) {

			int ii = i + 4 + m_size_atomlist_ligand;

			bonded_atom_1 = RbtWorkSpace::to_string(pSF->GetWorkSpace()->GetModel(1)->GetBondList()[i]->GetAtom1Ptr()->GetAtomId());
			bonded_atom_2 = RbtWorkSpace::to_string(pSF->GetWorkSpace()->GetModel(1)->GetBondList()[i]->GetAtom2Ptr()->GetAtomId());
			bond_order = RbtWorkSpace::to_string(pSF->GetWorkSpace()->GetModel(1)->GetBondList()[i]->GetFormalBondOrder());

			buf_stringstream.str("");
			buf_stringstream.clear(std::stringstream::goodbit);

			buf_stringstream << std::right << std::setw(3) << bonded_atom_1
				<< std::right << std::setw(3) << bonded_atom_2
				<< std::right << std::setw(3) << bond_order
				<< std::right << std::setw(3) << "0"
				<< std::right << std::setw(3) << "0"
				<< std::right << std::setw(3) << "0"
				<< std::right << std::setw(3) << "0";

			m_v_sdf_ligand.at(ii) = buf_stringstream.str();
		}

		m_v_sdf_ligand.at(4 + m_size_atomlist_ligand + m_size_bondlist_ligand) = "M  END";
		m_v_sdf_ligand.at(4 + m_size_atomlist_ligand + m_size_bondlist_ligand + 1) = "$$$$";

		buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/1_ligand_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".sd";
		std::ofstream ofs_SDF_file(buf_str.c_str());
		for (int i = 0; i < m_v_sdf_ligand.size(); ++i) {
			ofs_SDF_file << m_v_sdf_ligand.at(i);
			ofs_SDF_file << std::endl;
		}
		ofs_SDF_file.close();



		buf_str = "source ./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/shell_script_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".sh";
		system(buf_str.c_str());

		bool b_ligand_energy(true);
		//double ligand_energy;
		std::string ligand_energy;
		energy = "1e+5";
		counter = 0;
		buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/FINAL_RESULTS_MMPBSA_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".dat";
		std::ifstream ifs_ligand_energy(buf_str.c_str());
		while (ifs_ligand_energy && std::getline(ifs_ligand_energy, buf_str)) {

			buf_stringstream.str("");
			buf_stringstream.clear(std::stringstream::goodbit);

			buf_stringstream << buf_str;
			buf_stringstream >> buf_str;

			if (buf_str == "TOTAL") {
				counter++;
				if (counter == 3) {
					buf_stringstream >> buf_str;
					//std::cout << "ligand energy = " << buf_str << std::endl;
					//ligand_energy = std::atof(buf_str.c_str());
					ligand_energy = buf_str;
					//if (ligand_energy > std::atof(energy.c_str())) {
					if (std::atof(ligand_energy.c_str()) > std::atof(energy.c_str())) {
						b_ligand_energy = false;
					}
					break;
				}

			}
		}
		ifs_ligand_energy.close();

		buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/FINAL_RESULTS_MMPBSA_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".dat";
		std::ifstream ifs_energy(buf_str.c_str());
		if (b_ligand_energy) {
			while (ifs_energy && std::getline(ifs_energy, buf_str)) {

				buf_stringstream.str("");
				buf_stringstream.clear(std::stringstream::goodbit);

				buf_stringstream << buf_str;
				buf_stringstream >> buf_str;

				if (buf_str == "DELTA") {
					buf_stringstream >> buf_str;
					if (buf_str == "TOTAL") {
						buf_stringstream >> energy;
						//std::cout << "energy = " << energy << std::endl;
						break;
					}
				}
			}
		}
		else {
			//energy = RbtWorkSpace::to_string(ligand_energy);
			energy = ligand_energy;
		}
		ifs_energy.close();

#if 1
		buf_str = "./dir_proc_" + RbtWorkSpace::to_string(iMPIRank) + "/FINAL_RESULTS_MMPBSA_proc_" + RbtWorkSpace::to_string(iMPIRank) + ".dat";
		std::remove(buf_str.c_str());
#endif

#if 0 // debug write
		std::cout << "after execution of MM-GB/PBSA" << std::endl;
		std::cout << "type any key to continue" << std::endl;
		std::cin >> buf_str;
#endif // temp



		// put the MM-GB/PBSA energy into "m_score" as the docking score.
		double MMGBPBSA_energy = 0;
		MMGBPBSA_energy = std::atof(energy.c_str());
		MMGBPBSA_energy = -MMGBPBSA_energy;

		// add the cavity restraint term onto the score.
		double score_cavity_restraint_Ad = 0;
		counter_SF = 0;
		score_cavity_restraint_Ad = -pSF->Score();
		m_score = MMGBPBSA_energy + score_cavity_restraint_Ad;

#if 1 // debug write
		buf_str = "temp_score_proc_" + RbtWorkSpace::to_string(iMPIRank);
		std::ofstream ofs_temp_score(buf_str.c_str(), std::ios::app);
		static int static_int = 0;
		if (static_int == 0) {
			ofs_temp_score << "number of calculations, energy in total, MM-GB/PBSA energy, cavity restraint term, ligand energy" << std::endl;
			static_int = 1;
		}
		ofs_temp_score << counter_mmgbpbsa_calculations << ", " << m_score << ", " << MMGBPBSA_energy << ", " << score_cavity_restraint_Ad << ", " << ligand_energy << std::endl;
		ofs_temp_score.close();
		time_end_mmgbpbsa = time(NULL);
		//std::cout << "duration of the mmgbpbsa calculation = "
			//<< time_end_mmgbpbsa - time_start_mmgbpbsa << " sec" << std::endl;
#endif

  }
	else {
		m_score = -pSF->Score();
	}
#endif // AdvanceSoft (2018)

  }
  else {
    m_score = 0.0;
  }
  SetRWFitness(0.0, 0.0);
}

RbtDouble RbtGenome::SetRWFitness(RbtDouble sigmaOffset, RbtDouble partialSum) {
    //Apply sigma truncation to the raw score
    m_RWFitness = std::max(0.0, GetScore()-sigmaOffset);
    //The fitness value we store is the partial sum from previous RbtGenome elements
    //Not normalised at this point.
    m_RWFitness += partialSum;
    return m_RWFitness;
}

void RbtGenome::NormaliseRWFitness(RbtDouble total) {
    if (total > 0.0) {
        m_RWFitness /= total;
    }
}

void RbtGenome::Print(ostream& s) const {
  s << *m_chrom << endl;
  s << "Score: " << GetScore() << "; RWFitness: " << GetRWFitness() << endl;
}

ostream& operator<<(ostream& s, const RbtGenome& g) {
  g.Print(s);
  return s;
}


