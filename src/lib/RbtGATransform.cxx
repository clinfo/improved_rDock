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

#include "RbtGATransform.h"
#include "RbtWorkSpace.h"
#include "RbtPopulation.h"
#include "RbtSFRequest.h"
#include <iomanip>
using std::setw;

#if 1 // AdvanceSoft (2018)
#include <mpi.h>
#endif // AdvanceSoft(2018)

RbtString RbtGATransform::_CT("RbtGATransform");
RbtString RbtGATransform::_NEW_FRACTION("NEW_FRACTION");
RbtString RbtGATransform::_PCROSSOVER("PCROSSOVER");
RbtString RbtGATransform::_XOVERMUT("XOVERMUT");
RbtString RbtGATransform::_CMUTATE("CMUTATE");
RbtString RbtGATransform::_STEP_SIZE("STEP_SIZE");
RbtString RbtGATransform::_EQUALITY_THRESHOLD("EQUALITY_THRESHOLD");
RbtString RbtGATransform::_NCYCLES("NCYCLES");
RbtString RbtGATransform::_NCONVERGENCE("NCONVERGENCE");
RbtString RbtGATransform::_HISTORY_FREQ("HISTORY_FREQ");

#if 1 // AdvanceSoft (2018)
RbtString RbtGATransform::_NUM_CHROM_MIGRATION("NUM_CHROM_MIGRATION");
RbtString RbtGATransform::_CYCLE_MIGRATION("CYCLE_MIGRATION");
#endif // AdvanceSoft(2018)

#if 1 // AdvanceSoft (2018)
int counter_iCycle = 0;
bool b_islandGA_w_migration(false);
bool initial_b_islandGA_w_migration(false);
#endif // AdvanceSoft(2018)

RbtGATransform::RbtGATransform(const RbtString& strName) :
                               RbtBaseBiMolTransform(_CT,strName),
                               m_rand(Rbt::GetRbtRand())
{
  AddParameter(_NEW_FRACTION, 0.5);
  AddParameter(_PCROSSOVER, 0.4);
  AddParameter(_XOVERMUT, true);
  AddParameter(_CMUTATE, false);
  AddParameter(_STEP_SIZE, 1.0);
  //AddParameter(_EQUALITY_THRESHOLD, 0.1); // AdvanceSoft (2018)
  //AddParameter(_NCYCLES, 100); // AdvanceSoft (2018)
  //AddParameter(_NCONVERGENCE, 6); // AdvanceSoft (2018)
  AddParameter(_HISTORY_FREQ, 0);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);

#if 1 // AdvanceSoft (2018)
  AddParameter(_EQUALITY_THRESHOLD, 0.1);
  AddParameter(_NCONVERGENCE, 6);
  AddParameter(_NCYCLES, 100);
  //AddParameter(_NCYCLES, 10);
  AddParameter(_NUM_CHROM_MIGRATION, 2);
  AddParameter(_CYCLE_MIGRATION, 5);
#endif // AdvanceSoft(2018)

}

RbtGATransform::~RbtGATransform() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtGATransform::SetupReceptor(){}

void RbtGATransform::SetupLigand(){}

void RbtGATransform::SetupTransform(){}

void RbtGATransform::Execute() {

#if 1 // debug write
	std::cout << std::endl;
	std::cout << "***************************************************" << std::endl;
	std::cout << "The program is running in RbtGATransform::Execute()" << std::endl;
	std::cout << "***************************************************" << std::endl;
	std::cout << std::endl;
#endif // debug write

#if 1 // AdvanceSoft (2018)
	std::string buf_str;
	std::stringstream buf_stringstream;
	int counter = 0;
	std::vector<double> buf_vec_double;

	RbtInt iMPIRank, nMPISize;
	MPI_Comm_rank(MPI_COMM_WORLD, &iMPIRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
	//MPI_Status status;

	int r_proc_Ad = 0;
	if (iMPIRank % m_iIsland != 0) {
		r_proc_Ad = iMPIRank - (iMPIRank % m_iIsland);
	}

	b_islandGA_w_migration = initial_b_islandGA_w_migration;
#endif // AdvanceSoft(2018)

	RbtWorkSpace* pWorkSpace = GetWorkSpace();
	if (pWorkSpace == NULL) {
		return;
	}
	RbtBaseSF* pSF = pWorkSpace->GetSF();
	if (pSF == NULL) {
		return;
	}
	RbtPopulationPtr pop = pWorkSpace->GetPopulation();
	if (pop.Null() || (pop->GetMaxSize() < 1)) {
		return;
	}

	//Remove any partitioning from the scoring function
	//Not appropriate for a GA
	pSF->HandleRequest(new RbtSFPartitionRequest(0.0));
	//This forces the population to rescore all the individuals in case
	//the scoring function has changed
	pop->SetSF(pSF);

	RbtDouble newFraction = GetParameter(_NEW_FRACTION);
	RbtDouble pcross = GetParameter(_PCROSSOVER);
	RbtBool xovermut = GetParameter(_XOVERMUT);
	RbtBool cmutate = GetParameter(_CMUTATE);
	RbtDouble relStepSize = GetParameter(_STEP_SIZE);
	RbtDouble equalityThreshold = GetParameter(_EQUALITY_THRESHOLD);
	RbtInt nCycles = GetParameter(_NCYCLES);
	RbtInt nConvergence = GetParameter(_NCONVERGENCE);
	RbtInt nHisFreq = GetParameter(_HISTORY_FREQ);

	RbtInt popsize = pop->GetMaxSize();
	RbtInt nrepl = newFraction * popsize;
	RbtBool bHistory = nHisFreq > 0;
	RbtInt iTrace = GetTrace();

	RbtDouble bestScore = pop->Best()->GetScore();
	//Number of consecutive cycles with no improvement in best score
	RbtInt iConvergence = 0;

#if 1 // AdvanceSoft (2018)
	counter_GA++;
	RbtInt num_chrom_migration = GetParameter(_NUM_CHROM_MIGRATION);
	RbtInt cycle_migration = GetParameter(_CYCLE_MIGRATION);
#endif // AdvanceSoft(2018)

#if 1 // debug write
	std::cout << "newFraction = " << newFraction << std::endl;
	std::cout << "pcross = " << pcross << std::endl;
	std::cout << "xovermut = " << xovermut << std::endl;
	std::cout << "cmutate = " << cmutate << std::endl;
	std::cout << "relStepSize = " << relStepSize << std::endl;
	std::cout << "equalityThreshold = " << equalityThreshold << std::endl;
	std::cout << "nCycles = " << nCycles << std::endl;
	std::cout << "nConvergence = " << nConvergence << std::endl;
	std::cout << "nHisFreq = " << nHisFreq << std::endl;
	std::cout << "popsize = " << popsize << std::endl;
	std::cout << "nrepl = " << nrepl << std::endl;
	std::cout << "bHistory = " << bHistory << std::endl;
	std::cout << "iTrace = " << iTrace << std::endl;
	if (b_islandGA_w_migration) {
		std::cout << "num_chrom_migration = " << num_chrom_migration << std::endl;
		std::cout << "cycle_migration = " << cycle_migration << std::endl;
	}
#endif 

	if (iTrace > 0) {
		cout.precision(3);
		cout.setf(ios_base::fixed, ios_base::floatfield);
		cout.setf(ios_base::right, ios_base::adjustfield);
		cout << endl
			<< setw(5) << "CYCLE"
			<< setw(5) << "CONV"
			<< setw(10) << "BEST"
			<< setw(10) << "MEAN"
			<< setw(10) << "VAR"
			<< endl;

		cout << endl
			<< setw(5) << "Init"
			<< setw(5) << "-"
			<< setw(10) << bestScore
			<< setw(10) << pop->GetScoreMean()
			<< setw(10) << pop->GetScoreVariance()
			<< endl;
	}

	for (RbtInt iCycle = 0;
		(iCycle < nCycles) && (iConvergence < nConvergence);
		++iCycle) {
		if (bHistory && ((iCycle % nHisFreq) == 0)) {
			pop->Best()->GetChrom()->SyncToModel();
			pWorkSpace->SaveHistory(true);
		}

		pop->GAstep(nrepl, relStepSize, equalityThreshold, pcross, xovermut, cmutate);

		RbtDouble score = pop->Best()->GetScore();
		if (score > bestScore) {
			bestScore = score;
			iConvergence = 0;
		}
		else {
			iConvergence++;
		}
		if (iTrace > 0) {
			cout << setw(5) << iCycle
				<< setw(5) << iConvergence
				<< setw(10) << score
				<< setw(10) << pop->GetScoreMean()
				<< setw(10) << pop->GetScoreVariance()
				<< endl;
		}



#if 1 // AdvanceSoft (2018)
		if (b_islandGA_w_migration) {
			if (m_iIsland != 1) {

				RbtGenomeList& m_pop_Ad = pop->GetGenomeList_Ad();
				int chrom_length = m_pop_Ad[0]->GetChrom()->GetLength();

				std::vector<int> migration_direction;
				migration_direction.resize(m_iIsland);

				for (int i = 0; i < m_iIsland; ++i) {
					migration_direction[i] = (i + 1) % m_iIsland;
				}

				int destination = migration_direction[iMPIRank % m_iIsland];
				if (iMPIRank >= m_iIsland) {
					destination += m_iIsland * (iMPIRank / m_iIsland);
				}
				int source = 0;
				for (int i = 0; i < m_iIsland; ++i) {
					if (migration_direction[i] == (iMPIRank % m_iIsland)) {
						source = i;
						if (iMPIRank >= m_iIsland) {
							source += m_iIsland * (iMPIRank / m_iIsland);
						}
					}
				}

				std::vector<int> v_b_islandGA_w_migration;
				v_b_islandGA_w_migration.resize(m_iIsland);

				if (iConvergence == nConvergence) {
					b_islandGA_w_migration = 0;
				}
				else {
					b_islandGA_w_migration = 1;
				}

				if ((iMPIRank % m_iIsland) == 0) {
					v_b_islandGA_w_migration[0] = b_islandGA_w_migration;
				}

				if (iMPIRank % m_iIsland != 0) {
					MPI_Send(&b_islandGA_w_migration, 1, MPI_INT, r_proc_Ad, iMPIRank, MPI_COMM_WORLD);
				}
				else {
					for (int i = 1; i < m_iIsland; ++i) {
						MPI_Status status_v_b_islandGA_w_migration;
						MPI_Recv(&v_b_islandGA_w_migration[i], 1, MPI_INT, iMPIRank + i, iMPIRank + i, MPI_COMM_WORLD, &status_v_b_islandGA_w_migration);
					}
				}

				if ((iMPIRank % m_iIsland) == 0) {
					for (int i = 1; i < m_iIsland; ++i) {
						MPI_Send(&v_b_islandGA_w_migration[0], m_iIsland, MPI_INT, iMPIRank + i, iMPIRank, MPI_COMM_WORLD);
					}
				}
				else {
					MPI_Status status_v_b_islandGA_w_migration_2;
					MPI_Recv(&v_b_islandGA_w_migration[0], m_iIsland, MPI_INT, r_proc_Ad, r_proc_Ad, MPI_COMM_WORLD, &status_v_b_islandGA_w_migration_2);
				}

				for (int i = 0; i < m_iIsland; ++i) {
					if (v_b_islandGA_w_migration[i] == 0) {
						b_islandGA_w_migration = 0;
					}
				}

				// migration starts hereafter
				if (!b_islandGA_w_migration) {
				}
				else {
					if ((iCycle % cycle_migration) == 0) {

						std::vector<std::vector<double> > v_send_chrom;
						v_send_chrom.resize(num_chrom_migration);
						for (int i = 0; i < num_chrom_migration; ++i) {
							v_send_chrom[i].resize(chrom_length);
						}
						std::vector<std::vector<double> > v_receive_chrom;
						v_receive_chrom.resize(num_chrom_migration);
						for (int i = 0; i < num_chrom_migration; ++i) {
							v_receive_chrom[i].resize(chrom_length);
						}
						RbtChromElementList send_chrom_element_list;
						send_chrom_element_list.resize(num_chrom_migration);
						RbtChromElementList receive_chrom_element_list;
						receive_chrom_element_list.resize(num_chrom_migration);

						for (int i = 0; i < num_chrom_migration; ++i) {
							receive_chrom_element_list[i] = m_pop_Ad[(m_pop_Ad.size() - i) - 1]->GetChrom();
						}

						for (int i = 0; i < num_chrom_migration; ++i) {
							send_chrom_element_list[i] = m_pop_Ad[i]->GetChrom();
							buf_vec_double.clear();
							send_chrom_element_list[i]->GetVector(buf_vec_double);
							v_send_chrom[i] = buf_vec_double;
						}

						for (int i = 0; i < num_chrom_migration; ++i) {
							MPI_Send(&v_send_chrom[i][0], chrom_length, MPI_DOUBLE, destination, iMPIRank, MPI_COMM_WORLD);
						}
						for (int i = 0; i < num_chrom_migration; ++i) {
							MPI_Status status_receive_chrom;
							MPI_Recv(&v_receive_chrom[i][0], chrom_length, MPI_DOUBLE, source, source, MPI_COMM_WORLD, &status_receive_chrom);
						}

						for (int i = 0; i < num_chrom_migration; ++i) {
							buf_vec_double.clear();
							buf_vec_double = v_receive_chrom[i];
							receive_chrom_element_list[i]->SetVector(buf_vec_double);
						}

						pop->SetSF(pSF);
						std::stable_sort(m_pop_Ad.begin(), m_pop_Ad.end(), Rbt::GenomeCmp_Score());
					} // if ((iCycle % cycle_migration) == 0)
				} // if (!b_islandGA_w_migration)
			} // if (m_iIsland != 1)
		} // if (b_islandGA_w_migration)

		counter_iCycle++;
	} // for (RbtInt iCycle = 0; (iCycle < nCycles) && (iConvergence < nConvergence); ++iCycle)

	if (m_iIsland != 1) {

		double score_Ad = pop->Best()->GetScore(); // get the best score in the island
		std::vector<double> v_score_Ad(m_iIsland); // vector of scores from each island

		if ((iMPIRank % m_iIsland) == 0) {
			v_score_Ad[0] = score_Ad;
		}

		if (iMPIRank % m_iIsland != 0) {
			MPI_Send(&score_Ad, 1, MPI_DOUBLE, r_proc_Ad, iMPIRank, MPI_COMM_WORLD);
		}
		else {
			for (int i = 1; i < m_iIsland; ++i) {
				MPI_Status status_score_Ad;
				MPI_Recv(&v_score_Ad[i], 1, MPI_DOUBLE, iMPIRank + i, iMPIRank + i, MPI_COMM_WORLD, &status_score_Ad);
			}
		}

		vector<double> v_best_chrom_Ad;
		v_best_chrom_Ad.clear();
		pop->Best()->GetChrom()->GetVector(v_best_chrom_Ad);
		int chrom_length = v_best_chrom_Ad.size();

		std::vector< std::vector<double> > vv_best_chrom_Ad;
		vv_best_chrom_Ad.resize(m_iIsland);
		for (int i = 0; i < m_iIsland; ++i) {
			vv_best_chrom_Ad[i].resize(chrom_length);
		}

		if ((iMPIRank % m_iIsland) == 0) {
			vv_best_chrom_Ad[0] = v_best_chrom_Ad;
		}

		if (iMPIRank % m_iIsland != 0) {
			MPI_Send(&v_best_chrom_Ad[0], chrom_length, MPI_DOUBLE, r_proc_Ad, iMPIRank, MPI_COMM_WORLD);
		}
		else {
			for (int i = 1; i < m_iIsland; ++i) {
				MPI_Status status_best_chrom_Ad;
				MPI_Recv(&vv_best_chrom_Ad[i][0], chrom_length, MPI_DOUBLE, iMPIRank + i, iMPIRank + i, MPI_COMM_WORLD, &status_best_chrom_Ad);
			}
		}

		vector<double> v_best_chrom_in_all_islands_Ad;
		v_best_chrom_in_all_islands_Ad.resize(chrom_length);

		if ((iMPIRank % m_iIsland) == 0) {

			int i_best_island_Ad = 0;
			int max_score_Ad = v_score_Ad[0];

			for (int i = 1; i < m_iIsland; ++i) {
				if (v_score_Ad[i] > max_score_Ad) {
					max_score_Ad = v_score_Ad[i];
					i_best_island_Ad = i;
				}
			}

			v_best_chrom_in_all_islands_Ad = vv_best_chrom_Ad[i_best_island_Ad];

			for (int i = 1; i < m_iIsland; ++i) {
				MPI_Send(&v_best_chrom_in_all_islands_Ad[0], chrom_length, MPI_DOUBLE, iMPIRank + i, iMPIRank, MPI_COMM_WORLD);
			}
		}
		else {
			MPI_Status status_best_chrom_in_all_islands_Ad;
			MPI_Recv(&v_best_chrom_in_all_islands_Ad[0], chrom_length, MPI_DOUBLE, r_proc_Ad, r_proc_Ad, MPI_COMM_WORLD, &status_best_chrom_in_all_islands_Ad);
		}

		pop->Best()->GetChrom()->SetVector(v_best_chrom_in_all_islands_Ad);
		pop->Best()->GetChrom()->SyncToModel();
		RbtInt ri = GetReceptor()->GetCurrentCoords();
		GetLigand()->SetDataValue("RI", ri);
	}
	else if (m_iIsland == 1) {

#if 1 // original code
		pop->Best()->GetChrom()->SyncToModel();
		RbtInt ri = GetReceptor()->GetCurrentCoords();
		GetLigand()->SetDataValue("RI", ri);
#endif // original code

	}
#endif // AdvanceSoft(2018)
}
