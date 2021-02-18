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

#include "RbtPopulation.h"
#include "RbtDebug.h"
#include "RbtDockingError.h"
#include <algorithm>

#if 1 // AdvanceSoft (2018)
#include <mpi.h>
#include <fstream>
#include <sstream>
bool b_parallel_mmgbpbsa(false);
#endif // AdvanceSoft (2018)

RbtString RbtPopulation::_CT("RbtPopulation"); 

RbtPopulation::RbtPopulation(RbtChromElement* pChr, RbtInt size, RbtBaseSF* pSF)
    throw (RbtError)
        : m_size(size), m_c(2.0), m_pSF(pSF), m_rand(Rbt::GetRbtRand()),
        m_scoreMean(0.0), m_scoreVariance(0.0)
{
  if (pChr == NULL) {
    throw RbtBadArgument(_WHERE_, "Null chromosome element passed to RbtPopulation constructor");
  }
  else if (size <= 0) {
    throw RbtBadArgument(_WHERE_, "Population size must be positive (non-zero)");
  }
  //Create a random population
  m_pop.reserve(m_size);
  for (RbtInt i = 0; i < m_size; ++i) {
    //The RbtGenome constructor clones the chromosome to create an independent copy
    RbtGenomePtr genome = new RbtGenome(pChr);
    genome->GetChrom()->Randomise();
    m_pop.push_back(genome);
  }
  //Calculate the scores and evaluate roulette wheel fitness
  SetSF(m_pSF);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtPopulation::~RbtPopulation() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

#if 0 // original code
//Sets the scoring function used for ranking genomes
void RbtPopulation::SetSF(RbtBaseSF* pSF) throw (RbtError) {
  if (pSF == NULL) {
    throw RbtBadArgument(_WHERE_, "Null scoring function passed to SetSF");
  }
  m_pSF = pSF;
  for (RbtGenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    (*iter)->SetScore(m_pSF);
  }
  std::stable_sort(m_pop.begin(), m_pop.end(), Rbt::GenomeCmp_Score());
  EvaluateRWFitness();
}
#endif // original code

#if 1 // AdvanceSoft (2018)
//Sets the scoring function used for ranking genomes
void RbtPopulation::SetSF(RbtBaseSF* pSF) throw (RbtError) {

	if (pSF == NULL) {
		throw RbtBadArgument(_WHERE_, "Null scoring function passed to SetSF");
	}

	std::string buf_str;
	std::stringstream buf_stringstream;
	int counter = 0;
	RbtInt iMPIRank, nMPISize;
	MPI_Comm_rank(MPI_COMM_WORLD, &iMPIRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);

#if 0 // degub write
	buf_stringstream.str("");
	buf_stringstream.clear(std::stringstream::goodbit);
	buf_stringstream << "out_1_initial_pop_proc_" << iMPIRank << ".dat";
	buf_str = buf_stringstream.str();
	std::ofstream ofs_initial_pop_proc(buf_str.c_str());
	ofs_initial_pop_proc << "m_pop.size() = " << m_pop.size() << std::endl;
	for (int i = 0; i < m_pop.size(); ++i) {
		ofs_initial_pop_proc << *m_pop[i] << std::endl;
	}
	ofs_initial_pop_proc.close();
#endif // debug write

	std::vector<double> v_score_Ad(m_pop.size(), 0); // all scores will be put in
	int num_chrom_per_process = m_pop.size() / nMPISize;
	std::vector<double> v_score_Ad_process(num_chrom_per_process);
	int num_residue_chrom = m_pop.size() % nMPISize;
	double score_Ad_residue = 0;
	std::vector<double> v_score_Ad_residue(nMPISize, 0);

#if 0 // temp
	if (b_parallel_mmgbpbsa) {
		for (int i = 0; i < m_pop.size(); ++i) {
			m_pop[i]->SetScore_Ad(0);
			m_pop[i]->SetRWFitness_Ad(0);
		}
	}
#endif // temp

	counter = 0;
	int counter_residue = 0;
	m_pSF = pSF;
	for (RbtGenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {

		if (nMPISize == 1) {
			(*iter)->SetScore(m_pSF);
		}
		else if (!b_parallel_mmgbpbsa) {
			(*iter)->SetScore(m_pSF);
		}
		else {
			if (counter < (iMPIRank + 1) * num_chrom_per_process) {
				if (counter >= iMPIRank * num_chrom_per_process) {
					(*iter)->SetScore(m_pSF);
					v_score_Ad_process[counter % num_chrom_per_process] = (*iter)->GetScore();
				}
			}
			if (counter >= nMPISize * num_chrom_per_process) {
				if (counter_residue == iMPIRank) {
					(*iter)->SetScore(m_pSF);
					score_Ad_residue = (*iter)->GetScore();
				}
				counter_residue++;
			}
			counter++;
		}
	}

#if 0 // degub write
	buf_stringstream.str("");
	buf_stringstream.clear(std::stringstream::goodbit);
	buf_stringstream << "out_2_pop_after_SetScore_proc_" << iMPIRank << ".dat";
	buf_str = buf_stringstream.str();
	std::ofstream ofs_pop_after_SetScore(buf_str.c_str());
	ofs_pop_after_SetScore << "m_pop.size() = " << m_pop.size() << std::endl;
	for (int i = 0; i < m_pop.size(); ++i) {
		ofs_pop_after_SetScore << "pop No.: " << i << std::endl;
		ofs_pop_after_SetScore << *m_pop[i] << std::endl;
	}
	ofs_pop_after_SetScore.close();
#endif // debug write

	if (nMPISize == 1) {}
	else if (!b_parallel_mmgbpbsa) {}
	else {
		MPI_Gather(&v_score_Ad_process[0], num_chrom_per_process, MPI_DOUBLE, &v_score_Ad[0], num_chrom_per_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#if 0 // degub write
		buf_stringstream.str("");
		buf_stringstream.clear(std::stringstream::goodbit);
		buf_stringstream << "out_v_score_Ad_process_" << iMPIRank << ".dat";
		buf_str = buf_stringstream.str();
		std::ofstream ofs_v_score_Ad_process(buf_str.c_str());
		ofs_v_score_Ad_process << "v_score_Ad_process.size() = " << v_score_Ad_process.size() << std::endl;
		for (int i = 0; i < v_score_Ad_process.size(); ++i) {
			ofs_v_score_Ad_process << v_score_Ad_process[i] << std::endl;
		}
		ofs_v_score_Ad_process.close();
#endif // debug write

#if 0 // degub write
		buf_stringstream.str("");
		buf_stringstream.clear(std::stringstream::goodbit);
		buf_stringstream << "out_v_score_Ad_" << iMPIRank << ".dat";
		buf_str = buf_stringstream.str();
		std::ofstream ofs_v_score_Ad(buf_str.c_str());
		ofs_v_score_Ad << "v_score_Ad.size() = " << v_score_Ad.size() << std::endl;
		for (int i = 0; i < v_score_Ad.size(); ++i) {
			ofs_v_score_Ad << v_score_Ad[i] << std::endl;
		}
		ofs_v_score_Ad.close();
#endif // debug write

		MPI_Gather(&score_Ad_residue, 1, MPI_DOUBLE, &v_score_Ad_residue[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (iMPIRank == 0) {
			for (int i = 0; i < counter_residue; ++i) {
				v_score_Ad[m_pop.size() - (counter_residue - i)] = v_score_Ad_residue[i];
			}
		}

#if 0 // degub write
		buf_stringstream.str("");
		buf_stringstream.clear(std::stringstream::goodbit);
		buf_stringstream << "out_score_Ad_residue_" << iMPIRank << ".dat";
		buf_str = buf_stringstream.str();
		std::ofstream ofs_score_Ad_residue(buf_str.c_str());
		ofs_score_Ad_residue << "score_Ad_residue = " << score_Ad_residue << std::endl;
		ofs_score_Ad_residue.close();
#endif // debug write

#if 0 // degub write
		buf_stringstream.str("");
		buf_stringstream.clear(std::stringstream::goodbit);
		buf_stringstream << "out_v_score_Ad_residue_" << iMPIRank << ".dat";
		buf_str = buf_stringstream.str();
		std::ofstream ofs_v_score_Ad_residue(buf_str.c_str());
		ofs_v_score_Ad_residue << "v_score_Ad_residue.size() = " << v_score_Ad_residue.size() << std::endl;
		for (int i = 0; i < v_score_Ad_residue.size(); ++i) {
			ofs_v_score_Ad_residue << v_score_Ad_residue[i] << std::endl;
		}
		ofs_v_score_Ad_residue.close();
#endif // debug write

		MPI_Bcast(&v_score_Ad[0], v_score_Ad.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

		for (int i = 0; i < m_pop.size(); ++i) {
			m_pop[i]->SetScore_Ad(v_score_Ad[i]);
		}

#if 0 // debug write
		static int static_int = 0;
		if (static_int == 0) {
			if (iMPIRank == 0) {
				std::ofstream ofs_v_score_Ad("out_v_score_Ad.dat", std::ios_base::out);
				for (int i = 0; i < v_score_Ad.size(); ++i) {
					ofs_v_score_Ad << v_score_Ad[i] << std::endl;
				}
				ofs_v_score_Ad.close();
			}
			static_int = 1;
		}
		else {
			if (iMPIRank == 0) {
				std::ofstream ofs_v_score_Ad("out_v_score_Ad.dat", std::ios_base::app);
				for (int i = 0; i < v_score_Ad.size(); ++i) {
					ofs_v_score_Ad << v_score_Ad[i] << std::endl;
				}
				ofs_v_score_Ad.close();
			}
		}
#endif // debug write

	}


#if 0 // degub write
	buf_stringstream.str("");
	buf_stringstream.clear(std::stringstream::goodbit);
	buf_stringstream << "out_3_pop_after_communication_proc_" << iMPIRank << ".dat";
	buf_str = buf_stringstream.str();
	std::ofstream ofs_pop_after_communication(buf_str.c_str());
	ofs_pop_after_communication << "m_pop.size() = " << m_pop.size() << std::endl;
	for (int i = 0; i < m_pop.size(); ++i) {
		ofs_pop_after_communication << "pop No.: " << i << std::endl;
		ofs_pop_after_communication << *m_pop[i] << std::endl;
	}
	ofs_pop_after_communication.close();
#endif // debug write



	std::stable_sort(m_pop.begin(), m_pop.end(), Rbt::GenomeCmp_Score());



#if 0 // original code
	m_pSF = pSF;
	for (RbtGenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
		(*iter)->SetScore(m_pSF);
	}
	std::stable_sort(m_pop.begin(), m_pop.end(), Rbt::GenomeCmp_Score());
#endif // original code



	EvaluateRWFitness();



#if 0 // degub write
	buf_stringstream.str("");
	buf_stringstream.clear(std::stringstream::goodbit);
	buf_stringstream << "out_4_pop_after_sort_proc_" << iMPIRank << ".dat";
	buf_str = buf_stringstream.str();
	std::ofstream ofs_pop_after_sort(buf_str.c_str());
	ofs_pop_after_sort << "m_pop.size() = " << m_pop.size() << std::endl;
	for (int i = 0; i < m_pop.size(); ++i) {
		ofs_pop_after_sort << "pop No.: " << i << std::endl;
		ofs_pop_after_sort << *m_pop[i] << std::endl;
	}
	ofs_pop_after_sort.close();
#endif // debug write

}
#endif // AdvanceSoft (2018)


void RbtPopulation::GAstep(RbtInt nReplicates, RbtDouble relStepSize,
                            RbtDouble equalityThreshold, RbtDouble pcross,
                            RbtBool xovermut, RbtBool cmutate) throw (RbtError)
{
  if (nReplicates <= 0) {
    throw RbtBadArgument(_WHERE_, "nReplicates must be positive (non-zero)");
  }
  RbtGenomeList newPop;
  newPop.reserve(nReplicates);
  for (RbtInt i = 0 ; i < nReplicates / 2 ; i++) {
    RbtGenomePtr mother = RouletteWheelSelect();
    RbtGenomePtr father = RouletteWheelSelect();
    //Check that mother and father are not the same genome
    //The check is on the pointers, not that the chromosomes are near-equal
    //If we repeatedly get the same genomes selected, this must mean
    //the population lacks diversity
    RbtInt j = 0;
    while (father == mother) {
      father = RouletteWheelSelect();
      if (j > 100)  
        throw RbtDockingError(_WHERE_, 
				            "Population failure - not enough diversity");
      j++;
    }
    RbtGenomePtr child1 = new RbtGenome(*mother);
    RbtGenomePtr child2 = new RbtGenome(*father);
    //Crossover
    if (m_rand.GetRandom01() < pcross) {
      Rbt::Crossover(father->GetChrom(), mother->GetChrom(), 
                     child1->GetChrom(), child2->GetChrom());
      //Cauchy mutation following crossover
      if (xovermut) {
        child1->GetChrom()->CauchyMutate(0.0, relStepSize);
        child2->GetChrom()->CauchyMutate(0.0, relStepSize);
      }
    }
    //Mutation
    else {
      //Cauchy mutation
      if (cmutate) {
        child1->GetChrom()->CauchyMutate(0.0, relStepSize);
        child2->GetChrom()->CauchyMutate(0.0, relStepSize);
      }
      //Regular mutation
      else {
        child1->GetChrom()->Mutate(relStepSize);
        child2->GetChrom()->Mutate(relStepSize);
      }
    }
    newPop.push_back(child1);
    newPop.push_back(child2);
  }
  //check if one more is needed (odd nReplicates).
  if (nReplicates % 2) {
    RbtGenomePtr mother = RouletteWheelSelect();
    RbtGenomePtr child = new RbtGenome(*mother);
    child->GetChrom()->CauchyMutate(0.0, relStepSize);
    newPop.push_back(child);
  }
  MergeNewPop(newPop, equalityThreshold);
  EvaluateRWFitness();
}

RbtGenomePtr RbtPopulation::Best() const {
  return m_pop.empty() ? RbtGenomePtr() : m_pop.front();
}

void RbtPopulation::Print(ostream& s) const {
  s << m_pop.size() << endl;
  for (RbtGenomeListConstIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    (*iter)->Print(s);
  }
}

ostream& operator<<(ostream& s, const RbtPopulation& p) {
  p.Print(s);
  return s;
}


RbtGenomePtr RbtPopulation::RouletteWheelSelect() const {
  RbtDouble cutoff = m_rand.GetRandom01();
  RbtInt size = m_pop.size();
  RbtInt lower = 0;
  RbtInt upper = size - 1;
  while (upper >= lower) {
    RbtInt i = lower + (upper - lower) / 2;
    if (m_pop[i]->GetRWFitness() > cutoff)
      upper = i - 1;
    else
      lower = i + 1;
  } 
  //make sure lower is a number between 0 and size - 1
  lower = std::min(size-1, lower);
  lower = std::max(0, lower);
  return (m_pop[lower]);
}

#if 0 // original code
void RbtPopulation::MergeNewPop(RbtGenomeList& newPop, RbtDouble equalityThreshold) {
  //Assume newPop needs scoring and sorting
  for (RbtGenomeListIter iter = newPop.begin(); iter != newPop.end(); ++iter) {
    (*iter)->SetScore(m_pSF);
  }
  std::stable_sort(newPop.begin(), newPop.end(), Rbt::GenomeCmp_Score());
     
  RbtGenomeList mergedPop;
  mergedPop.reserve(m_pop.size() + newPop.size());
  //Merge pops by score
  std::merge(m_pop.begin(), m_pop.end(),
                newPop.begin(), newPop.end(), 
                std::back_inserter(mergedPop), Rbt::GenomeCmp_Score());
  //Remove neighbouring duplicates by equality of chromosome element values
  RbtGenomeListIter end = std::unique(mergedPop.begin(), mergedPop.end(),
                                        Rbt::isGenome_eq(equalityThreshold));
  mergedPop.erase(end, mergedPop.end());
  m_pop.clear();
  end = (mergedPop.size() > m_size) ? (mergedPop.begin() + m_size) : mergedPop.end();
  std::copy(mergedPop.begin(), end, back_inserter(m_pop));
}
#endif // original code

#if 1 // AdvanceSoft (2018)
void RbtPopulation::MergeNewPop(RbtGenomeList& newPop, RbtDouble equalityThreshold) {

	std::string buf_str;
	std::stringstream buf_stringstream;
	int counter = 0;
	RbtInt iMPIRank, nMPISize;
	MPI_Comm_rank(MPI_COMM_WORLD, &iMPIRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);

#if 0 // degub write
	buf_stringstream.str("");
	buf_stringstream.clear(std::stringstream::goodbit);
	buf_stringstream << "out_1_initial_pop_proc_" << iMPIRank << ".dat";
	buf_str = buf_stringstream.str();
	std::ofstream ofs_initial_pop_proc(buf_str.c_str());
	for (int i = 0; i < newPop.size(); ++i) {
		ofs_initial_pop_proc << *newPop[i] << std::endl;
	}
	ofs_initial_pop_proc.close();
#endif // debug write

	std::vector<double> v_score_Ad(newPop.size(), 0); // all scores will be put in
	int num_chrom_per_process = newPop.size() / nMPISize;
	std::vector<double> v_score_Ad_process(num_chrom_per_process);
	int num_residue_chrom = newPop.size() % nMPISize;
	double score_Ad_residue = 0;
	std::vector<double> v_score_Ad_residue(nMPISize, 0);

#if 0 // temp
	if (b_parallel_mmgbpbsa) {
		for (int i = 0; i < newPop.size(); ++i) {
			newPop[i]->SetScore_Ad(0);
			newPop[i]->SetRWFitness_Ad(0);
		}
	}
#endif // temp

	counter = 0;
	int counter_residue = 0;
	//m_pSF = pSF;
	for (RbtGenomeListIter iter = newPop.begin(); iter != newPop.end(); ++iter) {

		if (nMPISize == 1) {
			(*iter)->SetScore(m_pSF);
		}
		else if (!b_parallel_mmgbpbsa) {
			(*iter)->SetScore(m_pSF);
		}
		else {
			if (counter < (iMPIRank + 1) * num_chrom_per_process) {
				if (counter >= iMPIRank * num_chrom_per_process) {
					(*iter)->SetScore(m_pSF);
					v_score_Ad_process[counter % num_chrom_per_process] = (*iter)->GetScore();
				}
			}
			if (counter >= nMPISize * num_chrom_per_process) {
				if (counter_residue == iMPIRank) {
					(*iter)->SetScore(m_pSF);
					score_Ad_residue = (*iter)->GetScore();
				}
				counter_residue++;
			}
			counter++;
		}
	}

#if 0 // degub write
	buf_stringstream.str("");
	buf_stringstream.clear(std::stringstream::goodbit);
	buf_stringstream << "out_2_pop_after_SetScore_proc_" << iMPIRank << ".dat";
	buf_str = buf_stringstream.str();
	std::ofstream ofs_pop(buf_str.c_str());
	ofs_pop << "newPop.size() = " << newPop.size() << std::endl;
	for (int i = 0; i < newPop.size(); ++i) {
		ofs_pop << "pop No.: " << i << std::endl;
		ofs_pop << *newPop[i] << std::endl;
	}
	ofs_pop.close();
#endif // debug write

	if (nMPISize == 1) {}
	else if (!b_parallel_mmgbpbsa) {}
	else {
		MPI_Gather(&v_score_Ad_process[0], num_chrom_per_process, MPI_DOUBLE, &v_score_Ad[0], num_chrom_per_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#if 0 // degub write
		buf_stringstream.str("");
		buf_stringstream.clear(std::stringstream::goodbit);
		buf_stringstream << "out_v_score_Ad_process_" << iMPIRank << ".dat";
		buf_str = buf_stringstream.str();
		std::ofstream ofs_v_score_Ad_process(buf_str.c_str());
		ofs_v_score_Ad_process << "v_score_Ad_process.size() = " << v_score_Ad_process.size() << std::endl;
		for (int i = 0; i < v_score_Ad_process.size(); ++i) {
			ofs_v_score_Ad_process << v_score_Ad_process[i] << std::endl;
		}
		ofs_v_score_Ad_process.close();
#endif // debug write

#if 0 // degub write
		buf_stringstream.str("");
		buf_stringstream.clear(std::stringstream::goodbit);
		buf_stringstream << "out_v_score_Ad_" << iMPIRank << ".dat";
		buf_str = buf_stringstream.str();
		std::ofstream ofs_v_score_Ad(buf_str.c_str());
		ofs_v_score_Ad << "v_score_Ad.size() = " << v_score_Ad.size() << std::endl;
		for (int i = 0; i < v_score_Ad.size(); ++i) {
			ofs_v_score_Ad << v_score_Ad[i] << std::endl;
		}
		ofs_v_score_Ad.close();
#endif // debug write

		MPI_Gather(&score_Ad_residue, 1, MPI_DOUBLE, &v_score_Ad_residue[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (iMPIRank == 0) {
			for (int i = 0; i < counter_residue; ++i) {
				v_score_Ad[newPop.size() - (counter_residue - i)] = v_score_Ad_residue[i];
			}
		}

#if 0 // degub write
		buf_stringstream.str("");
		buf_stringstream.clear(std::stringstream::goodbit);
		buf_stringstream << "out_score_Ad_residue_" << iMPIRank << ".dat";
		buf_str = buf_stringstream.str();
		std::ofstream ofs_score_Ad_residue(buf_str.c_str());
		ofs_score_Ad_residue << "score_Ad_residue = " << score_Ad_residue << std::endl;
		ofs_score_Ad_residue.close();
#endif // debug write

#if 0 // degub write
		buf_stringstream.str("");
		buf_stringstream.clear(std::stringstream::goodbit);
		buf_stringstream << "out_v_score_Ad_residue_" << iMPIRank << ".dat";
		buf_str = buf_stringstream.str();
		std::ofstream ofs_v_score_Ad_residue(buf_str.c_str());
		ofs_v_score_Ad_residue << "v_score_Ad_residue.size() = " << v_score_Ad_residue.size() << std::endl;
		for (int i = 0; i < v_score_Ad_residue.size(); ++i) {
			ofs_v_score_Ad_residue << v_score_Ad_residue[i] << std::endl;
		}
		ofs_v_score_Ad_residue.close();
#endif // debug write

		MPI_Bcast(&v_score_Ad[0], v_score_Ad.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

		for (int i = 0; i < newPop.size(); ++i) {
			newPop[i]->SetScore_Ad(v_score_Ad[i]);
		}

#if 0 // debug write
		// output v_score_Ad
		if (iMPIRank == 0) {
			std::ofstream ofs_v_score_Ad("out_v_score_Ad.dat", std::ios_base::app);
			for (int i = 0; i < v_score_Ad.size(); ++i) {
				ofs_v_score_Ad << v_score_Ad[i] << std::endl;
			}
			ofs_v_score_Ad.close();
		}
#endif // debug write

	}

#if 0 // debug write
	buf_stringstream.str("");
	buf_stringstream.clear(std::stringstream::goodbit);
	buf_stringstream << "out_3_pop_after_communication_proc_" << iMPIRank << ".dat";
	buf_str = buf_stringstream.str();
	std::ofstream ofs_pop_2(buf_str.c_str());
	ofs_pop_2 << "newPop.size() = " << newPop.size() << std::endl;
	for (int i = 0; i < newPop.size(); ++i) {
		ofs_pop_2 << "pop No.: " << i << std::endl;
		ofs_pop_2 << *newPop[i] << std::endl;
	}
	ofs_pop_2.close();
#endif // debug write



	std::stable_sort(newPop.begin(), newPop.end(), Rbt::GenomeCmp_Score());

	RbtGenomeList mergedPop;
	mergedPop.reserve(m_pop.size() + newPop.size());
	//Merge pops by score
	std::merge(m_pop.begin(), m_pop.end(),
		newPop.begin(), newPop.end(),
		std::back_inserter(mergedPop), Rbt::GenomeCmp_Score());
	//Remove neighbouring duplicates by equality of chromosome element values
	RbtGenomeListIter end = std::unique(mergedPop.begin(), mergedPop.end(),
		Rbt::isGenome_eq(equalityThreshold));
	mergedPop.erase(end, mergedPop.end());
	m_pop.clear();
	end = (mergedPop.size() > m_size) ? (mergedPop.begin() + m_size) : mergedPop.end();
	std::copy(mergedPop.begin(), end, back_inserter(m_pop));
}
#endif // AdvaceSoft (2018)


void RbtPopulation::EvaluateRWFitness() {
  //Determine mean and variance of true scores
  RbtDouble sum(0.0);
  RbtDouble sumSq(0.0);
  for (RbtGenomeListConstIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    RbtDouble score = (*iter)->GetScore();
    sum += score;
    sumSq += score*score;
  }
  RbtDouble popSize = m_pop.size();
  m_scoreMean = sum / popSize;
  m_scoreVariance = (sumSq / popSize) - (m_scoreMean * m_scoreMean);
  RbtDouble sigma = sqrt(m_scoreVariance);
  // calculate scaled fitness values using sigma truncation
  // Goldberg page 124
  RbtDouble offset = m_scoreMean - m_c * sigma;
  RbtDouble partialSum = 0.0;
  //First set the unnormalised fitness values
  for (RbtGenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    partialSum = (*iter)->SetRWFitness(offset, partialSum);
  }
  //Now normalise so that fitness values run from 0 to 1, ready for
  //roulette wheel selection
  for (RbtGenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    (*iter)->NormaliseRWFitness(partialSum);
  }
}

