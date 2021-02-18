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

#include "RbtBiMolWorkSpace.h"
#include "RbtBaseSF.h"
#include "RbtAnnotationHandler.h"
#include "RbtBaseMolecularFileSource.h"

////////////////////////////////////////
//Constructors/destructors


#if 1 // AdvanceSoft (2018)
#include <mpi.h>
bool b_receptor(false);
#endif // AdvanceSoft(2018)

#if 1 // Xability (2017)
RbtBiMolWorkSpace::RbtBiMolWorkSpace() : RbtWorkSpace(2), 
                                         m_spLigandAtc(RbtModelPtr()) {
  _RBTOBJECTCOUNTER_CONSTR_("RbtBiMolWorkSpace");
}
#endif // Xability (2017)

RbtBiMolWorkSpace::~RbtBiMolWorkSpace() {
  _RBTOBJECTCOUNTER_DESTR_("RbtBiMolWorkSpace");
}
	
////////////////////////////////////////
//Public methods
////////////////

//Model handling
RbtModelPtr RbtBiMolWorkSpace::GetReceptor() const {return GetModel(RECEPTOR);}
RbtModelPtr RbtBiMolWorkSpace::GetLigand() const {return GetModel(LIGAND);}
RbtModelList RbtBiMolWorkSpace::GetSolvent() const {return GetModels(SOLVENT);}
RbtBool RbtBiMolWorkSpace::hasSolvent() const {
	return GetNumModels() > SOLVENT;
}

void RbtBiMolWorkSpace::SetReceptor(RbtModelPtr spReceptor) {
  SetModel(RECEPTOR,spReceptor);
  //DM 30 Oct 2001 - no need to clear the ligand population, there may be scenarios where we want
  //to continue a ligand docking with a revised receptor model
}

void RbtBiMolWorkSpace::SetLigand(RbtModelPtr spLigand) {
 //DM 30 Oct 2001 - invalidate the current population
  ClearPopulation();
  //DM 26 Aug 2002 - better to call SetModel AFTER defining the ligand as flexible
  SetModel(LIGAND,spLigand);
}


void RbtBiMolWorkSpace::SetSolvent(RbtModelList solventList) {
	if (hasSolvent()) {
		SetModels(SOLVENT,solventList);
	}
	else {
		AddModels(solventList);
	}
}

void RbtBiMolWorkSpace::RemoveSolvent() {
	if (hasSolvent()) {
		RemoveModels(SOLVENT);
	}
}

void RbtBiMolWorkSpace::UpdateModelCoordsFromChromRecords(RbtBaseMolecularFileSource* pSource,
                                            RbtInt iTrace) {
    RbtInt nModels = GetNumModels();
    for (RbtInt iModel = 0; iModel < nModels; iModel++) {
        //if (iModel == LIGAND) continue;
        RbtModelPtr spModel = GetModel(iModel);
        if (spModel.Ptr()) {
            RbtChromElementPtr spChrom = spModel->GetChrom();
            if (spChrom.Ptr()) {
                RbtInt chromLength = spChrom->GetLength();
                if (chromLength > 0) {
                    ostringstream ostr;
                    ostr << "CHROM." << iModel;
                    RbtString chromField = ostr.str();
                    if (pSource->isDataFieldPresent(chromField)) {
                        //TODO: Move this code to RbtVariant class
                        //Concatenate the multi-record value into a single string
                        RbtString chromRecord = Rbt::ConvertListToDelimitedString(
                                                        pSource->GetDataValue(chromField), ",");
                        //Now split into string values and convert to doubles
                         RbtStringList chromValues = Rbt::ConvertDelimitedStringToList(
                                                        chromRecord, ",");
                         if (chromLength == chromValues.size()) {
                            RbtDoubleList chromVec;
                            for (RbtStringListConstIter iter = chromValues.begin();
                                                        iter != chromValues.end();
                                                        ++iter) {
                                RbtString chromValue(*iter);
                                chromVec.push_back(atof(chromValue.c_str()));
                            }
                            spChrom->SetVector(chromVec);
                            spChrom->SyncToModel();
                            spModel->UpdatePseudoAtoms();
                         }
                         else if (iTrace > 0) {
                            cout << "Mismatched chromosome sizes for model #" << iModel << endl;
                            cout << chromField << " record in " << pSource->GetFileName()
                                 << " has " << chromValues.size() << " elements" << endl;
                            cout << "Expected number of elements is " << chromLength << endl;
                            cout << "Model chromosome not updated" << endl;
                         }
                    }
                }
            }
        }
    }
}

//Model I/O
//Saves ligand to file sink
void RbtBiMolWorkSpace::Save(RbtBool bSaveScores) {
  SaveLigand(GetSink(),bSaveScores);

#if 1 // AdvanceSoft (2018)
  if (b_receptor) {
	  SaveReceptor(GetSink_rec(), bSaveScores);
  }
#endif // AdvanceSoft(2018)

}

//Model I/O
//Saves ligand to history file sink
void RbtBiMolWorkSpace::SaveHistory(RbtBool bSaveScores) {
  SaveLigand(GetHistorySink(),bSaveScores);
}

//Private method to save ligand to file sink, optionally with score attached
void RbtBiMolWorkSpace::SaveLigand(RbtMolecularFileSinkPtr spSink, RbtBool bSaveScores) {
  if (spSink.Null()) //Check we have a valid sink
    return;
  RbtModelPtr spLigand(GetLigand());
  if (spLigand.Null()) //Check we have a ligand
    return;
  
  RbtBaseSF* pSF(GetSF());
  //If we have a scoring function, clear any current score data
  if (pSF) {
    RbtString strSFName(pSF->GetFullName());
    spLigand->ClearAllDataFields(strSFName);
    spLigand->ClearAllDataFields(RbtAnnotationHandler::_ANNOTATION_FIELD);
  }
  //Save current score components as model data fields
  if (bSaveScores && pSF) {
    RbtStringVariantMap scoreMap;
    pSF->ScoreMap(scoreMap);
    for (RbtStringVariantMapConstIter vIter = scoreMap.begin(); vIter != scoreMap.end(); vIter++) {
      spLigand->SetDataValue((*vIter).first,(*vIter).second);
    }
    //Save the chromosome values for all models for later retrieval
    spLigand->ClearAllDataFields("CHROM.");
    RbtInt nModels = GetNumModels();
    for (RbtInt iModel = 0; iModel < nModels; iModel++) {
        RbtModelPtr spModel = GetModel(iModel);
        if (spModel.Ptr()) {
            RbtChromElementPtr spChrom = spModel->GetChrom();
            if (spChrom.Ptr() && (spChrom->GetLength() > 0) ) {
                spChrom->SyncFromModel();
                RbtDoubleList chromVec;
                spChrom->GetVector(chromVec);
                ostringstream ostr;
                ostr << "CHROM." << iModel;
                spLigand->SetDataValue(ostr.str(),RbtVariant(chromVec,72,8));
            }
        }
    }
  }
  spSink->SetModel(spLigand);
  
#if 1 // Xability (2017)
  spSink->SetModelAtc(m_spLigandAtc);
#endif // Xability (2017)
  
  if (hasSolvent()) {
  	spSink->SetSolvent(GetSolvent());
  }
  else {
  	spSink->SetSolvent(RbtModelList());
  }
#if 1 // AdvanceSoft (2018)
  if (!b_parallel_mmgbpbsa) {
	  spSink->Render();
  }
  else {
	  int iMPIRank = 0;
	  MPI_Comm_rank(MPI_COMM_WORLD, &iMPIRank);
	  if (iMPIRank == 0) {
		  spSink->Render();
	  }
  }
#endif // AdvanceSoft (2018)
}

#if 1 // AdvanceSoft (2018)
void RbtBiMolWorkSpace::SaveReceptor(RbtMolecularFileSinkPtr spSink, RbtBool bSaveScores) {
	if (spSink.Null()) //Check we have a valid sink
		return;

	RbtModelPtr spLigand(GetLigand());

	RbtModelPtr spReceptor(GetReceptor());
	if (spReceptor.Null()) //Check we have a ligand
		return;

#if 0 // original code
	RbtModelPtr spLigand(GetLigand());
	if (spLigand.Null()) //Check we have a ligand
		return;
#endif // original code

	RbtBaseSF* pSF(GetSF());
	//If we have a scoring function, clear any current score data
	if (pSF) {
		RbtString strSFName(pSF->GetFullName());
		spReceptor->ClearAllDataFields(strSFName);
		spReceptor->ClearAllDataFields(RbtAnnotationHandler::_ANNOTATION_FIELD);
	}
	//Save current score components as model data fields
	if (bSaveScores && pSF) {
		RbtStringVariantMap scoreMap;
		pSF->ScoreMap(scoreMap);
		for (RbtStringVariantMapConstIter vIter = scoreMap.begin(); vIter != scoreMap.end(); vIter++) {
			spReceptor->SetDataValue((*vIter).first, (*vIter).second);
		}
		//Save the chromosome values for all models for later retrieval
		spReceptor->ClearAllDataFields("CHROM.");
		RbtInt nModels = GetNumModels();
		for (RbtInt iModel = 0; iModel < nModels; iModel++) {
			RbtModelPtr spModel = GetModel(iModel);
			if (spModel.Ptr()) {
				RbtChromElementPtr spChrom = spModel->GetChrom();
				if (spChrom.Ptr() && (spChrom->GetLength() > 0)) {
					spChrom->SyncFromModel();
					RbtDoubleList chromVec;
					spChrom->GetVector(chromVec);
					ostringstream ostr;
					ostr << "CHROM." << iModel;
					spReceptor->SetDataValue(ostr.str(), RbtVariant(chromVec, 72, 8));
				}
			}
		}
	}
	spSink->SetModel(spReceptor);
	spSink->SetModel_lig(spLigand);

	if (hasSolvent()) {
		spSink->SetSolvent(GetSolvent());
	}
	else {
		spSink->SetSolvent(RbtModelList());
	}
	if (!b_parallel_mmgbpbsa) {
		spSink->Render_rec();
	}
	else {
		int iMPIRank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &iMPIRank);
		if (iMPIRank == 0) {
			spSink->Render_rec();
		}
	}
}
#endif // AdvanceSoft(2018)