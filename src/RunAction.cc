//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "RunAction.hh"
#include "Analysis.hh"
#include "RecorderBase.hh"
//#include "G4MPImanager.hh"
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
RunAction::RunAction(RecorderBase* r)
  :recorder(r)
{}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
RunAction::~RunAction()
{}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void RunAction::BeginOfRunAction(const G4Run* aRun){
 
   Analysis* myana = Analysis::GetAnalysis();
  myana-> Clear();
 G4cout << "### Run " <<  aRun->GetRunID() << " start." << G4endl;
  if(recorder)recorder->RecordBeginOfRun(aRun);
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void RunAction::EndOfRunAction(const G4Run* aRun){

  //G4int rank = G4MPImanager::GetManager()-> GetRank();
    G4int rank=0;
  char str[64];
  G4int id;
  id=aRun->GetRunID();
   //if(id>6)id=id+1;
  sprintf(str, "%03d-%03d.root", id,rank);
  G4String fname(str);
  Analysis* myana = Analysis::GetAnalysis();
  myana-> Save(fname);

  if(recorder)recorder->RecordEndOfRun(aRun);
  G4cout<< ">>> end of run " <<id << G4endl;
}
