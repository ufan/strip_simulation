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
#include "PMTSD.hh"
#include "PMTHit.hh"
#include "DetectorConstruction.hh"
//#include "UserTrackInformation.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include <iostream>
#include <fstream>
using namespace std;
extern ofstream fout;
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
PMTSD::PMTSD(G4String name)
  :G4VSensitiveDetector(name),pmtHitCollection(0),pmtPositionsX(0)
  ,pmtPositionsY(0),pmtPositionsZ(0)
{
  collectionName.insert("pmtHitCollection");
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
PMTSD::~PMTSD()
{}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void PMTSD::Initialize(G4HCofThisEvent* HCE){
  pmtHitCollection = new PMTHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  //Store collection with event and keep ID
  static G4int HCID = -1;
  if(HCID<0){ 
    HCID = GetCollectionID(0); 
  }
  HCE->AddHitsCollection( HCID, pmtHitCollection );
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4bool PMTSD::ProcessHits(G4Step* ,G4TouchableHistory* ){
  return false;
}

//Generates a hit and uses the postStepPoint's mother volume replica number
//PostStepPoint because the hit is generated manually when the photon is
//absorbed by the photocathode
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4bool PMTSD::ProcessHits_constStep(const G4Step* aStep,
				       G4TouchableHistory* ){

  //G4int pmtThreshold=1;
  G4ThreeVector xyz;
  G4double time=-1.0*ns;
  time=aStep->GetPreStepPoint()->GetGlobalTime();

  //G4double energy=0;
  //energy=aStep->GetTrack()->GetKineticEnergy();
  //need to know if this is an optical photon
  if(aStep->GetTrack()->GetDefinition() 
     != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
 
  //User replica number 1 since photocathode is a daughter volume
  //to the pmt  pmt is a daughter to the guide  which was replicated
  G4int pmtNumber=
    aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(2);
  G4VPhysicalVolume* physVol=
    aStep->GetPostStepPoint()->GetTouchable()->GetVolume(2);

  //Find the correct hit collection
  G4int n=pmtHitCollection->entries();
  PMTHit* hit=NULL;
  for(G4int i=0;i<n;i++){
    if((*pmtHitCollection)[i]->GetPMTNumber()==pmtNumber){
      hit=(*pmtHitCollection)[i];
      break;
    }
  }
  
  if(hit==NULL){//this pmt wasnt previously hit in this event
    hit = new PMTHit(); //so create new hit
    hit->SetPMTNumber(pmtNumber);
    hit->SetPMTPhysVol(physVol);
    hit->InitPMTs();
    hit->SetPMTTime(time) ; 
    pmtHitCollection->insert(hit);
    xyz=G4ThreeVector((*pmtPositionsX)[pmtNumber],(*pmtPositionsY)[pmtNumber],
		   (*pmtPositionsZ)[pmtNumber]);
    hit->SetPMTPos(xyz);
  }

  hit->IncPhotonCount(); //increment hit for the selected pmt
 
  if(time<hit->GetPMTTime())
  hit->SetPMTTime(time) ; 
  hit->SetTofPs(time);
  hit->SetDrawit(true);
   
 // fout<<pmtNumber<<"\t"<<energy/eV<<"\t\t"<<time/ns<<G4endl;
    
  return true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void PMTSD::EndOfEvent(G4HCofThisEvent* ){
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void PMTSD::clear(){
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void PMTSD::DrawAll(){
} 

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void PMTSD::PrintAll(){
} 

