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

#ifndef PMTHit_h
#define PMTHit_h 1
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"



class PMTHit : public G4VHit
{
public:
  
  PMTHit();
  ~PMTHit();
  PMTHit(const PMTHit &right);

  const PMTHit& operator=(const PMTHit &right);
  G4int operator==(const PMTHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  
  void Draw();
  void Print();

  inline void SetDrawit(G4bool b){drawit=b;}
  inline G4bool GetDrawit(){return drawit;}

  inline void SetPhotonCount(G4int n) { photons = n; }
  inline void IncPhotonCount(){photons++;}
  inline G4int GetPhotonCount(){return photons;}

  inline void SetPMTNumber(G4int n) { pmtNumber = n; }
  inline G4int GetPMTNumber() { return pmtNumber; }

  inline void SetPMTPhysVol(G4VPhysicalVolume* physVol){this->physVol=physVol;}
  inline G4VPhysicalVolume* GetPMTPhysVol(){return physVol;}

  inline void SetPMTPos(G4ThreeVector xyz){
    pos=xyz;
  }
  
  inline G4ThreeVector GetPMTPos(){return pos;}
  
  inline void SetPMTTime(G4double t) { time = t; }
  inline G4double GetPMTTime() { return time; }

  //Initialize the arrays to store time of photons  
  inline void InitPMTs(){
    NumT=80;
    for(G4int i=0;i<NumT;i++)
     {
       TofPs[i]=-1.0;
     }
  }
  //Store a photon time
  inline void SetTofPs(G4double tofp){
      if(TofPs)
      {
       for(G4int j=0;j<NumT;j++)
          if( TofPs[j]<0.0)
             {TofPs[j]=tofp;break;}
          else  if( TofPs [j]>tofp)
                  {
                    for(G4int k=NumT-1;k>j;k--)
                     { 
                         TofPs[k]=TofPs[k-1];
                     }
                   TofPs[j]=tofp ;break;
                  }
      }      
    }
  inline void SetTscint(G4double t) { Tscint = t; }
  inline G4double GetTscint(){
         return TofPs[NumT-1];
      }
 
private:
  G4int pmtNumber;
  G4int photons;
  G4ThreeVector pos;
  G4VPhysicalVolume* physVol;
  G4bool drawit;
  G4double time;
  G4double Tscint;
  G4double TofPs[110];
  G4int NumT;
};


typedef G4THitsCollection<PMTHit> PMTHitsCollection;

extern G4Allocator<PMTHit> PMTHitAllocator;

inline void* PMTHit::operator new(size_t){
  void *aHit;
  aHit = (void *) PMTHitAllocator.MallocSingle();
  return aHit;
}

inline void PMTHit::operator delete(void *aHit){
  PMTHitAllocator.FreeSingle((PMTHit*) aHit);
}

#endif


