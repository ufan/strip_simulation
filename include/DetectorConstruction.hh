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
#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Box;
class G4Tubs;
class PMTSD;
class ScintSD;
class G4Sphere;
class DetectorMessenger;

#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();

  //Functions to modify the geometry

  void SetDefaults();
  void SetGuideLongth(G4double avalue){guide_long=avalue;}

  //Get values
  
  G4double GetHousingThickness(){return d_mtl;}
  
  
    
  //rebuild the geometry based on changes. must be called
  void UpdateGeometry();
  G4bool GetUpdated(){return updated;}

 
private:

  void DefineMaterials();
  G4VPhysicalVolume* ConstructDetector();
  G4bool updated;
   
  G4Box* experimentalHall_box;
  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;

  //Materials & Elements
 
  G4Material* Al;
  G4Element* N;
  G4Element* O;
  G4Material* Air;
  G4Material* Vacuum;
  G4Element* C;
  G4Element* H;
  G4Material* Glass;
  G4Material* Scintillator;


  //Geometry
  
  G4double d_mtl;
  G4double guide_long;
  G4double outerRadius_pmt;

  G4MaterialPropertiesTable *myMPT1;
  G4MaterialPropertiesTable *myMPT2 ;
  G4MaterialPropertiesTable *myMPT3;

  //Sensitive Detectors
  static ScintSD* scint_SD;
  static PMTSD* pmt_SD;

  DetectorMessenger* fDetectorMessenger;

};

#endif




