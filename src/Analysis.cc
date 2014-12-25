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
/// @file Analysis.cc
/// @brief Define histograms

#include "Analysis.hh"
#include "AnaMessenger.hh"
#include "G4SystemOfUnits.hh"

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

Analysis* Analysis::myanalysis = new Analysis();

// --------------------------------------------------------------------------
Analysis::Analysis()
{
  fMessenger = new AnaMessenger(this);
  fAnaName   = "bar";
  p1=-1;p2=-1;p3=-1;p4=-1;
  // ROOT style
//gROOT-> Reset();

  // define tree
  mytree=new TTree("T","An photon number tree with 4 branches");
  mytree->Branch("p1",&p1,"p1/I");
  mytree->Branch("p2",&p2,"p2/I");
  mytree->Branch("p3",&p3,"p3/I");
  mytree->Branch("p4",&p4,"p4/I");
  mytree->Branch("ed",&edpos,"ed/D");
}

// --------------------------------------------------------------------------
Analysis::~Analysis()
{
  delete mytree;
  delete fMessenger;
  myanalysis = NULL;
}

// --------------------------------------------------------------------------
Analysis* Analysis::GetAnalysis()
{
  if ( myanalysis == NULL ) {
    myanalysis = new Analysis();
  }

  return myanalysis;
}

// --------------------------------------------------------------------------
void Analysis::Update()
{
  return;
}

// --------------------------------------------------------------------------
void Analysis::Clear()
{
  mytree-> Reset();
 
  return;
}

// --------------------------------------------------------------------------
void Analysis::Save(const G4String& fname)
{
  G4String name;
  name="./edpos/"+fAnaName+fname;
  //name=fAnaName+fname;
  TFile* file = new TFile(name.c_str(),
                          "RECREATE", "Geant4 ROOT analysis");
  mytree-> Write();
  
  file-> Close();

  delete file;

  return;
}

// --------------------------------------------------------------------------
void Analysis::FillIncident(G4int ph1,G4int ph2,G4int ph3,G4int ph4,G4double ed)
{
  p1=ph1;p2=ph2;p3=ph3;p4=ph4;edpos=ed;
  mytree->Fill();
  p1=-1;p2=-1;p3=-1;p4=-1;edpos=-1.0;

}

void Analysis::SetFileName(const G4String& nam) 
{
  fAnaName = nam;
}
