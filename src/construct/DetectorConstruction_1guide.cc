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
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "PMTSD.hh"
#include "ScintSD.hh"


#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4UImanager.hh"
#include "G4NistManager.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"


ScintSD* DetectorConstruction::scint_SD;
PMTSD* DetectorConstruction::pmt_SD;

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
DetectorConstruction::DetectorConstruction()

{
 fDetectorMessenger = new DetectorMessenger(this);
 //fDetectorMessenger =0;
  SetDefaults();
  
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
DetectorConstruction::~DetectorConstruction(){
 delete fDetectorMessenger;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void DetectorConstruction::DefineMaterials(){
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;

  //***Elements
  H = new G4Element("H", "H", z=1., a=1.01*g/mole);
  C = new G4Element("C", "C", z=6., a=12.01*g/mole);
  N = new G4Element("N", "N", z=7., a= 14.01*g/mole);
  O = new G4Element("O"  , "O", z=8., a= 16.00*g/mole);
  
  //***Materials
  //Aluminum
  Al = new G4Material("Al",z=13.,a=26.98*g/mole,density=2.7*g/cm3);
  //Vacuum
  Vacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
			  density=universe_mean_density,kStateGas,0.1*kelvin,
			  1.e-19*pascal); 
  //Air
  Air = new G4Material("Air", density= 1.29*mg/cm3, 2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);
  //Glass
  Glass = new G4Material("Glass", density=1.032*g/cm3,2);
  Glass->AddElement(C,91.533*perCent);
  Glass->AddElement(H,8.467*perCent);

  //***Material properties tables

  //add by llx

  
  Scintillator = new G4Material("Scintillator",1.032*g/cm3,2);
  Scintillator->AddElement(C,9);
   Scintillator->AddElement(H,10);
  Scintillator->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
  
 const G4int NUMENTRIES = 6;

  G4double PP[NUMENTRIES] =
    { 2.0E-9*GeV, 2.783E-9*GeV,  2.924E-9*GeV, 2.995E-9*GeV, 3.064E-9*GeV ,3.14*eV};

  G4double RINDEX1[NUMENTRIES]= 
    { 1.58, 1.58, 1.58, 1.58, 1.58 ,1.58};

   G4double RINDEX2[NUMENTRIES] =
    { 1.00, 1.00, 1.00, 1.00, 1.00,1.0};

  G4double RINDEX3[NUMENTRIES] =
    { 1.49, 1.49, 1.49, 1.49, 1.49,1.49};

  G4double ABSORPTION1[NUMENTRIES] = 
    { 3.8*m,3.8*m,3.8*m,3.8*m,3.8*m,3.8*m};

  G4double Glass_AbsLength[NUMENTRIES]=
    {420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,4.2*m};

  G4double INTENSITY[NUMENTRIES] =
    {1., 1., 1., 1., 1.,1.0};

  myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX", PP, RINDEX1, NUMENTRIES);
  myMPT1->AddProperty("SLOWCOMPONENT", PP, INTENSITY, NUMENTRIES);
  myMPT1->AddProperty("ABSLENGTH",PP, ABSORPTION1, NUMENTRIES);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD",10000./MeV );
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 0.9*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",14.2*ns);
  myMPT1->AddConstProperty("YIELDRATIO",1.0);
  Scintillator->SetMaterialPropertiesTable(myMPT1);

  myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", PP, RINDEX2, NUMENTRIES);
  Air->SetMaterialPropertiesTable(myMPT2);

  myMPT3 = new G4MaterialPropertiesTable();
  myMPT3->AddProperty("ABSLENGTH",PP,Glass_AbsLength,NUMENTRIES); 
  myMPT3->AddProperty("RINDEX", PP, RINDEX3, NUMENTRIES);
  Glass->SetMaterialPropertiesTable(myMPT3);

  
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4VPhysicalVolume* DetectorConstruction::Construct(){
  DefineMaterials();
  return ConstructDetector();
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
  
  //Create experimental hall
  experimentalHall_box
    = new G4Box("expHall_box",0.4*m,0.4*m,1.5*m);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,
                                             Air,"expHall_log",0,0,0);
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
			      experimentalHall_log,"expHall",0,false,0);

  experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);

//paddle 
  G4Box* paddle_box = new G4Box("paddle_box",15*cm,15*cm,
			    140*cm);
  G4LogicalVolume* paddle_log = new G4LogicalVolume(paddle_box,
				       G4Material::GetMaterial("Air"),
				      "paddle_log",0,0,0); 
  G4PVPlacement* paddle_phys = new G4PVPlacement(0,G4ThreeVector(),paddle_log,"paddle",
				   experimentalHall_log,false,0);  

  paddle_log->SetVisAttributes(G4VisAttributes::Invisible);
 //guide
 ///////////////////////////////////////////////////////////////////////////
 //----------------------- head_Con

  G4double  pRmin1 = 0*mm;
  G4double  pRmax1 =  sqrt(18)*5*mm;//sqrt(109)*5*mm;
  G4double  pRmin2 = 0 *mm;
  G4double  pRmax2 = 26*mm;
  G4double  pDConz =guide_long/2.0;
  G4Cons* head_Con=new G4Cons("head_Con",pRmin1,pRmax1,pRmin2,pRmax2,pDConz,0.0,360.0*deg);

 //----------------------- head_Trd

  G4double  dx1 = 15*mm;
  G4double  dx2 = 26*mm;
  G4double  dy1 = 15*mm;
  G4double  dy2 = 26*mm;
  G4double  pDTrdz = pDConz;
  G4Trd* head_Trd=new G4Trd("head_Trd",dx1,dx2,dy1,dy2,pDTrdz);
//-------------------------------Con*Trd

  G4RotationMatrix* sRot = new G4RotationMatrix;
  sRot->rotateZ(0.*rad);
  G4ThreeVector sTrans(0.,0.,0.);
  G4IntersectionSolid* head = new G4IntersectionSolid("head", head_Con, head_Trd, sRot, sTrans);

//------------------------------ end

  G4double  pDcz = 2.5 *mm;
  G4double  pSPhi = 0. *deg;
  G4double  pDPhi = 360. *deg;
  G4Tubs* end_Cyl=new G4Tubs("end_Cyl",pRmin2,pRmax2,pDcz,pSPhi,pDPhi);

 //---------------------------- head+end

  G4double mPos_x = 0. *mm;
  G4double mPos_y = 0. *mm;
  G4double mPos_z =pDConz+pDcz;

  G4RotationMatrix* mRot = new G4RotationMatrix;
  mRot->rotateZ(0.*rad);
  G4ThreeVector mTrans(mPos_x,mPos_y,mPos_z);
  G4UnionSolid* guide= new G4UnionSolid("head_end", head, end_Cyl, mRot, mTrans);

//----------------------------scintillator_middle
  G4double  scintz = 15.0*mm;
  G4Box* scint_middle = new G4Box("scint_middle",dx1,dy1,scintz);

//--------------------------scintillator_middle+all=scintillator_union

  G4double suPos_x = 0. *mm;
  G4double suPos_y = 0. *mm;
  G4double suPos_z = pDConz+scintz-0.001*mm;

  G4RotationMatrix* suRot = new G4RotationMatrix;
  suRot->rotateX(0.*deg);
  G4ThreeVector suTrans(suPos_x,suPos_y,suPos_z);
  G4UnionSolid* scint_union= new G4UnionSolid("scint_union", scint_middle, guide, suRot, suTrans);

  G4LogicalVolume* scint_log = new G4LogicalVolume(scint_union,Scintillator,
				    "scint_log",0,0,0);
  G4PVPlacement* scint_phys = new G4PVPlacement(0,G4ThreeVector(),scint_log,"scintillator",
				   paddle_log,false,0); 
  //-----------------------  guide

  G4double guidePos_x = 0. *mm;
  G4double guidePos_y = 0. *mm;
  G4double guidePos_z = suPos_z;
  G4ThreeVector  gpos=G4ThreeVector(guidePos_x,guidePos_y,guidePos_z);

  G4RotationMatrix gRot;
  gRot.rotateX(0.*deg);
  G4LogicalVolume* guide_log = new G4LogicalVolume(guide,G4Material::GetMaterial("Glass"),"guide_log",0,0,0);
  G4PVPlacement* guideup_phys= new G4PVPlacement(G4Transform3D(gRot,gpos),
                                                 guide_log,"guide",scint_log,false,0);
  
 /////////////////////////////////////////////////////////////
   //****************** Build PMTs
    G4double innerRadius_pmt = 0.*cm;
    G4double outerRadius_pmt = 23*mm;
    G4double height_pmt = d_mtl/2.;
    G4double startAngle_pmt = 0.*deg;
    G4double spanningAngle_pmt = 360.*deg;
    
    G4Tubs* pmt = new G4Tubs("pmt_tube",innerRadius_pmt,outerRadius_pmt,
		     height_pmt,startAngle_pmt,spanningAngle_pmt);
    
    //the "photocathode" is a metal slab at the back of the glass that
    //is only a very rough approximation of the real thing since it only
    //absorbs or detects the photons based on the efficiency set below
    G4Tubs* photocath = new G4Tubs("photocath_tube",innerRadius_pmt,outerRadius_pmt,
			   height_pmt/2.0,startAngle_pmt,spanningAngle_pmt);
    
    G4LogicalVolume* pmt_log = new G4LogicalVolume(pmt,G4Material::GetMaterial("Glass"),
				  "pmt_log");
    G4LogicalVolume* photocath_log = new G4LogicalVolume(photocath,
					G4Material::GetMaterial("Al"),
					"photocath_log");

    
     new G4PVPlacement(0,G4ThreeVector(0,0,height_pmt/2),
				       photocath_log,"photocath",
				       pmt_log,false,0);
    
    G4ThreeVector  phpos=G4ThreeVector(0.*mm,0.0*mm,pDConz+pDcz);

    G4RotationMatrix phRot;
    phRot.rotateX(0.*deg);

    G4PVPlacement* pmt_phys = new G4PVPlacement(G4Transform3D(phRot,phpos),pmt_log,"PMT",guide_log,false,0);

   

    G4VisAttributes* pmt_va = new G4VisAttributes();
    pmt_va->SetForceSolid(true);
    pmt_log->SetVisAttributes(pmt_va);

    //---pmt sensitive detector
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    if(!pmt_SD){
      pmt_SD = new PMTSD("/Det/pmtSD");
      SDman->AddNewDetector(pmt_SD);
      //Created here so it exists as pmts are being placed
    }
   G4ThreeVector xyz;
    pmt_SD->InitPMTs(4); //let pmtSD know # of pmts
    pmt_SD->SetPMTPos(0,0,0,0);
    //pmt_SD->SetPMTPos(1,0*mm,0,0);
    //**********Setup Sensitive Detectors***************
    if(!scint_SD){//determine if it has already been created
      scint_SD = new ScintSD("/Det/scintSD");
      SDman->AddNewDetector(scint_SD);    
    }
   // scint_log->SetSensitiveDetector(scint_SD);
    
    //sensitive detector is not actually on the photocathode.
    //processHits gets done manually by the stepping action.
    //It is used to detect when photons hit and get absorbed&detected at the
    //boundary to the photocathode (which doesnt get done by attaching it to a
    //logical volume.
    //It does however need to be attached to something or else it doesnt get
    //reset at the begining of events
    photocath_log->SetSensitiveDetector(pmt_SD);

  const G4int num = 2;
  G4double Ephoton[num] = {2.0*eV, 7.14*eV};
  G4double Air_RIND[num] = { 1.,1.};
  G4double Glass_RIND[num]={1.49,1.49};

  //**Photocathode surface properties
  G4double photocath_EFF[num]={0.8,0.8}; //Enables 'detection' of photons
  G4double photocath_ReR[num]={1.92,1.92};
  G4double photocath_ImR[num]={1.69,1.69};
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY",Ephoton,photocath_EFF,num);
  photocath_mt->AddProperty("REALRINDEX",Ephoton,photocath_ReR,num);
  photocath_mt->AddProperty("IMAGINARYRINDEX",Ephoton,photocath_ImR,num);
  G4OpticalSurface* photocath_opsurf= 
    new G4OpticalSurface("photocath_opsurf",glisur,polished,
			 dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

  new G4LogicalSkinSurface("photocath_surf",photocath_log,photocath_opsurf);
//-------------------------------- the surface between Air and scint
//--------------------------------slice BC408
  
  G4double PSReflectivity[num] = {0.69,0.69};
  G4double PSEfficiency[num] = {0.,0.}; 
  G4double PSSpecularSpike[num]={0.1,0.1};
  G4double PSSpecularLobe[num]={0.9,0.9};
  G4double PSBackscatter[num]={0.,0.};
  G4MaterialPropertiesTable* mySTS = new G4MaterialPropertiesTable();
  mySTS->AddProperty("REFLECTIVITY", Ephoton, PSReflectivity, num);
  mySTS->AddProperty("EFFICIENCY", Ephoton, PSEfficiency, num);
  mySTS->AddProperty("RINDEX",Ephoton,Air_RIND,num);
  mySTS->AddProperty("SPECULARLOBECONSTANT",Ephoton,PSSpecularLobe,num);
  mySTS->AddProperty("SPECULARSPIKECONSTANT",Ephoton,PSSpecularSpike,num);
  mySTS->AddProperty("BACKSCATTERCONSTANT",Ephoton,PSBackscatter,num);
  G4OpticalSurface* PSSurface = new G4OpticalSurface("PSSurface");
  G4double PSSigma_alpha=0.07;
  PSSurface->SetMaterialPropertiesTable(mySTS);
  PSSurface->SetType(dielectric_LUT);
  PSSurface->SetFinish(polishedteflonair);
  PSSurface->SetModel(LUT);
  PSSurface->SetSigmaAlpha(PSSigma_alpha);
  new G4LogicalBorderSurface("PSSurface",paddle_phys,experimentalHall_phys,PSSurface); 
  new G4LogicalBorderSurface("PSSurface",scint_phys,paddle_phys,PSSurface); 
 
//--------------------------------scint guide
  G4double GuideReflectivity[num] = {0.69,0.69};
  G4double GuideEfficiency[num] = {0.0,0.0}; 
  G4double GuideSpecularSpike[num]={0.1,0.1};
  G4double GuideSpecularLobe[num]={0.9,0.9};
  G4double GuideBackscatter[num]={0.,0.};
  G4MaterialPropertiesTable* mySTG = new G4MaterialPropertiesTable();
  mySTG->AddProperty("REFLECTIVITY", Ephoton, GuideReflectivity, num);
  mySTG->AddProperty("EFFICIENCY", Ephoton, GuideEfficiency, num);
  mySTG->AddProperty("RINDEX",Ephoton,Air_RIND,num);
  mySTG->AddProperty("SPECULARLOBECONSTANT",Ephoton,GuideSpecularLobe,num);
  mySTG->AddProperty("SPECULARSPIKECONSTANT",Ephoton,GuideSpecularSpike,num);
  mySTG->AddProperty("BACKSCATTERCONSTANT",Ephoton,GuideBackscatter,num);
  G4OpticalSurface* GuideSurface = new G4OpticalSurface("GuideSurface");
  G4double GuideSigma_alpha=0.05;
  GuideSurface->SetMaterialPropertiesTable(mySTG);
  GuideSurface->SetType(dielectric_dielectric);
  GuideSurface->SetFinish(groundbackpainted);
  GuideSurface->SetModel(unified);
  GuideSurface->SetSigmaAlpha(GuideSigma_alpha);
  new G4LogicalBorderSurface("GuideSurfaceU",guideup_phys,scint_phys,GuideSurface); 
  
   //---------------------- between guide and pmt
  G4double glassReflectivity[num] = {0.1,0.1};
  G4double glassEfficiency[num] = {0.9,0.9};
  G4double glassSpecularSpike[num]={1,1};
  G4MaterialPropertiesTable* mySTGD = new G4MaterialPropertiesTable();
  mySTGD->AddProperty("REFLECTIVITY", Ephoton, glassReflectivity, num);
  mySTGD->AddProperty("EFFICIENCY", Ephoton, glassEfficiency, num);
  mySTGD->AddProperty("RINDEX", Ephoton, Glass_RIND, num); 
  mySTGD->AddProperty("SPECULARSPIKECONSTANT",Ephoton,glassSpecularSpike,num); 
  G4OpticalSurface* OppmtSurfaceGD = new G4OpticalSurface("GDSurface");
  OppmtSurfaceGD->SetMaterialPropertiesTable(mySTGD);
  OppmtSurfaceGD->SetType(dielectric_dielectric);
  OppmtSurfaceGD->SetFinish(polished);
  OppmtSurfaceGD->SetModel(unified);
  new G4LogicalBorderSurface("GDSurfaceU",guideup_phys,pmt_phys,OppmtSurfaceGD);
 
  return experimentalHall_phys;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void DetectorConstruction::SetDefaults(){
  //Resets to default values
  d_mtl=0.0635*cm;
  guide_long=30*mm;
  outerRadius_pmt = 2.3*cm;

  updated=true;
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void DetectorConstruction::UpdateGeometry(){

  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();
  G4SurfaceProperty::CleanSurfacePropertyTable();

  //define new one
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  updated=false;
}







