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
// $Id: DetectorConstruction.cc 77656 2013-11-27 08:52:57Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "ScoreParameterisation.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4SDManager.hh"

#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Construct materials
    ConstructMaterials();
    G4Material* air = G4Material::GetMaterial("G4_AIR");
    G4Material* lead = G4Material::GetMaterial("G4_Pb");
    G4Material* NaI = G4Material::GetMaterial("NaI");
    G4bool checkOverlaps = true;

    
    // geometries --------------------------------------------------------------
    // experimental hall (world volume)
    G4VSolid* worldSolid = new G4Box("worldBox",2.5*m,2.5*m,2.5*m);
    G4LogicalVolume* worldLogical= new G4LogicalVolume(worldSolid,air,"worldLogical");
    G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,false,0,checkOverlaps);
                          
    //dimensions of target disc
    
    G4double innerRadius = 0.0*cm;
    G4double outerRadius = 1.0*cm;
    G4double hz = 0.5*mm;
    G4double startAngle = 0.0*deg;
    G4double spanningAngle = 360.0*deg;      
    
    //create target disc volume
    
    G4Tubs* targetDisc = new G4Tubs("Target", innerRadius, outerRadius, hz, startAngle, spanningAngle);               
    
    G4LogicalVolume* targetLog = new G4LogicalVolume(targetDisc, lead, "Target");
    
    G4double pos_x =  0.0*meter;
	G4double pos_y =  0.0*meter;
	G4double pos_z =  0.0*meter;
	
	G4VPhysicalVolume* targetPhysical
		= new G4PVPlacement(0,
							G4ThreeVector(pos_x, pos_y, pos_z),
							targetLog,
							"Target"	,
							worldLogical,
							false,
							0);
							
							
	//dimensions of detector mother sphere
	
	G4double innerSphereRadius = 1.0*m;
	G4double outerSphereRadius = 1.01*m;
	G4double startingPhiSphere = 0.0*deg;
	G4double spanningPhiSphere = 360.0*deg;
	G4double startingThetaSphere = 0.0*deg;
	G4double spanningThetaSphere = 180.0*deg;
	
	G4VSolid* sphereSolid = new G4Sphere("SphereSolid", 
											innerSphereRadius, 
											outerSphereRadius,
											startingPhiSphere,
											spanningPhiSphere,
											startingThetaSphere,
											spanningThetaSphere);
											
	G4LogicalVolume* sphereLog = new G4LogicalVolume(sphereSolid, NaI, "SphereLogical");
	
	G4VPhysicalVolume* spherePhysical 
								= new G4PVPlacement(0,
								G4ThreeVector(pos_x,pos_y,pos_z),
								sphereLog,
								"Sphere",
								worldLogical,
								false,
								0);
										
	
	//Create scoring rings
	
	G4double deltaTheta = 6.0*deg;
	
	G4VSolid* scoreSolid = new G4Sphere("ScoreSolid",
										innerSphereRadius,
										outerSphereRadius,
										startingPhiSphere,
										spanningPhiSphere,
										startingThetaSphere,
										deltaTheta);
										
	G4LogicalVolume* scoringLog = new G4LogicalVolume(scoreSolid, NaI, "ScoringLogical");
	
	ScoreParameterisation* param = new ScoreParameterisation(innerSphereRadius,
															outerSphereRadius,
															deltaTheta);
															

	new G4PVParameterised ("ScorePhysical", scoringLog, sphereLog, kZAxis,
								30, param);

	
  
	SetupScoring(scoringLog);
  
	
    // return the world physical volume ----------------------------------------
    
    return worldPhysical;
}

void DetectorConstruction::SetupScoring(G4LogicalVolume* scoringVolume)
	{
  
    //create sensitive detector "rings"
    
    G4MultiFunctionalDetector* detector = new G4MultiFunctionalDetector("Rings");
    
    //get pointer to detector manager
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    
    //register detector with manager
    sdManager->AddNewDetector(detector);
   
    //attach detector to scoring volume
    
    scoringVolume->SetSensitiveDetector(detector);
    
    //create primitive scorer Edep
    
    G4PSEnergyDeposit* scoreEdep = new G4PSEnergyDeposit("scoreEdep");
    
    //register scorer with detector
    detector->RegisterPrimitive(scoreEdep);
    
    G4cout<<"Created G4MultiFunctionalDetector named "
        <<detector->GetName()<<", and a G4PSEnergyDeposit scorer "
        <<"named "<<scoreEdep->GetName()<<G4endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructMaterials()
{
    G4NistManager* nistManager = G4NistManager::Instance();

    // Air 
    nistManager->FindOrBuildMaterial("G4_AIR");
  
    // Argon gas
    nistManager->FindOrBuildMaterial("G4_Ar");

    // Scintillator
    // (PolyVinylToluene, C_9H_10)
    nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
    //create CSI and NaI for detector crystals
    
    G4Element* elI = new G4Element("Iodine","I", 53., 126.9*g/mole);
    G4Element* elCs = new G4Element("Cesium", "Cs", 55., 132.9*g/mole);
    G4Element* elNa = new G4Element("Sodium", "Na", 11.0, 22.9*g/mole);
    
    G4Material* CsI = new G4Material("CsI", 4.51*g/cm3, 2);
    
    CsI->AddElement(elCs,1);
    CsI->AddElement(elI,1);
    
    G4Material* NaI = new G4Material("NaI", 3.67*g/cm3, 2);
    
    NaI->AddElement(elNa,1);
    NaI->AddElement(elI,1);
    
    
     nistManager->FindOrBuildMaterial("G4_Pb");

    
    

    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
