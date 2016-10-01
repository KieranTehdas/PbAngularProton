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
// $Id: PrimaryGeneratorAction.cc 77781 2013-11-28 07:54:07Z gcosmo $
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),     
  fParticleGun(0)
{
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);
    
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    //fPositron = particleTable->FindParticle(particleName="e+");
    //fMuon = particleTable->FindParticle(particleName="mu+");
    //fPion = particleTable->FindParticle(particleName="pi+");
    //fKaon = particleTable->FindParticle(particleName="kaon+");
    //fProton = particleTable->FindParticle(particleName="proton");
    G4ParticleDefinition* photon = particleTable->FindParticle(particleName="gamma");
    fProton = particleTable->FindParticle(particleName="proton");
    
    G4double kineticEnergy = 50.0*MeV;
    
    // default particle kinematics
    fParticleGun->SetParticleDefinition(fProton);
	fParticleGun->SetParticleEnergy(kineticEnergy);
    
    // define commands for this class
    //DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 
 //set position of particle source   
  
	G4double source_x = 0.0*m;
	G4double source_y = 0.0*m;
	G4double source_z = 0.0*m;
  
	fParticleGun->SetParticlePosition(G4ThreeVector(source_x,source_y,source_z));

//create an isotropic source

	G4double cosTheta = -1.0 + 2.0*G4UniformRand();
	G4double phi = CLHEP::twopi*G4UniformRand();
	G4double sinTheta = sqrt(1.0 - cosTheta*cosTheta);

	fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(sinTheta*cos(phi),
								sinTheta*sin(phi),
							    cosTheta));

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
