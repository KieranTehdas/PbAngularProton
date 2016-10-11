//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id:$
// GEANT4 tag $Name:$
//
// T. Aso        Original author
//
#include "ScoreParameterisation.hh"
#include "G4Sphere.hh"
#include "G4VPhysicalVolume.hh"

ScoreParameterisation::ScoreParameterisation(const G4double minRadius, 
							     const G4double maxRadius,
							     const G4double deltaTheta)
  :G4VPVParameterisation()
  ,fMinRadius(minRadius)
  ,fMaxRadius(maxRadius)
  ,fDeltaTheta(deltaTheta)
{}

ScoreParameterisation::~ScoreParameterisation() {}

void ScoreParameterisation::ComputeDimensions(G4Sphere& sphere, 
						      const G4int copyNo, 
						      const G4VPhysicalVolume*) const
{
  G4double startTheta = copyNo*fDeltaTheta;
  sphere.SetInsideRadius(fMinRadius);
  sphere.SetOuterRadius(fMaxRadius);
  sphere.SetStartPhiAngle(0.*CLHEP::deg);
  sphere.SetDeltaPhiAngle(360.*CLHEP::deg);
  sphere.SetStartThetaAngle(startTheta);
  sphere.SetDeltaThetaAngle(fDeltaTheta);
}





