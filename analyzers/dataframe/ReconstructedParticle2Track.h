
#ifndef  RECONSTRUCTEDPARTICLE2TRACK_ANALYZERS_H
#define  RECONSTRUCTEDPARTICLE2TRACK_ANALYZERS_H

#include <cmath>
#include <vector>

#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/TrackState.h"
#include "TLorentzVector.h"




/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2TRK_D0 (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<edm4hep::TrackState> tracks);

/// Return the Z0 of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2TRK_Z0 (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<edm4hep::TrackState> tracks);

/// Return the Phi of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2TRK_phi (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<edm4hep::TrackState> tracks);

/// Return the omega of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2TRK_omega (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<edm4hep::TrackState> tracks);

/// Return the tanLambda of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2TRK_tanLambda (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<edm4hep::TrackState> tracks);


/*
struct D0s {
   D0s() ;
   std::vector<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> KPi_from_D0,
      ROOT::VecOps::RVec<float> KPi_from_D0_track_phi,
              ROOT::VecOps::RVec<float> KPi_from_D0_track_tanLambda,
              ROOT::VecOps::RVec<float> KPi_from_D0_track_omega) ;
};
*/

#endif
