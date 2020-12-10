#include "ReconstructedParticle2Track.h"


ROOT::VecOps::RVec<float> getRP2TRK_D0(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,  ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).D0);
    else result.push_back(std::nan(""));
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2TRK_Z0(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,  ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).Z0);
    else result.push_back(std::nan(""));
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2TRK_phi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,  ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).phi);
    else result.push_back(std::nan(""));
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2TRK_omega(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,  ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).omega);
    else result.push_back(std::nan(""));
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2TRK_tanLambda(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,  ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).tanLambda);
    else result.push_back(std::nan(""));
  }
  return result;
}

/*
D0s::D0s() {}
std::vector<edm4hep::ReconstructedParticleData> D0s::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> KPi_from_D0,
              ROOT::VecOps::RVec<float> KPi_from_D0_track_phi,
              ROOT::VecOps::RVec<float> KPi_from_D0_track_tanLambda,
              ROOT::VecOps::RVec<float> KPi_from_D0_track_omega) {

///  this was  for  a test
         
               
  std::vector<edm4hep::ReconstructedParticleData> result;
         
  int n = KPi_from_D0.size();
  int nJ = n / 2;
  
  float  Bfield = 2. ;
  float mass1 = 0.494 ;  // K mass
  float mass2 = 0.140 ;  // Pi mass

  for (int i=0; i < nJ; i++) {
     int i1 = i*2 ;
     int i2 = i1 + 1;
     TLorentzVector reso_lv;
     edm4hep::ReconstructedParticleData reso;
     reso.charge = KPi_from_D0[i1].charge + KPi_from_D0[i2].charge ;
  
     float pt1 = std::abs( 1e-3 * 0.3 * Bfield / KPi_from_D0_track_omega[i1] );   /// omega is in mm-1
     //tan(theta) = 1 / tan(lambda); eta = -log( tan ( theta/2) ) );
     float eta1 = - log( std::abs( tan ( 0.5 / KPi_from_D0_track_tanLambda[i1] )) );   // eta is not correct
     float phi1 = KPi_from_D0_track_phi[i1];
     TLorentzVector leg1 ;
     leg1.SetPtEtaPhiM ( pt1, eta1, phi1, mass1 );
     std::cout << " a reco K leg  pt = " << pt1 << " eta = " << eta1 << " phi = " << phi1 << std::endl;
     
     TLorentzVector leg2 ;
     float pt2 = std::abs( 1e-3 * 0.3 * Bfield / KPi_from_D0_track_omega[i2] );
     float eta2 = - log( std::abs( tan ( 0.5 / KPi_from_D0_track_tanLambda[i2] ) ));
     float phi2 = KPi_from_D0_track_phi[i2];
     leg2.SetPtEtaPhiM ( pt2, eta2, phi2, mass2 );
     std::cout << " a reco Pi leg  pt = " << pt2 << " eta = " << eta2 << " phi = " << phi2 << std::endl;

     
     reso_lv = leg1 + leg2;
      reso.momentum.x = reso_lv.Px();
      reso.momentum.y = reso_lv.Py();
      reso.momentum.z = reso_lv.Pz();
      reso.mass = reso_lv.M();
      result.emplace_back(reso);
      
  }  // end loop over the JPsis
  return result;
      
}     

*/
