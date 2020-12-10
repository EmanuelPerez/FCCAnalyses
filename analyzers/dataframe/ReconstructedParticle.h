
#ifndef  RECONSTRUCTEDPARTICLE_ANALYZERS_H
#define  RECONSTRUCTEDPARTICLE_ANALYZERS_H

#include <cmath>
#include <vector>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"

#include <random>
#include <chrono>

struct JPsis {
   JPsis() ;
   std::vector<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> muons_from_JPsis) ;
};

/*
struct Dimuons {
  Dimuons();
  std::vector<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs);
};
*/

struct Pairs {
  Pairs( bool same) ;
  float m_same;
  std::vector<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs1,
						             ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs2);
};

/// TO BE MOVED LATER
struct ResonanceBuilder {
  int m_resonance_pdgid;
  float m_resonance_mass;
  ResonanceBuilder(int arg_resonance_pdgid, float arg_resonance_mass);
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs);
};

struct recoil {
  recoil(float arg_sqrts);
  float m_sqrts = 240.0;
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) ;
};

/// select ReconstructedParticles with transverse momentum greater than a minimum value [GeV]
struct selRP_pT {
  selRP_pT(float arg_min_pt);
  float m_min_pt = 20; //> transverse momentum threshold [GeV]
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
};

struct selRP_E {
  selRP_E(float arg_min_e);
  float m_min_e = 2; //> energy threshold [GeV]
  std::vector<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
};


struct  selRP_mass {
  selRP_mass( float arg_mass_min, float arg_mass_max) ;
  float m_mass_min = 0. ;
  float m_mass_max = 10. ;
  std::vector<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
};

struct selRP_Fakes {
  selRP_Fakes( float arg_fakeRate, float  arg_mass );
  float m_fakeRate = 1e-3;
  float m_mass = 0.106;  // muon mass
  std::default_random_engine m_generator;
  std::uniform_real_distribution<float> m_flat;
  std::vector<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
};

/// return reconstructed particles
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> getRP(ROOT::VecOps::RVec<int> index, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the transverse momenta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_pt(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the momenta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_p(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the momenta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_px(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the momenta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_py(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the momenta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_pz(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the pseudo-rapidity of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_eta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the rapidity of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_y(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the theta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_theta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the phi of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_phi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the energy of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_e(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the masses of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_mass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in); 

/// return the charges of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_charge(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in); 

/// return the TlorentzVector of the input ReconstructedParticles
ROOT::VecOps::RVec<TLorentzVector> getRP_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// concatenate both input vectors and return the resulting vector
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> mergeParticles(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> y);

std::vector<edm4hep::ReconstructedParticleData> my_mergeParticles(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> y);

/// return the size of the input collection
int getRP_n(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

#endif
