
#ifndef  ALGORITHM_ANALYZERS_H
#define  ALGORITHM_ANALYZERS_H

#include <cmath>
#include <vector>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include "edm4hep/Vector3f.h"
#include "edm4hep/Vector3d.h"
#include "edm4hep/Vector2i.h"
#include "TFitter.h"


struct getRP_combination{
  getRP_combination(int arg_n, int arg_charge, bool arg_abs);
  int  m_n;
  int  m_charge;
  bool m_abs;
  ROOT::VecOps::RVec<int> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
};


struct sphericityFit {
  sphericityFit(ROOT::VecOps::RVec<float> arg_px, ROOT::VecOps::RVec<float> arg_py, ROOT::VecOps::RVec<float> arg_pz);
  ROOT::VecOps::RVec<float> m_px;
  ROOT::VecOps::RVec<float> m_py;
  ROOT::VecOps::RVec<float> m_pz;
  float operator()(const double *par);
};


struct minimize_sphericity {
  minimize_sphericity(std::string arg_minname, std::string arg_algoname);
  char const *m_minname  = "Minuit2";
  char const *m_algoname = "";
  ROOT::VecOps::RVec<float> operator()(ROOT::VecOps::RVec<float> px, ROOT::VecOps::RVec<float> py, ROOT::VecOps::RVec<float> pz);
};



struct thrustFit {
  thrustFit(ROOT::VecOps::RVec<float> arg_px, ROOT::VecOps::RVec<float> arg_py, ROOT::VecOps::RVec<float> arg_pz);
  ROOT::VecOps::RVec<float> m_px;
  ROOT::VecOps::RVec<float> m_py;
  ROOT::VecOps::RVec<float> m_pz;
  float operator()(const double *par);
};


struct minimize_thrust {
  minimize_thrust(std::string arg_minname, std::string arg_algoname);
  char const *m_minname  = "Minuit2";
  char const *m_algoname = "";
  ROOT::VecOps::RVec<float> operator()(ROOT::VecOps::RVec<float> px, ROOT::VecOps::RVec<float> py, ROOT::VecOps::RVec<float> pz);
};



ROOT::VecOps::RVec<float> axisCosTheta(ROOT::VecOps::RVec<float> axis, ROOT::VecOps::RVec<float> px, ROOT::VecOps::RVec<float> py, ROOT::VecOps::RVec<float> pz);

struct getAxisCharge {
  getAxisCharge(bool arg_pos);
  bool m_pos = 0;
  float operator() (ROOT::VecOps::RVec<float> angle, ROOT::VecOps::RVec<float> charge, ROOT::VecOps::RVec<float> px, ROOT::VecOps::RVec<float> py, ROOT::VecOps::RVec<float> pz);
};

struct getAxisMass {
  getAxisMass(bool arg_pos);
  bool m_pos = 0;
  float operator() (ROOT::VecOps::RVec<float> angle, ROOT::VecOps::RVec<float> energy, ROOT::VecOps::RVec<float> px, ROOT::VecOps::RVec<float> py, ROOT::VecOps::RVec<float> pz);
};

/// Get the mass from a list of reconstructed particles
float getMass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

#endif
