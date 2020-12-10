#include "ReconstructedParticle.h"

// -- Build dimuon objects from a collection of RecoParticles
//Dimuons::Dimuons() { }
//std::vector<edm4hep::ReconstructedParticleData> Dimuons::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs) {

Pairs::Pairs( bool same) { 
   m_same = same;
}

std::vector<edm4hep::ReconstructedParticleData> Pairs::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs1,
                                                                    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs2) {
  std::vector<edm4hep::ReconstructedParticleData> result;
  int n1 = legs1.size();
  int n2 = legs2.size();
  //std::cout << " in pairs  : n1 = " << n1 << " n2 = " << n2 << std::endl;

  int n1_min = 1;
  if ( m_same) n1_min = 2;
  if (n1 >= n1_min ) {
/*
    ROOT::VecOps::RVec<bool> v(n);
    std::fill(v.end() - 2, v.end(), true);
    do {
      edm4hep::ReconstructedParticleData reso;
      //reso.pdg = m_resonance_pdgid;
      TLorentzVector reso_lv;
      for (int i = 0; i < n; ++i) {
          if (v[i]) {
            reso.charge += legs[i].charge;
            TLorentzVector leg_lv;
            leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
            reso_lv += leg_lv;
          }
      }
      reso.momentum.x = reso_lv.Px();
      reso.momentum.y = reso_lv.Py();
      reso.momentum.z = reso_lv.Pz();
      reso.mass = reso_lv.M();
      result.emplace_back(reso);
    } while (std::next_permutation(v.begin(), v.end()));
*/
   // float muon_mass = 0.106 ; 
	// the muon mass is hardcoded below instead of using the mass of the recoed particles,
	// because this is  used too  for fake (Pi, K's) muons.

   edm4hep::ReconstructedParticleData reso;
   TLorentzVector reso_lv;
   for (int i=0; i < n1; i++) {
      TLorentzVector leg1 ;
      //leg1.SetXYZM( legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
      leg1.SetXYZM( legs1[i].momentum.x, legs1[i].momentum.y, legs1[i].momentum.z, legs1[i].mass);
      int  jmin = i+1;
      int  jmax = n1;
      if ( ! m_same ) {
         jmin = 0; 
         jmax = n2 ;
      }
      for (int j=jmin; j < jmax; j++) {
         TLorentzVector leg2;
         //leg2.SetXYZM( legs[j].momentum.x, legs[j].momentum.y, legs[j].momentum.z, legs[j].mass);
         leg2.SetXYZM( legs2[j].momentum.x, legs2[j].momentum.y, legs2[j].momentum.z, legs2[j].mass );
         reso.charge = legs1[i].charge + legs2[j].charge ;
         reso_lv = leg1+ leg2;
         reso.momentum.x = reso_lv.Px();
         reso.momentum.y = reso_lv.Py();
         reso.momentum.z = reso_lv.Pz();
         reso.mass = reso_lv.M();
         result.emplace_back(reso);
         //std::cout << " Leg1: mass = " << legs1[i].mass << " charge = " << legs1[i].charge << std::endl;
         //std::cout << " Leg2: mass = " << legs2[j].mass << " charge = " << legs2[j].charge << std::endl;
      }
   }
  }
  return result;
}


JPsis::JPsis() { } 
std::vector<edm4hep::ReconstructedParticleData> JPsis::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> muons_from_JPsis) {

// input = the muons from selMuons_JPsimatch
// by construction, pairs shouldbe made from muons 0 and 1, (2 and 3) etc,in case > 1 JPsi

  std::vector<edm4hep::ReconstructedParticleData> result;

  int n = muons_from_JPsis.size();
  int nJ = n /2;	// number of JPsis to be built

  for (int i=0; i < nJ; i++) {
     int i1 = i*2 ;
     int i2 = i1 + 1;
     TLorentzVector reso_lv;
     edm4hep::ReconstructedParticleData reso;
     reso.charge = muons_from_JPsis[i1].charge + muons_from_JPsis[i2].charge ;
     TLorentzVector leg1 ;
     leg1.SetXYZM( muons_from_JPsis[i1].momentum.x, muons_from_JPsis[i1].momentum.y, muons_from_JPsis[i1].momentum.z, muons_from_JPsis[i1].mass);
     TLorentzVector leg2 ;
     leg2.SetXYZM( muons_from_JPsis[i2].momentum.x, muons_from_JPsis[i2].momentum.y, muons_from_JPsis[i2].momentum.z, muons_from_JPsis[i2].mass);
     reso_lv = leg1 + leg2;
     ///std::cout <<" masses of the JPsi legs = " << muons_from_JPsis[i1].mass << " " << muons_from_JPsis[i2].mass << std::endl;

     //std::cout << " in JPsis: leg1 mass = " << muons_from_JPsis[i1].mass << " charge " <<  muons_from_JPsis[i1].charge  << std::endl;
     //std::cout << " in JPsis: leg2 mass = " << muons_from_JPsis[i2].mass << " charge " <<  muons_from_JPsis[i2].charge  << std::endl;

      reso.momentum.x = reso_lv.Px();
      reso.momentum.y = reso_lv.Py();
      reso.momentum.z = reso_lv.Pz();
      reso.mass = reso_lv.M();
      result.emplace_back(reso);

  }  // end loop over the JPsis
  return result;

}



//TOBEMOVED LATER
ResonanceBuilder::ResonanceBuilder(int arg_resonance_pdgid, float arg_resonance_mass) {m_resonance_pdgid = arg_resonance_pdgid; m_resonance_mass = arg_resonance_mass;}
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> ResonanceBuilder::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  int n = legs.size();
  if (n >1) {
    ROOT::VecOps::RVec<bool> v(n);
    std::fill(v.end() - 2, v.end(), true);
    do {
      edm4hep::ReconstructedParticleData reso;
      //reso.pdg = m_resonance_pdgid;
      TLorentzVector reso_lv; 
      for (int i = 0; i < n; ++i) {
          if (v[i]) {
            reso.charge += legs[i].charge;
            TLorentzVector leg_lv;
            leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
            reso_lv += leg_lv;
          }
      }
      reso.momentum.x = reso_lv.Px();
      reso.momentum.y = reso_lv.Py();
      reso.momentum.z = reso_lv.Pz();
      reso.mass = reso_lv.M();
      result.emplace_back(reso);
    } while (std::next_permutation(v.begin(), v.end()));
  }
  if (result.size() > 1) {
    auto resonancesort = [&] (edm4hep::ReconstructedParticleData i ,edm4hep::ReconstructedParticleData j) { return (abs( m_resonance_mass -i.mass)<abs(m_resonance_mass-j.mass)); };
    std::sort(result.begin(), result.end(), resonancesort);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator first = result.begin();
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator last = result.begin() + 1;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> onlyBestReso(first, last);
    return onlyBestReso;
  } else {
    return result;
  }
}



recoil::recoil(float arg_sqrts) : m_sqrts(arg_sqrts) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  recoil::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  auto recoil_p4 = TLorentzVector(0, 0, 0, m_sqrts);
  for (auto & v1: in) {
    TLorentzVector tv1;
    tv1.SetXYZM(v1.momentum.x, v1.momentum.y, v1.momentum.z, v1.mass);
    recoil_p4 -= tv1;
  }
  auto recoil_fcc = edm4hep::ReconstructedParticleData();
  recoil_fcc.momentum.x = recoil_p4.Px();
  recoil_fcc.momentum.y = recoil_p4.Py();
  recoil_fcc.momentum.z = recoil_p4.Pz();
  recoil_fcc.mass = recoil_p4.M();
  result.push_back(recoil_fcc);
  return result;
};

ROOT::VecOps::RVec<float> getRP_pt(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
 ROOT::VecOps::RVec<float> result;
 for (size_t i = 0; i < in.size(); ++i) {
   result.push_back(sqrt(in[i].momentum.x * in[i].momentum.x + in[i].momentum.y * in[i].momentum.y));
 }
 return result;
}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> mergeParticles(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> y) {
  //to be keept as ROOT::VecOps::RVec
  std::vector<edm4hep::ReconstructedParticleData> result;
  result.reserve(x.size() + y.size());
  result.insert( result.end(), x.begin(), x.end() );
  result.insert( result.end(), y.begin(), y.end() );
  return ROOT::VecOps::RVec(result);
}


// wih the mergeParticles above, no  collection  is written to the ntuple ??
std::vector<edm4hep::ReconstructedParticleData> my_mergeParticles(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> y) {
      std::vector<edm4hep::ReconstructedParticleData> result;
  //std::cout << " enter in my_mergeParticles sizes = " << x.size() << " " << y.size() <<  std::endl;
  result.reserve(x.size() + y.size());
  result.insert( result.end(), x.begin(), x.end() );
  result.insert( result.end(), y.begin(), y.end() );
  //std::cout << " size of the merged collection " << result.size() <<  std::endl;
  return result;
}


ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> getRP(ROOT::VecOps::RVec<int> index, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  for (size_t i = 0; i < index.size(); ++i) {
    if (index[i]>-1)
      result.push_back(in.at(index[i]));
    //else
    //  std::cout << "electron index negative " << index[i]<<std::endl;
  }  
  return result;
}


ROOT::VecOps::RVec<float> getRP_mass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.mass);
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_eta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Eta());
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_phi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Phi());
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_e(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.energy);
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_p(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.P());
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_px(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.x);
  }
  return result;
}


ROOT::VecOps::RVec<float> getRP_py(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.y);
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_pz(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.z);
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_charge(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.charge);
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_y(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Rapidity());
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_theta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Theta());
  }
  return result;
}

ROOT::VecOps::RVec<TLorentzVector> getRP_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<TLorentzVector> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv);
  }
  return result;
}


int getRP_n(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x) {
  int result =  x.size();
  return result;
}

selRP_pT::selRP_pT(float arg_min_pt) : m_min_pt(arg_min_pt) {};

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  selRP_pT::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if (std::sqrt(std::pow(p.momentum.x,2) + std::pow(p.momentum.y,2)) > m_min_pt) {
      result.emplace_back(p);
    }
  }
  return result;
}


selRP_E::selRP_E(float arg_min_e) : m_min_e(arg_min_e) {};

std::vector<edm4hep::ReconstructedParticleData>  selRP_E::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  std::vector<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if ( p.energy > m_min_e) {
      result.emplace_back(p);
    }
  }
  return result;
}

selRP_mass::selRP_mass( float arg_mass_min, float  arg_mass_max) : m_mass_min( arg_mass_min ), m_mass_max( arg_mass_max ) { } ;
std::vector<edm4hep::ReconstructedParticleData> selRP_mass::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  std::vector<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if ( p.mass > m_mass_min  && p.mass < m_mass_max) {
      result.emplace_back(p);
    }
  }
  return result;
}


// -------------------------------------------------------------------------------------------------

// -- Fakes, more general :
// -- Randomly select RecoParticles according to a given fake rate (misid efficiency)

selRP_Fakes::selRP_Fakes( float arg_fakeRate, float  arg_mass ) : m_fakeRate(arg_fakeRate), m_mass( arg_mass)  {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  m_generator = generator;
  std::uniform_real_distribution<float> flatdis(0.,1.);
  m_flat.param( flatdis.param() );
};

std::vector<edm4hep::ReconstructedParticleData> selRP_Fakes::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {

  std::vector<edm4hep::ReconstructedParticleData> result;

  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    float arandom =  m_flat (m_generator );
    if ( arandom <= m_fakeRate) {
       edm4hep::ReconstructedParticleData reso;
       reso.momentum.x = p.momentum.x ;
       reso.momentum.y = p.momentum.y ;
       reso.momentum.z = p.momentum.z ;
       reso.mass = m_mass;
       reso.charge = p.charge;
       result.push_back( reso );
    }
  }

  return result;
}

