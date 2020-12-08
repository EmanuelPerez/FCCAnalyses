#include "ReconstructedParticle2MC.h"

ROOT::VecOps::RVec<float> getRP2MC_p(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);

  for (unsigned int i=0; i<recind.size();i++) {
    TLorentzVector tlv;
    tlv.SetXYZM(mc.at(mcind.at(i)).momentum.x,mc.at(mcind.at(i)).momentum.y,mc.at(mcind.at(i)).momentum.z,mc.at(mcind.at(i)).mass);
    result[recind.at(i)]=tlv.P();
  }
  return result;
}

ROOT::VecOps::RVec<TLorentzVector> getRP2MC_tlv(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<TLorentzVector> result;
  result.resize(reco.size(),TLorentzVector());

  for (unsigned int i=0; i<recind.size();i++) {
    TLorentzVector tlv;
    tlv.SetXYZM(mc.at(mcind.at(i)).momentum.x,mc.at(mcind.at(i)).momentum.y,mc.at(mcind.at(i)).momentum.z,mc.at(mcind.at(i)).mass);
    result[recind.at(i)]=tlv;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_px(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).momentum.x;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_py(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).momentum.y;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_pz(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).momentum.z;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_pdg(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).PDG;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_charge(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).charge;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_mass(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).mass;
  }
  return result;
}

ROOT::VecOps::RVec<int> getRP2MC_index(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco) {
  ROOT::VecOps::RVec<int> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mcind.at(i);
  }
  return result;
}

ROOT::VecOps::RVec<int> getRP2MC_parentid (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> parents){
  ROOT::VecOps::RVec<int> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    if (mc.at(mcind.at(i)).parents_begin!=mc.at(mcind.at(i)).parents_end){
      result[recind.at(i)]=parents.at(mc.at(mcind.at(i)).parents_begin);
    }
  }

  /*  if (recind.size()>reco.size()){ 
    std::cout << recind.size() <<"========="<<reco.size()<<std::endl;
    for (unsigned int i=0; i<recind.size();i++) {
      if (i<recind.size()-1 && recind[i]==recind[i+1]){
	
	TLorentzVector tlv;
	tlv.SetXYZM(mc.at(mcind.at(i)).momentum.x,mc.at(mcind.at(i)).momentum.y,mc.at(mcind.at(i)).momentum.z,mc.at(mcind.at(i)).mass);
	TLorentzVector tlv2;
	tlv2.SetXYZM(reco.at(recind.at(i)).momentum.x,reco.at(recind.at(i)).momentum.y,reco.at(recind.at(i)).momentum.z,reco.at(recind.at(i)).mass);
	std::cout << "n mc " << mc.size() << " rec ind " << recind.at(i) << " reco P "<< tlv2.P()<< "  mc ind " << mcind.at(i) << " truth P " << tlv.P() << " pdg_id " << mc.at(mcind.at(i)).PDG  << "  parent id " << parents.at(mc.at(mcind.at(i)).parents_begin) << " parent pdg id " << mc.at(parents.at(mc.at(mcind.at(i)).parents_begin)).PDG << std::endl;

	tlv.SetXYZM(mc.at(mcind.at(i+1)).momentum.x,mc.at(mcind.at(i+1)).momentum.y,mc.at(mcind.at(i+1)).momentum.z,mc.at(mcind.at(i+1)).mass);
	tlv2.SetXYZM(reco.at(recind.at(i+1)).momentum.x,reco.at(recind.at(i+1)).momentum.y,reco.at(recind.at(i+1)).momentum.z,reco.at(recind.at(i+1)).mass);
	std::cout << "n mc " << mc.size() << " rec ind " << recind.at(i+1) << " reco P "<< tlv2.P()<< "  mc ind " << mcind.at(i+1) << " truth P " << tlv.P() << " pdg_id " << mc.at(mcind.at(i+1)).PDG  << "  parent id " << parents.at(mc.at(mcind.at(i+1)).parents_begin) << " parent pdg id " << mc.at(parents.at(mc.at(mcind.at(i+1)).parents_begin)).PDG << std::endl;
	}
    }
    }*/
  return result;
}


ROOT::VecOps::RVec<float>  getRP2MC_p_func::operator() (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);

  for (unsigned int i=0; i<recind.size();i++) {
    TLorentzVector tlv;
    tlv.SetXYZM(mc.at(mcind.at(i)).momentum.x,mc.at(mcind.at(i)).momentum.y,mc.at(mcind.at(i)).momentum.z,mc.at(mcind.at(i)).mass);
    result[recind.at(i)]=tlv.P();
  }

  if (recind.size()>reco.size()){ 
    std::cout << recind.size() <<"========="<<reco.size()<<std::endl;
     for (unsigned int i=0; i<recind.size();i++) {
       TLorentzVector tlv;
       tlv.SetXYZM(mc.at(mcind.at(i)).momentum.x,mc.at(mcind.at(i)).momentum.y,mc.at(mcind.at(i)).momentum.z,mc.at(mcind.at(i)).mass);
       TLorentzVector tlv2;
       tlv2.SetXYZM(reco.at(recind.at(i)).momentum.x,reco.at(recind.at(i)).momentum.y,reco.at(recind.at(i)).momentum.z,reco.at(recind.at(i)).mass);
       std::cout << "n mc " << mc.size() << " rec ind " << recind.at(i) << " reco P "<< tlv2.P()<< "  mc ind " << mcind.at(i) << " truth P " << tlv.P() << " pdg_id " << mc.at(mcind.at(i)).PDG << " parent_begin " <<  mc.at(mcind.at(i)).parents_begin << " parent_end " <<  mc.at(mcind.at(i)).parents_end << " daut_begin " <<  mc.at(mcind.at(i)).daughters_begin << " daut_end " <<  mc.at(mcind.at(i)).daughters_end <<std::endl;
     }
  }
  return result;
}




// -------------------------------------------------------------------------------------------------

// -- select RecoParticles associated with MC muons
// -- ( for muons from JPsi, can not use the Muon collection because it oontains
// -- only the isolated muons)

selRP_PDG::selRP_PDG( int arg_pdg ): m_PDG(arg_pdg) {} ;
std::vector<edm4hep::ReconstructedParticleData>  selRP_PDG::operator() (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {

  std::vector<edm4hep::ReconstructedParticleData> result;

  for (int i=0; i<recind.size();i++) {
      int reco_idx = recind.at(i);
      int mc_idx = mcind.at(i);
      int pdg = mc.at(mc_idx).PDG ;
      if ( std::abs( pdg ) == std::abs( m_PDG)  ) {
         result.push_back( reco.at( reco_idx ) ) ;
      }
  }
  return result;
}


// -- select the reco'ed particles associated with MC muons that come from the
// -- decay of a J/Psi

selMuons_JPsimatch::selMuons_JPsimatch( int arg_dum ) : m_dummy(arg_dum) {};
    
std::vector<edm4hep::ReconstructedParticleData> selMuons_JPsimatch::operator() (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> mcdaughters) {
    
  std::vector<edm4hep::ReconstructedParticleData> result;
  std::vector< std::array<int, 2> >  MCmuons_from_JPsis = get_MC_muons_from_JPsis( mc, mcdaughters) ;
  int nJPsis = MCmuons_from_JPsis.size() ;

  if ( nJPsis <1) return result ;


  for( int ijpsi=0; ijpsi < nJPsis; ijpsi ++) {

      std::array<int, 2>  MCmuons_JPsi = MCmuons_from_JPsis[ ijpsi ] ;
      if (MCmuons_JPsi[0]<0 || MCmuons_JPsi[1]<0) continue;

      int nlegs =0;
      for (int i=0; i<recind.size();i++) {
          int reco_idx = recind.at(i);
          int mc_idx = mcind.at(i);
          if ( mc_idx == MCmuons_JPsi[0] || mc_idx == MCmuons_JPsi[1]) {
             result.push_back( reco.at( reco_idx ) ) ;	 
             nlegs ++;
          }
      }

      if ( nlegs == 1 ) {   // one of the muons from this JPsi didnot make a RecoParticle
			    // e.g. outside of the tracker acceptance
			    // in which case, remove the other leg if it was found
	result.pop_back() ;
      }

  }  // loop over the JPsis

  return result;
}




