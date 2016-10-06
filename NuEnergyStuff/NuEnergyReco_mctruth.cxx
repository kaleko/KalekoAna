#ifndef NUENERGYRECO_MCTRUTH_CXX
#define NUENERGYRECO_MCTRUTH_CXX

#include "NuEnergyReco_mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"

#include <algorithm> //for std::sort

namespace larlite {

  bool NuEnergyReco_mctruth::initialize() {

    if(tree) delete tree;
    tree = new TTree("tree","tree");
    tree->Branch("true_nu_E",&true_nu_E,"true_nu_E/D");
    tree->Branch("nu_pdg",&nu_pdg,"nu_pdg/I");
    tree->Branch("nuance_intxn",&nuance_intxn,"nuance_intxn/I");
    tree->Branch("isCC",&isCC,"isCC/O");
    tree->Branch("mode",&_mode,"mode/I");
    tree->Branch("areOtherParticles",&areOtherParticles,"areOtherParticles/O");
    tree->Branch("otherParticles_KE",&otherParticles_KE,"otherParticles_KE/D");
    tree->Branch("E_ccqe",&E_ccqe,"E_ccqe/D");
    tree->Branch("lep_E",&lep_E,"lep_E/D");
    tree->Branch("n_lep",&n_lep,"n_lep/I");
    tree->Branch("n_protons",&n_protons,"n_protons/I");
    tree->Branch("tot_prot_KE",&tot_prot_KE,"tot_prot_KE/D");
    tree->Branch("n_neutrons",&n_neutrons,"n_neutrons/I");
    tree->Branch("tot_neutron_KE",&tot_neutron_KE,"tot_neutron_KE/D");
    tree->Branch("max_prot_KE",&max_prot_KE,"max_prot_KE/D");
    tree->Branch("second_prot_KE",&second_prot_KE,"second_prot_KE/D");
    tree->Branch("tot_pt",&tot_pt,"tot_pt/D");

    return true;
  }
  
  bool NuEnergyReco_mctruth::analyze(storage_manager* storage) {

    ResetTTreeVars();
    
    auto ev_tru = storage->get_data<event_mctruth >("generator");

    if(!ev_tru){
      std::cout<<"wtf"<<std::endl;
      return false;
    }
    
    auto mcnu_traj = ev_tru->at(0).GetNeutrino().Nu().Trajectory();

    true_nu_E = 1000.*mcnu_traj.at(0).E();
    nu_pdg = ev_tru->at(0).GetNeutrino().Nu().PdgCode();
    nuance_intxn = ev_tru->at(0).GetNeutrino().InteractionType();
    //CCNC returns 0 if CC, 1 if NC
    isCC = ev_tru->at(0).GetNeutrino().CCNC() == 0 ? true : false;
    _mode = ev_tru->at(0).GetNeutrino().Mode();
    n_lep = 0;
    std::vector<float> prot_energies;
    prot_energies.clear();
    std::vector<float> neut_energies;
    neut_energies.clear();

    
    auto mcpart_list = ev_tru->at(0).GetParticles();
    TLorentzVector total_momentum(0.,0.,0.,0.);
    for (auto const mcp : mcpart_list){
   
      if (mcp.StatusCode() == 1){

	double part_E  = 1000.*mcp.Trajectory().at(0).E();
	double part_KE = 1000.*(mcp.Trajectory().at(0).E()-mcp.Trajectory().at(0).Momentum().M());
	auto part_3mom = 1000.*mcp.Trajectory().at(0).Momentum().Vect();
	std::vector<double> part_4mom = {part_3mom.X(),part_3mom.Y(),part_3mom.Z(),part_E};

	//Store the mcparticle lepton energy while we're at it
	if(abs(mcp.PdgCode()) == 13 || abs(mcp.PdgCode()) == 11){
	  lep_E = part_E;
	  n_lep++;
	  //	  E_ccqe = 1000.*_calc.ComputeECCQE(part_4mom);
	  total_momentum += mcp.Trajectory().at(0).Momentum();
	}
	
	if(abs(mcp.PdgCode()) == 2212 && part_KE > _min_p_E_MEV){
	  total_momentum += mcp.Trajectory().at(0).Momentum();
	  prot_energies.push_back(part_KE);
	}

	if(abs(mcp.PdgCode()) == 2112 && part_KE > _min_n_E_MEV) 
	  neut_energies.push_back(part_KE);
   
	//Sometimes there are 311 (K0), 3222 (sigma+), 3122 (lambda+)
	//particles in the event...
	if (abs(mcp.PdgCode()) != 2212 && //proton
	    abs(mcp.PdgCode()) != 2112 && //neutron
	    abs(mcp.PdgCode()) != 11  &&  //electron
	    abs(mcp.PdgCode()) != 12  &&  //nue
	    abs(mcp.PdgCode()) != 2000000101 ){//bindino
	  areOtherParticles = true; 
	  total_momentum += mcp.Trajectory().at(0).Momentum();
	  otherParticles_KE += part_KE;
	}
	
      }
    }
    
    //this equals x component^2 + y component ^2, sqrt
    tot_pt = total_momentum.Vect().Perp();

    n_protons = prot_energies.size();
    n_neutrons = neut_energies.size();

    //sort energies vector from lowest to highest
    std::sort(prot_energies.begin(),prot_energies.end());
    std::sort(neut_energies.begin(),neut_energies.end());

    max_prot_KE = n_protons > 0 ? prot_energies.back() : -999.;
    second_prot_KE = n_protons > 1 ? prot_energies.end()[-2] : -999.;

    tot_prot_KE = 0.;
    tot_neutron_KE = 0.;

    for (double E : prot_energies)
      tot_prot_KE += E;

    for (double E : neut_energies)
      tot_neutron_KE += E;

    tree->Fill();
  
    return true;
  }

  bool NuEnergyReco_mctruth::finalize() {

    if(_fout) { _fout->cd(); tree->Write(); }
  
    return true;
  }

  void NuEnergyReco_mctruth::ResetTTreeVars(){
    true_nu_E = 0;
    nu_pdg = -1;    
    nuance_intxn = -1;
    isCC = false;
    areOtherParticles = false;
    otherParticles_KE = 0;
    _mode = -1;
    E_ccqe = 0;
    lep_E = 0;
    tot_pt = -1.;
  }
}
#endif
