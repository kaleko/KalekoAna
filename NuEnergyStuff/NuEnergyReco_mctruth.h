/**
 * \file NuEnergyReco_mctruth.h
 *
 * \ingroup scratch_ana
 * 
 * \brief Class def header for a class NuEnergyReco_mctruth
 *
 * @author davidkaleko
 */

/** \addtogroup scratch_ana

    @{*/

#ifndef LARLITE_NUENERGYRECO_MCTRUTH_H
#define LARLITE_NUENERGYRECO_MCTRUTH_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcpart.h"
#include "TTree.h"
//#include "/Users/davidkaleko/larlite/UserDev/LowEnergyExcess/Utilities/ECCQECalculator.h"

namespace larlite {
  /**
     \class NuEnergyReco_mctruth
     User custom analysis class made by davidkaleko
     This uses on mctruth to study the best way to compute the true
     neutrino energy from "reconstructed" (mctruth) quantities.
     It basically makes a ttree and stores it for TTree->Draw stuff.
   */
  class NuEnergyReco_mctruth : public ana_base{
  
  public:

    /// Default constructor
    NuEnergyReco_mctruth(){ 
      _name="NuEnergyReco_mctruth"; 
      _min_p_E_MEV = 30.;
      _min_n_E_MEV = -1.;
      tree=0; 
      _fout=0;
    };
    
    /// Default destructor
    virtual ~NuEnergyReco_mctruth(){};

    /** IMPLEMENT in NuEnergyReco_mctruth.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in NuEnergyReco_mctruth.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in NuEnergyReco_mctruth.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    
    void SetMinPE_MEV(double kaleko){ _min_p_E_MEV = kaleko; }
    void SetMinNE_MEV(double kaleko){ _min_n_E_MEV = kaleko; }

    protected:

    
    double _min_p_E_MEV;
    double _min_n_E_MEV;
    
    void ResetTTreeVars();

    TTree* tree;
    double true_nu_E;
    int nu_pdg;
    int nuance_intxn;
    bool isCC;
    int _mode;
    bool areOtherParticles;
    double otherParticles_KE;
    double E_ccqe;
    double lep_E;
    int n_lep;
    double tot_prot_KE;
    int n_protons;
    double tot_neutron_KE;
    int n_neutrons;
    double max_prot_KE;
    double second_prot_KE;
    double tot_pt;

    //    larlite::util::ECCQECalculator _calc;


  };
}
#endif
