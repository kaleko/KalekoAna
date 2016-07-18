/**
 * \file TestMCTruth.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class TestMCTruth
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_TESTMCTRUTH_H
#define LARLITE_TESTMCTRUTH_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAABox.h"
#include "FidVolBox.h"
#include "DataFormat/mctruth.h"
#include "TH2.h"
#include "TH1.h"
#include "TTree.h"
#include "GeoAlgo/GeoVector.h"
#include "NuEnergyCalc.h"

namespace larlite {
    /**
       \class TestMCTruth
       User custom analysis class made by SHELL_USER_NAME
     */
    class TestMCTruth : public ana_base {

    public:

        /// Default constructor
        TestMCTruth() { _name = "TestMCTruth"; _fout = 0; _h_mup = 0; _h_mup_diff = 0; _tree = 0; _Ecalc = 0;}

        /// Default destructor
        virtual ~TestMCTruth() {}

        virtual bool initialize();

        virtual bool analyze(storage_manager* storage);

        virtual bool finalize();

        void printTruthStuff(const event_mctruth &mctruth);
        
    protected:

        // Fiducial volume box
        geoalgo::AABox _fidvolBox;

        size_t evt_counter;
        size_t kept_events;

        TH2D* _h_mup;
        TH1F* _h_mup_diff;

        TTree *_tree;
        double absphidiff_minuspi;
        double reco_mu_p;
        double true_mu_p;
        double true_p_p;
        double premnant;
        double ptremnant;
        double Eremnant;
        double plane_Z_anglediff;
        double true_nu_E;
        double Eccqe_truemumom;
        double Eccqe_recomumom;
        double true_mu_TotE;
        double true_p_TotE;
        double true_p_KE;
        double true_n_TotE;

        double reco_mu_TotE;
        double reco_nu_E_4momentum;
        double reco_mu_TotE_quadratic;
        double reco_mu_p_quadratic;
        double thetap;
        double thetamu;
        double thetamup;

        double bindinoE;
        double bindinop;
        bool genieCCQE;
        double genieMomDiff;
        NuEnergyCalc *_Ecalc;
    };
}
#endif

//**************************************************************************
//
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group
