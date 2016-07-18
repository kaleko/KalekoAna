#ifndef LARLITE_TESTMCTRUTH_CXX
#define LARLITE_TESTMCTRUTH_CXX

#include "TestMCTruth.h"
#include "DataFormat/mctruth.h"
#include "TLorentzVector.h"

namespace larlite {

    bool TestMCTruth::initialize() {

        _fidvolBox = FidVolBox();

        evt_counter = 0;
        kept_events = 0;

        _h_mup = new TH2D("h_mup", "Reco vs True Mu Momentum (Reco From Proton Momentum & angles)", 100, 0, 3, 100, 0, 3);
        _h_mup_diff = new TH1F("h_mup_diff", "Reco - True / True Mu Momentum [GeV]", 100, -1, 1);

        _tree = new TTree("tree", "tree");
        _tree->Branch("true_mu_TotE", &true_mu_TotE, "true_mu_TotE/D");
        _tree->Branch("true_p_TotE", &true_p_TotE, "true_p_TotE/D");
        _tree->Branch("true_n_TotE", &true_n_TotE, "true_n_TotE/D");
        _tree->Branch("true_p_p", &true_p_p, "true_p_p/D");
        _tree->Branch("thetap", &thetap, "thetap/D");
        _tree->Branch("true_mu_p", &true_mu_p, "true_mu_p/D");
        _tree->Branch("thetamu", &thetamu, "thetamu/D");

        _tree->Branch("absphidiff_minuspi", &absphidiff_minuspi, "absphidiff_minuspi/D");
        _tree->Branch("reco_mu_p", &reco_mu_p, "reco_mu_p/D");


        _tree->Branch("premnant", &premnant, "premnant/D");
        _tree->Branch("ptremnant", &ptremnant, "ptremnant/D");
        _tree->Branch("plane_Z_anglediff", &plane_Z_anglediff, "plane_Z_anglediff/D");
        _tree->Branch("true_nu_E", &true_nu_E, "true_nu_E/D");
        _tree->Branch("Eremnant", &Eremnant, "Eremnant/D");
        _tree->Branch("Eccqe_truemumom", &Eccqe_truemumom, "Eccqe_truemumom/D");
        _tree->Branch("Eccqe_recomumom", &Eccqe_recomumom, "Eccqe_recomumom/D");

        _tree->Branch("true_p_KE", &true_p_KE, "true_p_KE/D");
        _tree->Branch("reco_mu_TotE", &reco_mu_TotE, "reco_mu_TotE/D");
        _tree->Branch("reco_nu_E_4momentum", &reco_nu_E_4momentum, "reco_nu_E_4momentum/D");
        _tree->Branch("reco_mu_TotE_quadratic", &reco_mu_TotE_quadratic, "reco_mu_TotE_quadratic/D");
        _tree->Branch("reco_mu_p_quadratic", &reco_mu_p_quadratic, "reco_mu_p_quadratic/D");


        _tree->Branch("thetamup", &thetamup, "thetamup/D");
        _tree->Branch("bindinoE", &bindinoE, "bindinoE/D");
        _tree->Branch("bindinop", &bindinop, "bindinop/D");

        _tree->Branch("genieCCQE", &genieCCQE, "genieCCQE/O");
        _tree->Branch("genieMomDiff", &genieMomDiff, "genieMomDiff/D");
        _Ecalc = new NuEnergyCalc();
        return true;
    }

    bool TestMCTruth::analyze(storage_manager* storage) {
        absphidiff_minuspi = -999.;
        reco_mu_p = -999.;
        true_mu_p = -999.;
        true_p_p = -999.;
        premnant = -999.;
        ptremnant = -999.;
        Eremnant = -999.;
        plane_Z_anglediff = -999.;
        true_nu_E = -999.;
        Eccqe_truemumom = -999.;
        Eccqe_recomumom = -999.;
        true_mu_TotE = -999.;
        true_p_TotE = -999.;
        true_p_KE = -999.;
        reco_mu_TotE = -999.;
        reco_nu_E_4momentum = -999.;
        reco_mu_TotE_quadratic = -999.;
        reco_mu_p_quadratic = -999.;
        bindinoE = -999.;
        bindinop = -999.;
        genieCCQE = false;
        genieMomDiff = -999.;


        double mumass_MEV = 105.6583715;
        double pmass_MEV = 938.272;
        double nmass_MEV = 939.565;


        std::cout << " --- New event --- " << std::endl;
        evt_counter++;
        auto ev_mctruth = storage->get_data<event_mctruth>("generator");
        if (!ev_mctruth) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
            return false;
        }
        if (ev_mctruth->size() != 1) {
            // print(larlite::msg::kERROR, __FUNCTION__, Form("MCTruth size is not equal to 1... it equals %lu!", ev_mctruth->size()));
            return false;
        }


        /// First select a good event (this is basically a filter)

        // Require the neutrino interaction is inside the fiducial volume
        // if (!_fidvolBox.Contain(ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Position()))
        //     return false;

        // Require two particles in the final state
        TLorentzVector p_momentum;
        TLorentzVector mu_momentum;
        TLorentzVector remnant_momentum;
        TLorentzVector bindino_momentum;
        TLorentzVector nu_momentum;
        TLorentzVector neutron_momentum;

        size_t n_muons = 0;
        size_t n_protons = 0;
        for (auto const& particle : ev_mctruth->at(0).GetParticles()) {

            if ( particle.StatusCode() == 15 ) {
                remnant_momentum = particle.Trajectory().at(0).Momentum();
                premnant = particle.Trajectory().at(0).Momentum().Vect().Mag();
                ptremnant = premnant * sin(particle.Trajectory().at(0).Momentum().Vect().Theta());
                Eremnant = particle.Trajectory().at(0).E();
                true_n_TotE = pow(pow(premnant, 2) + pow(nmass_MEV / 1000., 2), 0.5);
            }

            if (particle.StatusCode() == 11 && particle.PdgCode() == 2112)
                neutron_momentum = particle.Trajectory().at(0).Momentum();

            if (abs(particle.PdgCode()) == 2000000101) {
                // std::cout<<"Found a bindino! particle energy is "
                // <<particle.Trajectory().at(0).E()
                // <<" mass is "<<particle.Mass()
                // <<" p is " << particle.Trajectory().at(0).Momentum().Vect().Mag()<<std::endl;
                bindinoE = particle.Trajectory().at(0).E();
                bindinop = particle.Trajectory().at(0).Momentum().Vect().Mag();
                bindino_momentum = particle.Trajectory().at(0).Momentum();
                // std::cout<<"bindino status code is "<<particle.StatusCode()<<std::endl;
            }
            // Only particles with status code 1 are relevant
            if ( particle.StatusCode() != 1) continue;

            // //Note: this KE is in units of GEV!
            // double KE = particle.Trajectory().at(0).E() - particle.Mass();

            // //Don't care about any particles with less than 20 MeV KE
            // if ( KE < 0.02 ) continue;

            //Count up the number of muons
            if ( particle.PdgCode() == 13 ) {
                mu_momentum = particle.Trajectory().front().Momentum();
                n_muons++;
            }

            //Count up the number of protons
            if ( particle.PdgCode() == 2212 ) {
                //only consider protons above 100 MeV
                // if ( KE < 0.100 ) continue;
                p_momentum = particle.Trajectory().front().Momentum();
                n_protons++;
            }

            //Skip this event if any other particles exist... this filter
            //is to pick out 1mu1pNn events!
            // 22  => gammas
            // 211 => charged pions
            // 111 => pi0s
            // 11  => electrons
            // 321 => kaons
            if ( abs(particle.PdgCode()) == 22  ||
                    abs(particle.PdgCode()) == 211 ||
                    abs(particle.PdgCode()) == 111 ||
                    abs(particle.PdgCode()) == 11  ||
                    abs(particle.PdgCode()) == 321 ||
                    abs(particle.PdgCode()) == 2112 ||
                    //abs(particle.PdgCode()) == 2000000101 || //kbindino
                    abs(particle.PdgCode()) == 14 ||
                    abs(particle.PdgCode()) == 12 ||
                    particle.PdgCode() == -13 ||
                    //particle.PdgCode() == 1000180400 || //argon 40
                    abs(particle.PdgCode()) == 311 ||
                    abs(particle.PdgCode()) == 3222 || //strange baryons
                    abs(particle.PdgCode()) == 3122 ) //lambda baryon
                return false;
            // else {
            //     if (particle.PdgCode() != 13 && particle.PdgCode() != 2212)
            //         std::cout << "Found particle with PDG code " << particle.PdgCode() << std::endl;
            // }

        }

        //Skip this event if the number of muons is not 1 or # protons is not 1
        // std::cout<<"n_muons is "<<n_muons<<" and n_protons is "<<n_protons<<std::endl;
        if ( n_muons != 1 || n_protons != 1 )
            return false;

        kept_events++;
        /// OK now we have a good event
        // std::cout << "Found a good event!" << std::endl;
        // std::cout << "evt_counter is " << evt_counter << std::endl;
        true_nu_E = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().E();
        nu_momentum = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Momentum();
        genieCCQE = !ev_mctruth->at(0).GetNeutrino().CCNC() && !ev_mctruth->at(0).GetNeutrino().Mode();//ev_mctruth->at(0).GetNeutrino().InteractionType() == 1001;
        // std::cout<<        ev_mctruth->at(0).GetNeutrino().InteractionType()<<std::endl;

        // TVector3 zdir(0, 0, 1.);
        TVector3 muvec = mu_momentum.Vect();
        TVector3 pvec = p_momentum.Vect();
        thetamu = muvec.Theta();
        // std::cout << "thetamu with anglezdir is " << thetamu << std::endl;
        // std::cout << "meanwhile thetamu with muvec.Theta() is " << muvec.Theta() << std::endl;
        thetap  = pvec.Theta();
        thetamup = thetamu + thetap; //muvec.Angle(pvec);
        // std::cout << "in main, thetamu, thetap, thetamup are "
        //           << thetamu << ","
        //           << thetap << ","
        //           << thetamup << std::endl;
        true_mu_TotE = mu_momentum.E();
        true_p_TotE = p_momentum.E();
        true_p_KE = p_momentum.E() - pmass_MEV / 1000.;

        // std::cout<<"thetamu is "<<thetamu<<", thetap is "<<thetap<<std::endl;
        // std::cout<<"proton TOTAL ENERGY is pvec mag squared plus proton mass squared sqrt is "
        // <<std::pow(std::pow(pvec.Mag(),2) + std::pow(pmass_MEV/1000.,2),0.5)<<std::endl;
        // std::cout<<"p energy from TLorentzVector is "<<p_momentum.E()<<std::endl;
        reco_mu_p = pvec.Mag() * sin(thetap) / sin(thetamu);
        true_mu_p = muvec.Mag();
        true_p_p = pvec.Mag();

        reco_nu_E_4momentum = _Ecalc->ComputeEnu1mu1p(muvec, true_mu_p * 1000., pvec, pvec.Mag() * 1000.);
        reco_mu_TotE = pow(pow(reco_mu_p, 2) + pow(mumass_MEV / 1000., 2), 0.5);
        Eccqe_recomumom = _Ecalc->ComputeECCQE(reco_mu_TotE * 1000., muvec, false);
        Eccqe_truemumom = _Ecalc->ComputeECCQE(mu_momentum.E() * 1000., muvec, false);
        absphidiff_minuspi = std::abs(::geoalgo::Vector(pvec).Phi() - ::geoalgo::Vector(muvec).Phi()) - 3.14159265;
        // std::cout<<"Reco mu E is "<<reco_mu_TotE<<", true mu E is "<<mu_momentum.E()<<std::endl;
        // std::cout << "Reco mu p is " << reco_mu_p << ", true mu p is " << true_mu_p << std::endl;
        _h_mup->Fill(reco_mu_p, muvec.Mag());
        _h_mup_diff->Fill((reco_mu_p - muvec.Mag()) / muvec.Mag());

        TVector3 cross = muvec.Cross(pvec);

        plane_Z_anglediff = cross.Theta();

        // std::cout << "Computing muon total energy ... " << std::endl;
        // std::cout << "Note that true mu total energy is " << true_mu_TotE * 1000. << " MEV." << std::endl;
        // reco_mu_TotE_quadratic = _Ecalc->ComputeEmu1mu1pQuadraticNumeric(muvec, pvec, pvec.Mag() * 1000.);
        // reco_mu_TotE_quadratic = _Ecalc->ComputeEmu1mu1pQuadratic(muvec, pvec, pvec.Mag() * 1000.);
        // reco_mu_TotE_quadratic = _Ecalc->ComputeEmu1mu1pQuadraticIterative(muvec, pvec, pvec.Mag() * 1000.) / 1000.;
        // std::cout<<"Using muon direction, proton direction, proton momentum mag... iterative yields "<<reco_mu_TotE_quadratic<<std::endl;
        reco_mu_p_quadratic = std::sqrt(std::pow(reco_mu_TotE_quadratic, 2) - std::pow(mumass_MEV / 1000., 2));
        //consider true_n_TotE*1000. below for first argument to check for "exact matching" in the future
        reco_mu_TotE_quadratic = _Ecalc->ComputeEmu1mu1pQuadraticIterative(938. + 25., true_p_TotE * 1000., thetamu, thetap, true_p_p * 1000.) / 1000.;
        // std::cout << "Beginning with the following: (En = " << 938 + 25
        //           << ", Ep = " << true_p_TotE * 1000.
        //           << ", thetamu = " << thetamu
        //           << ", thetap = " << thetap
        //           << ", pp = " << true_p_p * 1000.
        //           << ")" << std::endl;
        // std::cout << "Iterative technique has found this result for Emu: " << reco_mu_TotE_quadratic * 1000. << std::endl;

        double avg_quad_mu_E = 0.;
        size_t ctr = 0;
        for (size_t i = 0; i <= 10; ++i) {
            double n_KE = 10. + i * 4.;
            double computed_E = _Ecalc->ComputeEmu1mu1pQuadraticIterative(938. + n_KE, true_p_TotE * 1000.,
                                thetamu, thetap,
                                true_p_p * 1000.);
            if (!std::isnan(computed_E)) {
                ctr++;
                avg_quad_mu_E += computed_E;
                // std::cout << " n_KE = " << n_KE << ", computed E = " << computed_E << std::endl;
            }
        }
        avg_quad_mu_E /= (double)ctr;
        // std::cout << "avg quad E is " << avg_quad_mu_E << std::endl;













        // std::cout << "--- Done computing muon total energy ---" << std::endl;
        // << ", reco_mu_TotE_quadratic is " << reco_mu_TotE_quadratic << std::endl;
        // TVector3 totalf = muvec + pvec;
        // TVector3 totali = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Momentum().Vect();
        // std::cout<<"Final momentum (mu + p) minus initial (nu) is :"<<std::endl;
        // (totalf-totali).Print();
        // std::cout << " -- Muon angle is " << thetamu << std::endl;
        // std::cout << "True muon TOTAL ENERGY is " << mu_momentum.E() << std::endl;
        // std::cout << "Muon momentum magnitude is " << muvec.Mag() << std::endl;
        // std::cout << "MUON transverse momentum component is " << muvec.Mag()*sin(thetamu) << std::endl;
        // std::cout << "Parallel component is " << muvec.Mag()*cos(thetamu) << std::endl;
        //std::cout << "(total energy) magnitude squared plus mass square sqrt is "
        //           << std::pow(std::pow(muvec.Mag(), 2) + std::pow(mumass_MEV / 1000., 2), 0.5) << std::endl;
        // std::cout << "PROTON transverse momentum component is " << pvec.Mag()*sin(thetap) << std::endl;
        printTruthStuff(*ev_mctruth);




        std::cout << "nu_momentum is :" << std::endl;
        nu_momentum.Print();
        std::cout << "bindino_momentum is :" << std::endl;
        bindino_momentum.Print();
        std::cout << "remnant_momentum is :" << std::endl;
        remnant_momentum.Print();
        std::cout << "from remnant, my articially constructed neutron 4 momentum is :" << std::endl;
        TLorentzVector neutron_4mom(
            -1.*remnant_momentum.Vect().X(),
            -1.*remnant_momentum.Vect().Y(),
            -1.*remnant_momentum.Vect().Z(),
            std::sqrt(remnant_momentum.Vect().Mag2() + 0.939565 * 0.939565)
        );
        std::cout << "mctruth neutron_momentum is :" << std::endl;
        neutron_momentum.Print();
        std::cout << "mu_momentum is :" << std::endl;
        mu_momentum.Print();
        std::cout << "p_momentum is :" << std::endl;
        p_momentum.Print();
        std::cout << "neutrino 4momentum plus mctruth neutron 4momentum minus (p + mu + bindino 4momentum) :" << std::endl;
        (nu_momentum + neutron_momentum - (p_momentum + mu_momentum + bindino_momentum) ).Print();
        std::cout << "neutrino 4momentum plus artificial neutron 4momentum minus (p + mu + bindino 4momentum) :" << std::endl;
        (nu_momentum + neutron_4mom - (p_momentum + mu_momentum + bindino_momentum) ).Print();
        genieMomDiff = (nu_momentum - (p_momentum + mu_momentum + remnant_momentum + bindino_momentum) ).Vect().Mag();

        // std::cout<<"TRUE MUON ENERGY (which below should match): "<<true_mu_TotE*1000.<<std::endl;
        // double poop = _Ecalc->ComputeEmu1mu1pQuadraticIterative(muvec, pvec, pvec.Mag() * 1000.);
        _tree->Fill();
        return true;
    }

    void TestMCTruth::printTruthStuff(const event_mctruth &mctruth) {
        // std::cout << "TRUE NEUTRINO ENERGY IS "
        //           << mctruth.at(0).GetNeutrino().Nu().Trajectory().front().E() << " GEV." << std::endl;
        std::cout << "----- LISTING ALL MCTRUTH PARTICLES IN EVENT:" << std::endl;
        TLorentzVector numom(0., 0., 0., 0.);
        TLorentzVector daughtmoms(0., 0., 0., 0.);
        for (auto const& particle : mctruth.at(0).GetParticles()) {

            // Only particles with status code 1 are relevant
            // if ( particle.StatusCode() != 1 && particle.StatusCode() != 0 && particle.StatusCode() != 15    ) continue;

            // //Note: this KE is in units of GEV!
            // double KE = particle.Trajectory().at(0).E() - particle.Mass();

            std::cout << " -- Particle with Status " << particle.StatusCode()
                      << " and PDG " << particle.PdgCode()
                      << " has TOTAL ENERGY " << particle.Trajectory().at(0).E() << " GEV."
                      << " and MASS " << particle.Mass() << std::endl;
            std::cout << " This particle's 4-Momentum is " << std::endl;
            particle.Trajectory().at(0).Momentum().Print();
            // std::cout << " This particle's 3-Momentum magnitude is "
            //           << particle.Trajectory().at(0).Momentum().Vect().Mag() << std::endl;

            if (particle.StatusCode() == 0 && particle.PdgCode() == 14)
                numom = particle.Trajectory().at(0).Momentum();
            if (particle.StatusCode() == 1 || particle.StatusCode() == 15)
                daughtmoms += particle.Trajectory().at(0).Momentum();
            // particle.Trajectory().at(0).Momentum().Vect().Print();
            // if(particle.StatusCode() == 15){
            //   std::cout<<"Particle with status code 15 has this momentum:"<<std::endl;
            //   particle.Trajectory().at(0).Momentum().Vect().Print();
            // }
        }
        std::cout << "Summary: neutrino 4momentum:" << std::endl;
        numom.Print();
        std::cout << "sum of status 1 and 15 4momenta:" << std::endl;
        daughtmoms.Print();

        std::cout << "----- DONE PRINTING TRUTH STUFF FOR THIS EVENT" << std::endl;
    }

    bool TestMCTruth::finalize() {

        std::cout << "End of TestMCTruth. Total events is "
                  << evt_counter << " while # of events with 1 mu and 1 p is "
                  << kept_events << std::endl;

        if (_fout) {
            _fout->cd();
            _h_mup->Write();
            _h_mup_diff->Write();
            _tree->Write();
        }
        return true;
    }

}
#endif
