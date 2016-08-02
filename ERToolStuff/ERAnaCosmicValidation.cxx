#ifndef ERTOOL_ERANACOSMICVALIDATION_CXX
#define ERTOOL_ERANACOSMICVALIDATION_CXX

#include "ERAnaCosmicValidation.h"

namespace ertool {

  ERAnaCosmicValidation::ERAnaCosmicValidation(const std::string& name)
    : AnaBase(name)
    , _result_tree(nullptr)
  {
    PrepareTreeVariables();

    // set default energy cut (for counting) to 0
    _eCut = 0;

    for (size_t i = 0; i < 4; ++i)
      _total_particles.push_back(0);

    _potential_secondary_nodeIDs.clear();

    evtno = 0;

    _tpc_box.Min( .3, -115., .5 );
    _tpc_box.Max( 256., 115., 1036. );

  }


  bool ERAnaCosmicValidation::Analyze(const EventData &data, const ParticleGraph &ps)
  {
    std::cout << " ----------- EVENT " << evtno << std::endl;
    evtno++;

    // First loop through the reconstructed cosmic primary particles just to get a list of all
    // MC secondaries that could potentially be associated with that particle
    // (this will be the denominator in our secondary tagger efficiency plot)
    _potential_secondary_nodeIDs.clear();
    ertool::Particle found_mcparticle;
    for ( auto const &pnodeid : ps.GetPrimaryNodes() ) {

      ertool::Particle primary = ps.GetParticle(pnodeid);

      if (primary.ProcessType() != ProcessType_t::kCosmic) continue;

      try { found_mcparticle = FindMCParticleGraphParticle(data, primary); }
      catch (ertool::ERException& e) {
        std::cout << "Skipping particle (in building secondary particle list)" << std::endl;
        continue;
      }

      /// This function grabs the mc secondaries and adds them to a std::set to be used later
      TrackAllSecondaries(found_mcparticle);
    }



    /// Grab the reconstructed particle graph
    auto const& particles = ps.GetParticleArray();
    /// Grab the mc particle graph
    auto const& mc_graph = MCParticleGraph();
    /// Grab the mc event data
    auto const& mc_data = MCEventData();

    if (!mc_graph.GetParticleArray().size()) {
      std::cout << "WARNING: Size of mc particle graph is zero! Perhaps you forgot to include mctruth/mctrack/mcshower?" << std::endl;
      return false;
    }

    //Loop through ALL particles in the reconstructed particle graph, fill ttree for each.
    for ( auto const & p : particles ) {

      // Reset ttree variables (since we will be filling the TTree once per particle)
      ResetTreeVariables();

      _is_track = p.RecoType() == RecoType_t::kTrack;
      //Particle is recod "primary" if it is tagged as cosmic, and tagged as a primary
      _reco_primary = (p.ProcessType() == ProcessType_t::kCosmic) && p.Primary();
      //Particle is recod "secondary" if it is a descdendant, and its parent is a cosmic
      _reco_secondary = p.Descendant() && ps.GetParticle(p.Parent()).ProcessType() == ProcessType_t::kCosmic;
      //Particle is recod "orphan" (specific process string added for ophan algo)
      _reco_orphan = p.ProcessType() == ProcessType_t::kCosmicOrphan;

      /// Find the same particle in the mcparticle event data and fill some more ttree variables
      /// that have to come from the mcparticle event data.
      ertool::Particle found_mcparticle;
      try {

        found_mcparticle = FindMCParticleGraphParticle(data, p);

        // Check the generation of the identified mcparticlegraph node to see if it's truly primary
        _true_primary = IsMCParticlePrimary(found_mcparticle);

        // To check if this mcparticle is a "secondary", I need to define "secondary"
        // Secondary means this mcparticle is a child of the mcparticle tagged as a primary cosmic.
        // That way, even if the "primary" tagged particle wasn't actually primary, my ERAlgoCRSecondary has still
        // done what it's supposed to if it finds a particle that is a secondary of that incorrectly tagged "primary"
        _true_secondary = IsMCParticleSecondary(found_mcparticle);

        // For now, if the particle isn't primary or secondary, it's considered orphan
        _true_orphan = !_true_primary && !_true_secondary;

        // Increment counters
        if (_true_primary) _total_particles.at(0)++;
        else if (_true_secondary) _total_particles.at(1)++;
        else _total_particles.at(2)++;
      }

      catch (ertool::ERException& e) {
        // std::cout << e.msg() << std::endl;
        std::cout << "Skipping this particle (in main particle loop)." << std::endl;
        _total_particles.at(3)++;
        continue;
      }




      /// Now we have the mcparticle (found_mcparticle), and the reco particle (p)
      /// And access to mcparticlegraph and reco particle graph
      /// So we fill the ttree with information

      /// Parent information:

      //    if (_true_primary && !_reco_primary) {

      // if (_is_track) {
      // std::cout << "This track should have been tagged as primary, but it wasn't." << std::endl;
      // data.Track(p).PrintInfo();
      if (found_mcparticle.Descendant()) {
        // std::cout << "the MC Parent of this track has the following info:" << std::endl;
        auto const mcparent = mc_graph.GetParticle(found_mcparticle.Parent());
        // if (mcparent.RecoType() == RecoType_t::kTrack) {
        //   std::cout << "  (it's a track) " << std::endl;
        //   mc_data.Track(mcparent).PrintInfo();
        // }
        // else if (mcparent.RecoType() == RecoType_t::kShower) {
        //   std::cout << "  (it's a shower) " << std::endl;
        //   mc_data.Shower(mcparent).PrintInfo();
        // }
        // else {
        //   std::cout << "  (it's neither track nor shower, so i print particle info )" << std::endl;
        //   std::cout<<mcparent.Print();
        // }
        _pdg_parent = mcparent.PdgCode();
      }
      else _pdg_parent = -1;
      _pdg_self = found_mcparticle.PdgCode();
      //}
      // else
      //   data.Shower(p).PrintInfo();


      // if (_true_secondary && !_reco_secondary) {
      // std::cout << " This particle should have been tagged secondary and it wasn't" << std::endl;
      // std::cout << " Its MC pdg is " << found_mcparticle.PdgCode() << std::endl;
      //   if (p.RecoType() == RecoType_t::kTrack) {
      //     std::cout << " It's a track. The info is printed:" << std::endl;
      //     data.Track(p).PrintInfo();
      //   }
      //   else if (p.RecoType() == RecoType_t::kShower) {
      //     std::cout << " It's a shower. The info is printed:" << std::endl;
      //     data.Shower(p).PrintInfo();
      //   }
      //   else std::cout << "wtffffff" << std::endl;
      //   std::cout << "   Parent info:" << std::endl;
      //   if (found_mcparticle.Descendant()) {
      //     auto const mcparent = mc_graph.GetParticle(found_mcparticle.Parent());
      //     std::cout << "   Parent pdg:"
      //               << mcparent.PdgCode()
      //               << " and node iD "
      //               << mcparent.ID()
      //               << std::endl;

      //     if (mcparent.RecoType() == RecoType_t::kTrack) {
      //       std::cout << " Parent is track. The info is printed:" << std::endl;
      //       data.Track(mcparent).PrintInfo();
      //     }
      //     else if (mcparent.RecoType() == RecoType_t::kShower) {
      //       std::cout << " Parent is shower.. The info is printed:" << std::endl;
      //       data.Shower(mcparent).PrintInfo();
      //     }

      //   }
      //   else std::cout << "   particle wasn't a descendant, so i didn't look for parent!" << std::endl;
      // }

      /// Particle energy information
      if (p.RecoType() == RecoType_t::kTrack) {
        _energy = mc_data.Track(found_mcparticle)._energy;
        _trk_length = mc_data.Track(found_mcparticle).Length();
        _start_x_vtx = mc_data.Track(found_mcparticle).at(0).at(0);
        _start_y_vtx = mc_data.Track(found_mcparticle).at(0).at(1);
        _start_z_vtx = mc_data.Track(found_mcparticle).at(0).at(2);
      }
      else if (p.RecoType() == RecoType_t::kShower) {
        _energy = mc_data.Shower(found_mcparticle)._energy;
        _trk_length = -999.;
        _start_x_vtx = mc_data.Shower(found_mcparticle).Start().at(0);
        _start_y_vtx = mc_data.Shower(found_mcparticle).Start().at(1);
        _start_z_vtx = mc_data.Shower(found_mcparticle).Start().at(2);
      }
      /// Distance to top wall information:

      // Finding perpendicular distance from particle start point to top wall
      ::geoalgo::Vector start;
      ::geoalgo::Vector end;
      ::geoalgo::Vector dir;
      if (p.RecoType() == RecoType_t::kTrack) {
        start = mc_data.Track(found_mcparticle).at(0);
        end = mc_data.Track(found_mcparticle).back();
        dir = mc_data.Track(found_mcparticle).at(0).Dir();
      }
      else if (p.RecoType() == RecoType_t::kShower) {
        start = mc_data.Shower(found_mcparticle).Start();
        end = end = mc_data.Shower(found_mcparticle).Start(); // end meaningless for showers
        dir = mc_data.Shower(found_mcparticle).Dir();
      }
      else std::cout << "p was neither shower nor track. type = " << p.RecoType() << std::endl;

      _perp_dist_any_wall = std::min(std::sqrt(_geoalg.SqDist(start, _tpc_box)),
                                     std::sqrt(_geoalg.SqDist(end, _tpc_box)));
      _going_upwards = dir.at(1) > 0;

      if (!_tpc_box.Contain(start)) {
        _start_contained_in_TPC = false;
        _perp_dist_top_wall = -999.;
        _orphan_dist = -999.;

      }
      else {
        _start_contained_in_TPC = true;
        _perp_dist_top_wall = _tpc_box.Max()[1] - start[1];
        if (dir[1] > 0) _orphan_dist = -100.;
        else _orphan_dist = (_tpc_box.Max()[1] - start[1]) / (dir[1] * -1.);
      }
      _trackend_contained_in_TPC = _tpc_box.Contain(end);


      // Loop through all other particles find shortest distance between each and this particle
      double _tmp_min_dist = 999999.;
      _dist_to_closest_particle = 999999.;
      for ( auto const &pnodeid2 : ps.GetParticleNodes() ) {

        ertool::Particle p2 = ps.GetParticle(pnodeid2);

        // Don't compare yourself to yourself
        if (p2.ID() == p.ID()) continue;

        if (p.RecoType() == RecoType_t::kTrack && p2.RecoType() == RecoType_t::kTrack)
          _tmp_min_dist =
            _findRel.FindClosestApproach(data.Track(p), data.Track(p2));
        else if (p.RecoType() == RecoType_t::kTrack && p2.RecoType() == RecoType_t::kShower)
          _tmp_min_dist =
            _findRel.FindClosestApproach(data.Track(p), data.Shower(p2));
        else if (p.RecoType() == RecoType_t::kShower && p2.RecoType() == RecoType_t::kTrack)
          _tmp_min_dist =
            _findRel.FindClosestApproach(data.Shower(p), data.Track(p2));
        else if (p.RecoType() == RecoType_t::kShower && p2.RecoType() == RecoType_t::kShower)
          _tmp_min_dist =
            _findRel.FindClosestApproach(data.Shower(p), data.Shower(p2));
        else std::cout << "wtf" << std::endl;

        if (_tmp_min_dist < _dist_to_closest_particle)
          _dist_to_closest_particle = _tmp_min_dist;

      }// End inner loop over particles to compare all permutations

      // Fill ttree
      _result_tree->Fill();
    }// End loop over reco particle graph

    return true;

  } // end analyze function

/// Simple function that takes in a particle from mcparticlegraph
/// and returns whether that particle is "primary"
// (if it doesn't have a parent that deposits energy in the detector)
  bool ERAnaCosmicValidation::IsMCParticlePrimary(const ertool::Particle & myparticle) {

    // If the particle is MC primary (currently MCParticles are never set as primary):
    if (myparticle.Primary()) return true;

    // If the particle doesn't descend from another, it's primary
    if (!myparticle.Descendant()) return true;

    // If the particle descends from another (IE it has a parent), do some more checks
    // std::cout << "mcparticle is a descendant... " << std::endl;
    /// Grab the mc event data
    auto const& mc_data = MCEventData();
    /// Grab the mc particle graph
    auto const& mc_graph = MCParticleGraph();
    /// Grab the particle's parent
    auto const& parent = mc_graph.GetParticle(myparticle.Parent());
    /// If the parent doesn't deposit energy in the detector (or is invisible), this particle is "primary"
    if (parent.RecoType() == RecoType_t::kTrack) {

      if (mc_data.Track(parent)._energy > 0) {
        // std::cout << "but its parent is a track that deposits energy " << mc_data.Track(parent)._energy
        //           << " so it is NOT primary" << std::endl;
        return false;
      }
      else {
        // std::cout << " But parent doesn't deposit energy so it IS primary." << std::endl;
        return true;
      }
    }
    else if (parent.RecoType() == RecoType_t::kShower) {

      if (mc_data.Shower(parent)._energy > 1) {
        // std::cout << "but its parent is a shower that deposits energy " << mc_data.Shower(parent)._energy
        //           << " so it is NOT primary" << std::endl;
        return false;
      }
      else {
        // std::cout << " But parent doesn't deposit energy so it IS primary." << std::endl;
        return true;
      }
    }
    else if (parent.RecoType() == RecoType_t::kInvisible) {
      return false;
    }
    else std::cout << " ---- wtf parent isn't track, shower, or invisible... ---- " << std::endl;

    // If we haven't figured out the particle is primary yet, then it isn't.
    // Only caste left is if parent isn't track, shower, or invisible
    return false;

  }

  /// Simple function that takes in a particle from mcparticlegraph
  /// and returns whether that particle is "secondary"
  /// (meaning, that particle has a parent in the mcparticlegraph that deposits energy in the detector)
  bool ERAnaCosmicValidation::IsMCParticleSecondary(const ertool::Particle & myparticle) {

    // I have a list of total possible secondaries (mcparticlegraph nodes that are children of the
    // mcparticlegraph nodes corresponding to reconstructed primary particles)
    // So to tell if an mcparticle is a secondary, I just ask if it's in the list.
    if ( _potential_secondary_nodeIDs.find(myparticle.ID()) != _potential_secondary_nodeIDs.end())
      return true;
    return false;

  }

  /// Simple function that takes in a particle from mcparticlegraph
  /// and returns whether that particle is "orphan"
  bool ERAnaCosmicValidation::IsMCParticleOrphan(const ertool::Particle & myparticle) {

    /// Right now, my definition of orphan is the same as definition of primary
    /// So to look at performance of orphan algo, you should do a ttree draw with a cut
    /// like _true_orphan && !_reco_primary
    /// (orphans are just primaries that the primary algo missed)
    return IsMCParticlePrimary(myparticle);

  }

  void ERAnaCosmicValidation::ProcessEnd(TFile * fout)
  {
    std::cout << "ERAnaCosmicValidation: Total particles analyzed counter:" << std::endl;
    std::cout << "    # Primary   = " << _total_particles.at(0) << std::endl;
    std::cout << "    # Secondary = " << _total_particles.at(1) << std::endl;
    std::cout << "    # Other     = " << _total_particles.at(2) << std::endl;
    std::cout << "    # Skipped   = " << _total_particles.at(3) << std::endl;
    std::cout << std::endl;

    if (fout) {
      fout->cd();
      _result_tree->Write();
    }
    return;
  }

  void ERAnaCosmicValidation::PrepareTreeVariables() {

    if (_result_tree) { delete _result_tree; }

    _result_tree = new TTree("cos_tree", "cos_tree");
    _result_tree->Branch("_is_track", &_is_track, "_is_track/O");
    _result_tree->Branch("_true_primary", &_true_primary, "_true_primary/O");
    _result_tree->Branch("_true_secondary", &_true_secondary, "_true_secondary/O");
    _result_tree->Branch("_true_orphan", &_true_orphan, "_true_orphan/O");
    _result_tree->Branch("_reco_primary", &_reco_primary, "_reco_primary/O");
    _result_tree->Branch("_reco_secondary", &_reco_secondary, "_reco_secondary/O");
    _result_tree->Branch("_reco_orphan", &_reco_orphan, "_reco_orphan/O");
    _result_tree->Branch("_pdg_parent", &_pdg_parent, "_pdg_parent/I");
    _result_tree->Branch("_pdg_self", &_pdg_self, "_pdg_self/I");
    _result_tree->Branch("_perp_dist_top_wall", &_perp_dist_top_wall, "_perp_dist_top_wall/F");
    _result_tree->Branch("_orphan_dist", &_orphan_dist, "_orphan_dist/F");
    _result_tree->Branch("_energy", &_energy, "_energy/F");
    _result_tree->Branch("_trk_length", &_trk_length, "_trk_length/F");
    _result_tree->Branch("_start_x_vtx", &_start_x_vtx, "_start_x_vtx/F");
    _result_tree->Branch("_start_y_vtx", &_start_y_vtx, "_start_y_vtx/F");
    _result_tree->Branch("_start_z_vtx", &_start_z_vtx, "_start_z_vtx/F");
    _result_tree->Branch("_start_contained_in_TPC", &_start_contained_in_TPC, "_start_contained_in_TPC/O");
    _result_tree->Branch("_trackend_contained_in_TPC", &_trackend_contained_in_TPC, "_trackend_contained_in_TPC/O");
    _result_tree->Branch("_perp_dist_any_wall", &_perp_dist_any_wall, "_perp_dist_any_wall/F");
    _result_tree->Branch("_going_upwards", &_going_upwards, "_going_upwards/O");
    _result_tree->Branch("_dist_to_closest_particle", &_dist_to_closest_particle, "_dist_to_closest_particle/F");

    return;
  }

  void ERAnaCosmicValidation::ResetTreeVariables() {

    _is_track = false;
    _true_primary = false;
    _true_secondary = false;
    _true_orphan = false;
    _reco_primary = false;
    _reco_secondary = false;
    _reco_orphan = false;
    _pdg_parent = -1;
    _pdg_self = -1;
    _perp_dist_top_wall = -999.;
    _orphan_dist = -999.;
    _energy = -999.;
    _trk_length = -999.;
    _start_x_vtx = -999.;
    _start_y_vtx = -999.;
    _start_z_vtx = -999.;
    _start_contained_in_TPC = false;
    _trackend_contained_in_TPC = false;
    _perp_dist_any_wall = -999.;
    _going_upwards = false;
    _dist_to_closest_particle = -999.;

    return;
  }


  ertool::Particle ERAnaCosmicValidation::FindMCParticleGraphParticle(const EventData & data, const ertool::Particle & p) {

    /// Grab the mc particle graph
    auto const& mc_graph = MCParticleGraph();
    /// Grab the mc event data
    auto const& mc_data = MCEventData();

    /// Loop through the MCParticleGraph and find the node that corresponds to this particle
    bool found_match = false;

    /// Container to store the particle from mcparticlegraph that matches my reco particle
    ertool::Particle found_mcparticle;
    ertool::Particle candidate_mcparticle;
    // Looping over the mcgraph subset with only this reco type is much faster
    // as most of the mcparticlegraph particles aren't associated with tracks or showers (if they don't enter the detector)
    for ( auto const & mcp : mc_graph.GetParticleNodes(p.RecoType()) ) {

      candidate_mcparticle = mc_graph.GetParticle(mcp);

      if ( candidate_mcparticle.RecoID() == p.RecoID() ){
        found_mcparticle = candidate_mcparticle;
        found_match = true;
        break;
      }

      // //If we're looking for a track
      // if ( p.RecoType() == RecoType_t::kTrack ) {
      //   //Check that the track we're looking for matches exactly the mcparticlegraph's track:
      //   //Check the start x, y, z coordinates all agree

      //   if (mc_data.Track(candidate_mcparticle).at(0).at(0) != data.Track(p).at(0).at(0) ||
      //       mc_data.Track(candidate_mcparticle).at(0).at(1) != data.Track(p).at(0).at(1) ||
      //       mc_data.Track(candidate_mcparticle).at(0).at(2) != data.Track(p).at(0).at(2)) continue;
      //   //Check energy agrees (note, mc_data.Track(candidate_mcparticle)._time is unfilled (1e308), but data.Track(p)._time IS filled)
      //   if (mc_data.Track(candidate_mcparticle)._energy != data.Track(p)._energy) continue;
      //   found_match = true;
      // }
      // else if ( p.RecoType() == RecoType_t::kShower ) { //If we're looking for a shower
      //   //Check that the shower we're looking for matches exactly the mcparticlegraph's shower:
      //   //Check the start point coordinates all agree
      //   if (mc_data.Shower(candidate_mcparticle).Start() != data.Shower(p).Start()) continue;
      //   //Check energy, dedx agree
      //   //(note, mc_data.Shower(candidate_mcparticle)._time is unfilled (1e308), but data.Shower(p)._time IS filled)
      //   if (mc_data.Shower(candidate_mcparticle)._energy != data.Shower(p)._energy ||
      //       mc_data.Shower(candidate_mcparticle)._dedx != data.Shower(p)._dedx ) continue;
      //   found_match = true;
      // }
      // else std::cout << " wtf you want to look for an mc particle matching a reco particle that is neither track nor shower?" << std::endl;

      // /// Note I have verified that only one match is ever found, so breaking here is fine to save computing time
      // if (found_match) {
      //   found_mcparticle = candidate_mcparticle;
      //   break;
      // }
    } // End loop over mcparticlegraph to find the mcparticlegraph node corresponding to this particle

    if (found_match) return found_mcparticle;
    else {
      throw ertool::ERException(Form("FindMCParticleGraphParticle failed to find a match! Particle info: %s.", p.Print().c_str()));
      return ertool::Particle();
    }

  }

/// Function to take in a particle from the mcparticle graph, and find all the mc secondaries associated with
/// that particle, and add them to a std::set keeping track of them all.
  void ERAnaCosmicValidation::TrackAllSecondaries(const ertool::Particle & p) {

    //If you want just one layer deep of children, use p.Children
    //If you want all the layers below this particle, you need to use ParticleGraph::GetAllDescendantNodes
    // (meaning you need access to the mcparticlegraph here)

     /// Grab the mc particle graph
    auto const& mc_graph = MCParticleGraph();
    // for (auto const& childnodeID : p.Children())
    for (auto const& childnodeID : mc_graph.GetAllDescendantNodes(p.ID()))
      _potential_secondary_nodeIDs.insert( childnodeID );

  }


}

#endif
