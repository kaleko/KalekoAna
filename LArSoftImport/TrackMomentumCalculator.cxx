#ifndef KALEKO_TRACKMOMENTUMCALCULATOR_CXX
#define KALEKO_TRACKMOMENTUMCALCULATOR_CXX

#include "TrackMomentumCalculator.h"
// \file  TrackMomentumCalculator.cxx
//
// \author sowjanyag@phys.ksu.edu

 Double_t xmeas[30]; Double_t ymeas[30]; Double_t eymeas[30]; Int_t nmeas;

double my_mcs_chi2( const double *x ) 
    {
      Double_t result = 0.0;
      
      Double_t p = x[0]; 
      
      Double_t theta0 = x[1]; 
      
      for ( Int_t i=0; i<nmeas; i++ )
        {
          Double_t xx = xmeas[i]; 
          
          Double_t yy = ymeas[i]; 
          
          Double_t ey = eymeas[i]; 
          
          Double_t rad_length = 14.0;
          
          Double_t l0 = xx/rad_length;
          
          Double_t res1 = 0.0;
          
          if ( xx>0 && p>0 ) res1 = ( 13.6/p )*sqrt( l0 )*( 1.0+0.038*TMath::Log( l0 ) );
          
          res1 = sqrt( res1*res1+theta0*theta0 );
          
          Double_t diff = yy-res1;
          
          if ( ey==0 ) { std::cout << " Zero denominator in my_mcs_chi2 ! " << std::endl; return -1; }
          
          Double_t incre = pow( diff/ey, 2.0 );
          
          result += incre;
          
        }
      
      result += 2.0/( 4.6 )*theta0; // *TMath::Log( 1.0/14.0 );
      
      if ( isnan( float(result) ) || isinf( float(result) ) ) { 
        std::cout<<"TrackMomentumCalculator: Is nan in my_mcs_chi2 ! "<<std::endl;
        return -1; }
        
      return result;
      
    }

    
namespace larlite {

namespace kaleko {

	TrackMomentumCalculator::TrackMomentumCalculator()
	{

		_myspline = new TrackMomentumSplines();
		n = 0;

		n_reco = 0;

		seg_stop = -1.0; n_seg = 0;

		gr_reco_xyz = new TPolyLine3D(); gr_reco_xy = new TGraph(); gr_reco_yz = new TGraph(); gr_reco_xz = new TGraph();

		gr_seg_xyz = new TPolyLine3D(); gr_seg_xy = new TGraph(); gr_seg_yz = new TGraph(); gr_seg_xz = new TGraph();

		steps_size = 5.0; n_steps = 9; for ( Int_t i = 1; i <= n_steps; i++ ) { steps.push_back( steps_size * i ); }

		basex.SetXYZ( 1.0, 0.0, 0.0 ); basey.SetXYZ( 0.0, 1.0, 0.0 ); basez.SetXYZ( 0.0, 0.0, 1.0 );

		nmeas = 0; p_mcs = -1.0; p_mcs_e = -1.0; chi2 = -1.0;

		steps_size2 = 10;//8.5; (Default)
		max_len_to_analyze = 9999.;

		p_mcs_2 = -1.0; LLbf = -1.0;

		kcal = 0.0022;

		// This pretty much has to be 100. If you use 20, you will find many events with low neutrino energy
		// that get reconstructed as having high neutrino energy becuase the muon only has 30cm in the
		// fiducial volume, and MCS way overestimates the energy. (you get high pion contamination above energy cut)
		minLength = 100.;
		maxLength = 1350.0;
		
		_debug_tree = 0;
		_debug_tree = new TTree("TMC_debug_tree","TMC_debug_tree");
		_debug_tree->Branch("full_track_len",&_full_track_len,"full_track_len/D");
		_debug_tree->Branch("full_MCS_E",&_full_MCS_E,"full_MCS_E/D");
		_debug_tree->Branch("full_range_E",&_full_range_E,"full_range_E/D");
		_debug_tree->Branch("delta_theta_x",&_delta_theta_x,"delta_theta_x/D");
		_debug_tree->Branch("delta_theta_y",&_delta_theta_y,"delta_theta_y/D");
		_debug_tree->Branch("seg_end_x",&_seg_end_x,"seg_end_x/D");
		_debug_tree->Branch("seg_end_y",&_seg_end_y,"seg_end_y/D");
		_debug_tree->Branch("seg_end_z",&_seg_end_z,"seg_end_z/D");
		_debug_tree->Branch("n_traj_points",&_n_traj_points,"n_traj_points/I");
		_debug_tree->Branch("seg_theta",&_seg_theta,"seg_theta/D");
		_debug_tree->Branch("seg_phi",&_seg_phi,"seg_phi/D");
		_debug_tree->Branch("counter",&_counter,"counter/I");
		_debug_tree->Branch("true_segment_E",&_true_segment_E,"true_segment_E/D");
		_debug_tree->Branch("true_predicted_RMS",&_true_predicted_RMS,"true_predicted_RMS/D");
		_debug_tree->Branch("segment_E",&_segment_E,"segment_E/D");
		_debug_tree->Branch("predicted_RMS",&_predicted_RMS,"predicted_RMS/D");
		_debug_tree->Branch("segment_E_fromMCS",&_segment_E_fromMCS,"segment_E_fromMCS/D");
		_debug_tree->Branch("predicted_RMS_fromMCS",&_predicted_RMS_fromMCS,"predicted_RMS_fromMCS/D");
		_debug_tree->Branch("resid_dist",&_resid_dist,"resid_dist/D");
		_debug_tree->Branch("llbf",&_llbf,"llbf/D");
		_debug_tree->Branch("run",&_run,"run/I");
		_debug_tree->Branch("subrun",&_subrun,"subrun/I");
		_debug_tree->Branch("eventid",&_eventid,"eventid/I");
		_counter = 0;


		//X is delta_theta/RMS
		double weight_x[89] = 
		{-2.9663,-2.8989,-2.8315,-2.7640,-2.6966,-2.6292,-2.5618,-2.4944,-2.4270,-2.3596,
-2.2921,-2.2247,-2.1573,-2.0899,-2.0225,-1.9551,-1.8876,-1.8202,-1.7528,-1.6854,
-1.6180,-1.5506,-1.4831,-1.4157,-1.3483,-1.2809,-1.2135,-1.1461,-1.0787,-1.0112,
-0.9438,-0.8764,-0.8090,-0.7416,-0.6742,-0.6067,-0.5393,-0.4719,-0.4045,-0.3371,
-0.2697,-0.2022,-0.1348,-0.0674,0.0000,0.0674,0.1348,0.2022,0.2697,0.3371,
0.4045,0.4719,0.5393,0.6067,0.6742,0.7416,0.8090,0.8764,0.9438,1.0112,
1.0787,1.1461,1.2135,1.2809,1.3483,1.4157,1.4831,1.5506,1.6180,1.6854,
1.7528,1.8202,1.8876,1.9551,2.0225,2.0899,2.1573,2.2247,2.2921,2.3596,
2.4270,2.4944,2.5618,2.6292,2.6966,2.7640,2.8315,2.8989,2.9663};
		//Y is weight
		double weight_y[89]	= {1.2948,0.9531,0.8876,1.0950,1.0707,1.3776,1.0803,0.9767,1.0002,1.1140,
1.1392,0.9604,0.9440,0.8483,0.9565,0.9652,0.8497,1.0676,1.0675,1.0576,
0.9237,1.0168,1.0425,1.0277,0.8664,0.9676,0.9572,0.9914,0.9495,0.8838,
0.9582,0.9178,0.9754,0.8876,1.0000,1.0689,0.9648,0.9589,0.9842,1.0007,
1.0031,1.0850,1.0267,0.9189,0.8550,0.9668,1.0200,1.1274,1.0625,1.0505,
1.0808,1.0473,1.0191,0.9933,1.0230,1.0017,0.9370,0.9068,0.9211,0.9291,
1.0140,0.9653,0.9494,0.9658,0.9390,1.0275,1.0173,1.0064,1.0097,1.0828,
1.0959,1.1180,1.0526,0.9352,1.0393,1.0099,1.0762,1.2356,1.1929,1.1611,
1.0467,1.1004,1.0995,1.1981,1.1229,1.1621,1.4986,1.1541,1.3291};

	  weight_graph = new TGraph(89,weight_x,weight_y);

	  gaus_smear = new TF1("gaus_smear","TMath::Gaus(x,0,10,true)",-200,200);

	  _run = -99999;
	  _subrun = -99999;
	  _eventid = -99999;

	}


	Double_t TrackMomentumCalculator::GetMuMultiScatterLLHD2( const larlite::mctrack &trk )
	{
		Double_t LLHD = -1.0;

		std::vector<Float_t> recoX; std::vector<Float_t> recoY; std::vector<Float_t> recoZ;

		recoX.clear(); recoY.clear(); recoZ.clear();

		Int_t n_points = trk.size();

		//    std::cout<<"kaleko: n_points = "<<n_points<<std::endl;

		for ( Int_t i = 0; i < n_points; i++ )
		{
			const TVector3 &pos = trk.at(i).Position().Vect();

			recoX.push_back( pos.X() ); recoY.push_back( pos.Y() ); recoZ.push_back( pos.Z() );

			// std::cout << " posX, Y, Z : " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << std::endl;

		}

		Int_t my_steps = recoX.size();
		//        std::cout<<"kaleko: my_steps = "<<my_steps<<std::endl;
		if ( my_steps < 2 ) return -1.0;

		Int_t check0 = GetRecoTracks( recoX, recoY, recoZ );
		//        std::cout<<"kaleko: check0 = "<<check0<<std::endl;
		if ( check0 != 0 ) return -1.0;

		seg_size = steps_size2;
		//    std::cout<<"kaleko: seg_size = "<<seg_size<<std::endl;

		//this fills segL
		Int_t check1 = GetSegTracks2( recoX, recoY, recoZ );
		//        std::cout<<"kaleko: check1 = "<<check1<<std::endl;
		if ( check1 != 0 ) return -1.0;

		Int_t seg_steps = segx.size();
		//        std::cout<<"kaleko: seg_steps = "<<seg_steps<<std::endl;
		if ( seg_steps < 2 ) return -1;

		Int_t seg_steps0 = seg_steps - 1;

		Double_t recoL = segL.at(seg_steps0);
		//        std::cout<<"kaleko: recoL = "<<recoL<<std::endl;
		if ( recoL < 20.0 || recoL > 1350.0 ) return -1;

		Int_t check2 = GetDeltaThetaij( dEi, dEj, dthij, seg_size, ind );
		//        std::cout<<"kaleko: check2 = "<<check2<<std::endl;
		if ( check2 != 0 ) return -1.0;

		Double_t p_range = recoL * kcal;

		Double_t logL = my_mcs_llhd( p_range, 0.5 );

		LLHD = logL;

		return LLHD;

	}

	Int_t TrackMomentumCalculator::GetRecoTracks( const std::vector<Float_t> &xxx, const std::vector<Float_t> &yyy, const std::vector<Float_t> &zzz )
	{
		Int_t a1 = xxx.size(); Int_t a2 = yyy.size(); Int_t a3 = zzz.size();
 
		if ( ( a1 != a2 ) || ( a1 != a3 ) || ( a2 != a3 ) ) { std::cout << " ( Get reco tacks ) Error ! " << std::endl; return -1; }

		n_reco = 0;

		for ( Int_t i = 0; i < a1; i++ )
		{
			Double_t nowx = xxx.at( i );

			Double_t nowy = yyy.at( i );

			Double_t nowz = zzz.at( i );

			x_reco[n_reco] = nowx;

			y_reco[n_reco] = nowy;

			z_reco[n_reco] = nowz;

			n_reco++;

		}

		gr_reco_xyz = new TPolyLine3D( n_reco, z_reco, x_reco, y_reco );

		gr_reco_yz = new TGraph( n_reco, z_reco, y_reco ); gr_reco_xz = new TGraph( n_reco, z_reco, x_reco ); gr_reco_xy = new TGraph( n_reco, x_reco, y_reco );

		return 0;

	}

	Int_t TrackMomentumCalculator::GetSegTracks2( const std::vector<Float_t> &xxx, const std::vector<Float_t> &yyy, const std::vector<Float_t> &zzz )
	{

		// std::cout<<"Start of GETSEGTRACKS2"<<std::endl;
		Double_t stag = 0.0;

		Int_t a1 = xxx.size(); Int_t a2 = yyy.size(); Int_t a3 = zzz.size();

		if ( ( a1 != a2 ) || ( a1 != a3 ) || ( a2 != a3 ) ) { std::cout << " ( Digitize reco tacks ) Error ! " << std::endl; return -1; }

		Int_t stopper = seg_stop / seg_size;

		Int_t a4 = a1 - 1;

		segx.clear(); segy.clear(); segz.clear();

		segnx.clear(); segny.clear(); segnz.clear();

		segL.clear();

		Int_t ntot = 0;

		n_seg = 0;

		Double_t x0; Double_t y0; Double_t z0;

		Double_t x00 = xxx.at( 0 ); Double_t y00 = yyy.at( 0 ); Double_t z00 = zzz.at( 0 );

		Int_t indC = 0;

		std::vector<Float_t> vx; std::vector<Float_t> vy; std::vector<Float_t> vz;

		vx.clear(); vy.clear(); vz.clear();

		for ( Int_t i = 0; i <= a4; i++ )
		{
			x0 = xxx.at( i ); y0 = yyy.at( i ); z0 = zzz.at( i );

			Double_t RR0 = sqrt( pow(x00 - x0, 2.0) + pow(y00 - y0, 2.0) + pow(z00 - z0, 2.0) );

			if ( RR0 >= stag )
			{
				segx.push_back( x0 ); segy.push_back( y0 ); segz.push_back( z0 );

				segL.push_back( stag );

				x_seg[ n_seg ] = x0; y_seg[ n_seg ] = y0; z_seg[ n_seg ] = z0;

				n_seg++;

				vx.push_back( x0 ); vy.push_back( y0 ); vz.push_back( z0 );

				ntot++;

				indC = i + 1;

				break;

			}

		}

		for ( Int_t i = indC; i < a4; i++ )
		{
			Double_t x1 = xxx.at( i ); Double_t y1 = yyy.at( i );	Double_t z1 = zzz.at( i );

			Double_t dr1 = sqrt( pow( x1 - x0, 2 ) + pow( y1 - y0, 2) + pow( z1 - z0, 2 ) );

			Double_t x2 = xxx.at( i + 1 ); Double_t y2 = yyy.at( i + 1 ); Double_t z2 = zzz.at( i + 1 );

			Double_t dr2 = sqrt( pow( x2 - x0, 2 ) + pow( y2 - y0, 2) + pow( z2 - z0, 2 ) );

			if ( dr1 < seg_size )
			{
				vx.push_back( x1 ); vy.push_back( y1 ); vz.push_back( z1 );

				ntot++;

			}

			if ( dr1 < seg_size && dr2 > seg_size )
			{
				// ..

				// std::cout << " 1 " << std::endl;

				Double_t t = 0.0;

				Double_t dx = x2 - x1; Double_t dy = y2 - y1; Double_t dz = z2 - z1;

				Double_t dr = sqrt( dx * dx + dy * dy + dz * dz );

				if ( dr == 0 ) { std::cout << " ( Zero ) Error ! " << std::endl; return -1; }

				Double_t beta = 2.0 * ( (x1 - x0) * dx + (y1 - y0) * dy + (z1 - z0) * dz ) / ( dr * dr );

				Double_t gamma = ( dr1 * dr1 - seg_size * seg_size ) / ( dr * dr );

				Double_t delta = beta * beta - 4.0 * gamma;

				if ( delta < 0.0 ) { std::cout << " ( Discriminant ) Error ! " << std::endl; return -1; }

				Double_t lysi1 = ( -beta + sqrt( delta ) ) / 2.0;

				t = lysi1;

				Double_t xp = x1 + t * dx;

				Double_t yp = y1 + t * dy;

				Double_t zp = z1 + t * dz;

				segx.push_back( xp ); segy.push_back( yp ); segz.push_back( zp );

				segL.push_back( 1.0 * n_seg * 1.0 * seg_size + stag );

				x_seg[ n_seg ] = xp; y_seg[ n_seg ] = yp; z_seg[ n_seg ] = zp; n_seg++;

				x0 = xp; y0 = yp; z0 = zp;

				vx.push_back( x0 ); vy.push_back( y0 ); vz.push_back( z0 );

				ntot++;

		// std::cout<<"This cout is inside of \"for ( Int_t i = indC; i < a4; i++ )\", ";
		// std::cout<<"then inside of \"if ( dr1 < seg_size && dr2 > seg_size )\""<<std::endl;
		// std::cout<<"  vx is (";
		// for (auto const& david : vx)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// std::cout<<"  vy is (";
		// for (auto const& david : vy)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// std::cout<<"  vz is (";
		// for (auto const& david : vz)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// std::cout<<"size of vx, vy, vz are "<<vx.size()<<", "<<vy.size()<<", "<<vz.size()<<std::endl;
				points_per_segment.push_back(vx.size());

				Double_t na = vx.size();

				Double_t sumx = 0.0;

				Double_t sumy = 0.0;

				Double_t sumz = 0.0;

				for ( Int_t i = 0; i < na; i++ )
				{
					Double_t xxw1 = vx.at( i );

					Double_t yyw1 = vy.at( i );

					Double_t zzw1 = vz.at( i );

					sumx += xxw1; sumy += yyw1; sumz += zzw1;

				}

				sumx = sumx / na; sumy = sumy / na; sumz = sumz / na;

				std::vector<Double_t> mx;

				std::vector<Double_t> my;

				std::vector<Double_t> mz;

				TMatrixDSym m( 3 );

				for ( Int_t i = 0; i < na; i++ )
				{
					Double_t xxw1 = vx.at( i ); Double_t yyw1 = vy.at( i ); Double_t zzw1 = vz.at( i );

					mx.push_back( xxw1 - sumx ); my.push_back( yyw1 - sumy ); mz.push_back( zzw1 - sumz );

					Double_t xxw0 = mx.at( i ); Double_t yyw0 = my.at( i ); Double_t zzw0 = mz.at( i );

					m( 0, 0 ) += xxw0 * xxw0 / na; m( 0, 1 ) += xxw0 * yyw0 / na; m( 0, 2 ) += xxw0 * zzw0 / na;

					m( 1, 0 ) += yyw0 * xxw0 / na; m( 1, 1 ) += yyw0 * yyw0 / na; m( 1, 2 ) += yyw0 * zzw0 / na;

					m( 2, 0 ) += zzw0 * xxw0 / na; m( 2, 1 ) += zzw0 * yyw0 / na; m( 2, 2 ) += zzw0 * zzw0 / na;

				}

				TMatrixDSymEigen me(m);

				TVectorD eigenval = me.GetEigenValues();

				TMatrixD eigenvec = me.GetEigenVectors();

				Double_t max1 = -666.0;

				Double_t ind1 = 0;

				for ( Int_t i = 0; i < 3; i++)
				{
					Double_t p1 = eigenval( i );

					if ( p1 > max1 ) { max1 = p1; ind1 = i; }

				}

				// std::cout << ind1 << std::endl;

				Double_t ax = eigenvec( 0, ind1 );

				Double_t ay = eigenvec( 1, ind1 );

				Double_t az = eigenvec( 2, ind1 );

				if ( segx.at(n_seg - 1) - segx.at(n_seg - 2) > 0 ) ax = TMath::Abs( ax );

				else ax = -1.0 * TMath::Abs( ax );

				if ( segy.at(n_seg - 1) - segy.at(n_seg - 2) > 0 ) ay = TMath::Abs( ay );

				else ay = -1.0 * TMath::Abs( ay );

				if ( segz.at(n_seg - 1) - segz.at(n_seg - 2) > 0 ) az = TMath::Abs( az );

				else az = -1.0 * TMath::Abs( az );

				segnx.push_back( ax ); segny.push_back( ay ); segnz.push_back( az );

				/// these are the direction cosines
				// std::cout<<"ax is "<<ax<<std::endl;
				// std::cout<<"ay is "<<ay<<std::endl;
				// std::cout<<"az is "<<az<<std::endl;

				// Double_t angx = find_angle( 1.0, ax ); Double_t angy = find_angle( 1.0, ay );

				// std::cout << angx*0.001*180.0/3.14 << std::endl;

				ntot = 0;

				vx.clear(); vy.clear(); vz.clear();

				vx.push_back( x0 ); vy.push_back( y0 ); vz.push_back( z0 );

				ntot++;

			}

			else if ( dr1 >= seg_size )
			{
				// ..

				// std::cout << " 2 " << std::endl;

				Double_t t = 0.0;

				Double_t dx = x1 - x0; Double_t dy = y1 - y0; Double_t dz = z1 - z0;

				Double_t dr = sqrt( dx * dx + dy * dy + dz * dz );

				if ( dr == 0 ) { std::cout << " ( Zero ) Error ! " << std::endl; return -1; }

				if ( dr != 0 ) t = seg_size / dr;

				Double_t xp = x0 + t * dx;

				Double_t yp = y0 + t * dy;

				Double_t zp = z0 + t * dz;

				segx.push_back( xp ); segy.push_back( yp ); segz.push_back( zp );

				segL.push_back( 1.0 * n_seg * 1.0 * seg_size + stag );

				x_seg[ n_seg ] = xp; y_seg[ n_seg ] = yp; z_seg[ n_seg ] = zp; n_seg++;

				x0 = xp; y0 = yp; z0 = zp;

				i = ( i - 1 );

				// ..

				vx.push_back( x0 ); vy.push_back( y0 ); vz.push_back( z0 );

				ntot++;

				points_per_segment.push_back(vx.size());

				Double_t na = vx.size();

				Double_t sumx = 0.0;

				Double_t sumy = 0.0;

				Double_t sumz = 0.0;

				for ( Int_t i = 0; i < na; i++ )
				{
					Double_t xxw1 = vx.at( i );

					Double_t yyw1 = vy.at( i );

					Double_t zzw1 = vz.at( i );

					sumx += xxw1; sumy += yyw1; sumz += zzw1;

				}

				sumx = sumx / na; sumy = sumy / na; sumz = sumz / na;

				std::vector<Double_t> mx;

				std::vector<Double_t> my;

				std::vector<Double_t> mz;

				TMatrixDSym m( 3 );

				for ( Int_t i = 0; i < na; i++ )
				{
					Double_t xxw1 = vx.at( i ); Double_t yyw1 = vy.at( i ); Double_t zzw1 = vz.at( i );

					mx.push_back( xxw1 - sumx ); my.push_back( yyw1 - sumy ); mz.push_back( zzw1 - sumz );

					Double_t xxw0 = mx.at( i ); Double_t yyw0 = my.at( i ); Double_t zzw0 = mz.at( i );

					m( 0, 0 ) += xxw0 * xxw0 / na; m( 0, 1 ) += xxw0 * yyw0 / na; m( 0, 2 ) += xxw0 * zzw0 / na;

					m( 1, 0 ) += yyw0 * xxw0 / na; m( 1, 1 ) += yyw0 * yyw0 / na; m( 1, 2 ) += yyw0 * zzw0 / na;

					m( 2, 0 ) += zzw0 * xxw0 / na; m( 2, 1 ) += zzw0 * yyw0 / na; m( 2, 2 ) += zzw0 * zzw0 / na;

				}

				TMatrixDSymEigen me(m);

				TVectorD eigenval = me.GetEigenValues();

				TMatrixD eigenvec = me.GetEigenVectors();

				Double_t max1 = -666.0;

				Double_t ind1 = 0;

				for ( Int_t i = 0; i < 3; i++)
				{
					Double_t p1 = eigenval( i );

					if ( p1 > max1 ) { max1 = p1; ind1 = i; }

				}

				// std::cout << ind1 << std::endl;

				Double_t ax = eigenvec( 0, ind1 );

				Double_t ay = eigenvec( 1, ind1 );

				Double_t az = eigenvec( 2, ind1 );

				if ( segx.at(n_seg - 1) - segx.at(n_seg - 2) > 0 ) ax = TMath::Abs( ax );

				else ax = -1.0 * TMath::Abs( ax );

				if ( segy.at(n_seg - 1) - segy.at(n_seg - 2) > 0 ) ay = TMath::Abs( ay );

				else ay = -1.0 * TMath::Abs( ay );

				if ( segz.at(n_seg - 1) - segz.at(n_seg - 2) > 0 ) az = TMath::Abs( az );

				else az = -1.0 * TMath::Abs( az );

				segnx.push_back( ax ); segny.push_back( ay ); segnz.push_back( az );

				// Double_t angx = find_angle( 1.0, ax ); Double_t angy = find_angle( 1.0, ay );

				// std::cout << angx*0.001*180.0/3.14 << std::endl;

				ntot = 0;

				vx.clear(); vy.clear(); vz.clear();

				vx.push_back( x0 ); vy.push_back( y0 ); vz.push_back( z0 );

				ntot++;

			}

			if ( n_seg >= ( stopper + 1.0 ) && seg_stop != -1 ) break;

		}

		gr_seg_xyz = new TPolyLine3D( n_seg, z_seg, x_seg, y_seg );

		gr_seg_yz = new TGraph( n_seg, z_seg, y_seg ); gr_seg_xz = new TGraph( n_seg, z_seg, x_seg ); gr_seg_xy = new TGraph( n_seg, x_seg, y_seg );


		// std::cout<<"This cout is at the very end of getsegtracks2:"<<std::endl;
		// std::cout<<"  segx is (";
		// for (auto const& david : segx)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// std::cout<<"  segy is (";
		// for (auto const& david : segy)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// std::cout<<"  segz is (";
		// for (auto const& david : segz)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;

		// std::cout<<"  vx is (";
		// for (auto const& david : vx)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// std::cout<<"  vy is (";
		// for (auto const& david : vy)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// std::cout<<"  vz is (";
		// for (auto const& david : vz)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;

		// std::cout<<" End of GetSegTracks2. "<<std::endl;




		return 0;

	}

	Int_t TrackMomentumCalculator::GetDeltaThetaij( std::vector<Float_t> &ei, std::vector<Float_t> &ej, std::vector<Float_t> &th, Double_t thick, std::vector<Float_t> &ind )
	{
		// std::cout<<" start of GetDeltaThetaij"<<std::endl;


		Int_t a1 = segx.size(); Int_t a2 = segy.size(); Int_t a3 = segz.size();

		if ( ( a1 != a2 ) || ( a1 != a3 ) ) { std::cout << " ( Get thij ) Error ! " << std::endl; return -1.0; }

		Int_t tot = a1 - 1; Double_t thick1 = thick + 0.13;

		ei.clear(); ej.clear(); th.clear(); ind.clear();

		for ( Int_t i = 0; i < tot; i++ )
		{
			Double_t dx = segnx.at( i ); Double_t dy = segny.at( i ); Double_t dz = segnz.at( i );

			TVector3 vec_z( dx, dy, dz );

			TVector3 vec_x;

			TVector3 vec_y;

			Double_t switcher = basex.Dot( vec_z );

			if ( switcher <= 0.995 )
			{
				vec_y = vec_z.Cross( basex ); vec_y = vec_y.Unit();

				vec_x = vec_y.Cross( vec_z );

			}

			else
			{
				// std::cout << " It switched ! Isn't this lovely !!! " << std::endl; // getchar();

				vec_y = basez.Cross( vec_z ); vec_y = vec_y.Unit();

				vec_x = vec_y.Cross( vec_z );

			}

			TVector3 Rx;

			Double_t ex = vec_x.Dot( basex ); Rx.SetX( ex );

			ex = vec_x.Dot( basey ); Rx.SetY( ex );

			ex = vec_x.Dot( basez ); Rx.SetZ( ex );

			TVector3 Ry;

			Double_t ey = vec_y.Dot( basex ); Ry.SetX( ey );

			ey = vec_y.Dot( basey ); Ry.SetY( ey );

			ey = vec_y.Dot( basez ); Ry.SetZ( ey );

			TVector3 Rz;

			Double_t ez = vec_z.Dot( basex ); Rz.SetX( ez );

			ez = vec_z.Dot( basey ); Rz.SetY( ez );

			ez = vec_z.Dot( basez ); Rz.SetZ( ez );

			Double_t refL = segL.at( i );

			for ( Int_t j = i; j < tot; j++ )
			{
				Double_t L1 = segL.at( j );

				Double_t L2 = segL.at( j + 1 );

				Double_t dz1 = L1 - refL;

				Double_t dz2 = L2 - refL;

				if ( dz1 <= thick1 && dz2 > thick1 )
				{
					Double_t here_dx = segnx.at( j );

					Double_t here_dy = segny.at( j );

					Double_t here_dz = segnz.at( j );

					TVector3 here_vec; here_vec.SetXYZ( here_dx, here_dy, here_dz );

					TVector3 rot_here; rot_here.SetXYZ( Rx.Dot( here_vec ), Ry.Dot( here_vec ), Rz.Dot( here_vec ) );

					Double_t scx = rot_here.X();

					Double_t scy = rot_here.Y();

					Double_t scz = rot_here.Z();

					Double_t azy = find_angle( scz, scy );

					Double_t azx = find_angle( scz, scx );

					Double_t ULim = 10000.0;

					Double_t LLim = -10000.0;

					Double_t cL = kcal;

					Double_t Li = segL.at( i );

					Double_t Lj = segL.at( j );

					if ( azy <= ULim && azy >= LLim )
					{
						ei.push_back( Li * cL );

						ej.push_back( Lj * cL );

						th.push_back( azy );

						ind.push_back( 2 );

					}

					if ( azx <= ULim && azx >= LLim )
					{
						ei.push_back( Li * cL );

						ej.push_back( Lj * cL );

						th.push_back( azx );

						ind.push_back( 1 );

					}

					break; // of course !

				}

			}

		}


		// std::cout<<"  ei is (";
		// for (auto const& david : ei)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;

		// std::cout<<"  ej is (";
		// for (auto const& david : ej)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;

		// std::cout<<"  th is (";
		// for (auto const& david : th)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;

		// std::cout<<"  ind is (";
		// for (auto const& david : ind)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;

		// std::cout<<"thick is "<<thick<<std::endl;

		// std::cout<<" END of GetDeltaThetaij"<<std::endl;

		return 0;

	}

	Double_t TrackMomentumCalculator::my_mcs_llhd( Double_t x0, Double_t x1, bool reweight )
	{
		_delta_theta_x = -99999.;
		_delta_theta_y = -99999.;


		Double_t p = x0;

		Double_t theta0x = x1;

		Double_t result = 0.0;

		Double_t nnn1 = dEi.size();

		Double_t red_length = steps_size2 / 14.0; //( 10.0 ) / 14.0;

		Double_t addth = 0;


		size_t non_outlier_deflection_counter = 0;

		for ( Int_t i = 0; i < nnn1; i++ )
		{
			double mydthij = dthij.at( i );
			// if you want to use smeared stuff!
			// double mydthij = smeared_dthij.at( i );

			//mydthij *= 1.06;
			// kaleko trying cutting on the tails of the distribution
			//if (TMath::Abs( mydthij ) > 150. ) continue; //continue;
			// kaleko cutting out the very center bin of the distribution
			//if (TMath::Abs( mydthij ) < 0.005 ) continue;
			// kaleko doing some bullshit
			//if (TMath::Abs( mydthij ) < 10. ) mydthij = 10.;

			Double_t Ei = p - dEi.at( i );

			Double_t Ej = p - dEj.at( i );

			// this should be a break? Ej is energy of this segment, Ei is energy of previous segment
			// this means the particle as ranged out. why make it go backwards?!
			if ( Ei > 0 && Ej < 0 ) {
				// std::cout<<"SETTING RESULT ENORMOUS AND BREAKING"<<std::endl;
				//i set result enormous because we're stopping looping over segments prematurely
				//each segment ADDS positive things to the result, so stopping early makes result
				//artificially small.
				result = 9999999999.;
				break;
				addth = 3.14 * 1000.0;
			}

			Ei = TMath::Abs( Ei );

			Ej = TMath::Abs( Ej );

			Double_t tH0 = ( 13.6 / sqrt( Ei * Ej ) ) * ( 1.0 + 0.038 * TMath::Log( red_length ) ) * sqrt( red_length );

			// Kaleko adding this ... if deltatheta/RMS is more than 3 (3 standard deviations away), skip this
			// NOTE!!! this doesn't work. this will just end up telling you all of your muons are very high energy
			// (doesn't fix data/MC discrepancy either)
			//if(TMath::Abs( mydthij ) / tH0 > 3.) continue;

			Double_t rms = -1.0;

			// if ( ind.at( i ) == 1 ) //kaleko: why only use x- scatter? ( 1 is x, 2 is y )
			// {
				rms = sqrt( tH0 * tH0 + pow( theta0x, 2.0 ) );

				Double_t DT = mydthij;// + addth; //addth functionaliy depricated

				// If this deflection is too large (theta/RMS more than 2.3 sigma), skip it
				// if(TMath::Abs(DT/rms) > 2.3) continue;
				// non_outlier_deflection_counter++;


				Double_t prob = my_g( DT, 0., rms );
				// std::cout<<"addth is "<<addth<<" and DT is "<<DT<<" and prob is "<<prob<<std::endl;

				// std::cout<<"DT is "
				// <<DT<<", RMS is "
				// <<rms<<", ratio is "
				// <<DT/rms<<", weight is "
				// <<weight_graph->Eval(DT/tH0)
				// <<std::endl;
				//if(reweight && TMath::Abs(DT/tH0) < 3.) prob *= weight_graph->Eval(DT/tH0);

				// note prob is NEGATIVE so you're adding stuff to result making it more positive
				// ultimately you take the smallest (positive) result
				result = result - 2.0 * prob;

			// }
			
		}

		if ( isnan( float( result ) ) || isinf( float( result ) ) ) { std::cout << " Is nan ! my_mcs_llhd ( 1 ) ! " << std::endl; getchar(); }

		// If only one segment went into this calculation, return a huge result (low likelihood)
		// if(non_outlier_deflection_counter < 2) result = 9999999999.;
		// std::cout<<"my_mcs_llhd returning a result of "<<result<<" / "<<non_outlier_deflection_counter
		// <<" = "<<result/non_outlier_deflection_counter<<std::endl;
		// Divide the result by the # of segments that went into it
		// Kaleko is adding this because when skipping deflections that are outliers, this skews the "result"
		// so it needs to be normalized by the number of segments that were considered in the likelihood
		// result /= non_outlier_deflection_counter;
		// std::cout<<"my_mcs_llhd returning a result of "<<result<<std::endl;
		return result;

	}

	Double_t TrackMomentumCalculator::my_g( Double_t xx, Double_t Q, Double_t s )
	{
		Double_t arg = 0.0;

		if ( s != 0 ) arg = ( xx - Q ) / s;

		else std::cout << " Error : The code tries to divide by zero ! " << std::endl;

		// std::cout<<"arg is "<<arg<<std::endl;		

		Double_t result = 0.0;

		if ( s != 0 ) result = -0.5 * TMath::Log( 2.0 * TMath::Pi() ) - TMath::Log( s ) - 0.5 * arg * arg;

		if ( isnan( float( result ) ) || isinf( float( result ) ) ) { std::cout << " Is nan ! my_g ! " << - TMath::Log( s ) << ", " << s << std::endl; getchar(); }

		// std::cout<<"my_g returning "<<result<<std::endl;
		// note my_g returns something NEGATIVE
		return result;

	}

void TrackMomentumCalculator::GetDeltaThetaRMS( Double_t &mean, Double_t &rms, Double_t &rmse, Double_t thick )
 {
   mean = 0.0; rms = 0.0; rmse = 0.0; 
   
   Int_t a1 = segx.size(); Int_t a2 = segy.size(); Int_t a3 = segz.size();
       
   if ( ( a1!=a2 ) || ( a1!=a3 ) ) { std::cout << " ( Get RMS ) Error ! " << std::endl; return; }
       
   Int_t tot = a1-1;
   
   Double_t thick1 = thick+0.13;
   
   std::vector<Float_t> buf0; buf0.clear();
   
   for ( Int_t i=0; i<tot; i++ )
     {
       Double_t dx = segnx.at( i ); Double_t dy = segny.at( i ); Double_t dz = segnz.at( i );
       
       TVector3 vec_z( dx, dy, dz ); 
                       
       TVector3 vec_x; 
       
       TVector3 vec_y; 
       
       Double_t switcher = basex.Dot( vec_z ); 
       
       if ( switcher<=0.995 ) 
         {
           vec_y = vec_z.Cross( basex ); vec_y = vec_y.Unit();  
           
           vec_x = vec_y.Cross( vec_z );
                   
         }
       
       else 
         {
           // std::cout << " It switched ! Isn't this lovely !!! " << std::endl; // getchar();
           
           vec_y = basez.Cross( vec_z ); vec_y = vec_y.Unit();  
           
           vec_x = vec_y.Cross( vec_z );
                           
         }
       
       TVector3 Rx;
       
       Double_t ex = vec_x.Dot( basex ); Rx.SetX( ex );
       
       ex = vec_x.Dot( basey ); Rx.SetY( ex );
       
       ex = vec_x.Dot( basez ); Rx.SetZ( ex );
       
       TVector3 Ry;
       
       Double_t ey = vec_y.Dot( basex ); Ry.SetX( ey );
 
       ey = vec_y.Dot( basey ); Ry.SetY( ey );
       
       ey = vec_y.Dot( basez ); Ry.SetZ( ey );
       
       TVector3 Rz;
       
       Double_t ez = vec_z.Dot( basex ); Rz.SetX( ez );
       
       ez = vec_z.Dot( basey ); Rz.SetY( ez );
       
       ez = vec_z.Dot( basez ); Rz.SetZ( ez );
               
       Double_t refL = segL.at( i );
       
       for ( Int_t j=i; j<tot; j++ )
         {
           Double_t L1 = segL.at( j );
           
           Double_t L2 = segL.at( j+1 );
                           
           Double_t dz1 = L1 - refL; 
           
           Double_t dz2 = L2 - refL; 
           
           if ( dz1<=thick1 && dz2>thick1 )
             {
               Double_t here_dx = segnx.at( j );
               
               Double_t here_dy = segny.at( j );
               
               Double_t here_dz = segnz.at( j );
               
               TVector3 here_vec; here_vec.SetXYZ( here_dx, here_dy, here_dz );
                                               
               TVector3 rot_here; rot_here.SetXYZ( Rx.Dot( here_vec ), Ry.Dot( here_vec ), Rz.Dot( here_vec ) );
                                               
               Double_t scx = rot_here.X(); 
           
               Double_t scy = rot_here.Y();  
               
               Double_t scz = rot_here.Z();  
               
               Double_t azy = find_angle( scz, scy ); azy*=1.0;
               
               Double_t azx = find_angle( scz, scx );
               
               Double_t ULim = 10000.0; Double_t LLim = -10000.0;
               
               // if ( azy<=ULim && azy>=LLim ) { buf0.push_back( azy ); } // hRMS.Fill( azy ); }
               
               if ( azx<=ULim && azx>=LLim ) { buf0.push_back( azx ); } // hRMS.Fill( azx ); }

               // if ( azy<=ULim && azy>=LLim && azx<=ULim && azx>=LLim && thick==5.0 ) { hSCA.Fill( azx, azy ); }
               
               // i=j-1;
               
               break; // of course !
               
             }
         }
    
     }
   
   Int_t nmeas = buf0.size(); Double_t nnn = 0.0;
   
   for ( Int_t i=0; i<nmeas; i++ ) { mean += buf0.at( i ); nnn++; } mean = mean/nnn;
   
   for ( Int_t i=0; i<nmeas; i++ ) rms += ( ( buf0.at( i ) )*( buf0.at( i ) ) );
   
   rms = rms/( nnn ); rms = sqrt( rms ); rmse = rms/sqrt( 2.0*tot );
   
   Double_t rms1 = rms;
   
   rms = 0.0;
   
   Double_t ntot1 = 0.0;
   
   Double_t lev1 = 2.50;

   for ( Int_t i=0; i<nmeas; i++ ) 
     {
       Double_t amp = buf0.at( i );
               
       if ( amp<( mean+lev1*rms1 ) && amp>( mean-lev1*rms1 )  ) 
         {
           ntot1++;
           
           rms += ( amp*amp );
           
         }
     
     }
   
   rms = rms/( ntot1 ); rms = sqrt( rms ); rmse = rms/sqrt( 2.0*ntot1 );
           
   return;
   
 }

	Double_t TrackMomentumCalculator::find_angle( Double_t vz, Double_t vy )
	{
		Double_t thetayz = -999.0;

		if ( vz > 0 && vy > 0 ) { Double_t ratio = TMath::Abs( vy / vz ); thetayz = TMath::ATan( ratio ); }

		else if ( vz < 0 && vy > 0 ) { Double_t ratio = TMath::Abs( vy / vz ); thetayz = TMath::ATan( ratio ); thetayz = 3.14159 - thetayz; }

		else if ( vz < 0 && vy < 0 ) { Double_t ratio = TMath::Abs( vy / vz ); thetayz = TMath::ATan( ratio ); thetayz = thetayz + 3.14159; }

		else if ( vz > 0 && vy < 0 ) { Double_t ratio = TMath::Abs( vy / vz ); thetayz = TMath::ATan( ratio ); thetayz = 2.0 * 3.14159 - thetayz; }

		else if ( vz == 0 && vy > 0 ) { thetayz = 3.14159 / 2.0; }

		else if ( vz == 0 && vy < 0 ) { thetayz = 3.0 * 3.14159 / 2.0; }

		if ( thetayz > 3.14159 ) { thetayz = thetayz - 2.0 * 3.14159; }

		Double_t result = 1000.0 * thetayz;

		return result;

	}

	Double_t TrackMomentumCalculator::GetMomentumMultiScatterLLHD( const larlite::mctrack &trk )
	{
		Double_t p = -1.0;

		std::vector<Float_t> recoX; std::vector<Float_t> recoY; std::vector<Float_t> recoZ;

		recoX.clear(); recoY.clear(); recoZ.clear();

		Int_t n_points = trk.size();

		for ( Int_t i = 0; i < n_points; i++ )
		{
			const TVector3 &pos = trk[i].Position().Vect();

			recoX.push_back( pos.X() ); recoY.push_back( pos.Y() ); recoZ.push_back( pos.Z() );

			// std::cout << " posX, Y, Z : " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << std::endl; getchar();

		}

		Int_t my_steps = recoX.size();

		if ( my_steps < 2 ) return -1.0;

		Int_t check0 = GetRecoTracks( recoX, recoY, recoZ );

		if ( check0 != 0 ) return -1.0;

		seg_size = steps_size2;

		Int_t check1 = GetSegTracks2( recoX, recoY, recoZ );

		if ( check1 != 0 ) return -1.0;

		Int_t seg_steps = segx.size();

		if ( seg_steps < 2 ) return -1;

		Int_t seg_steps0 = seg_steps - 1;

		Double_t recoL = segL.at(seg_steps0);

		if ( recoL < minLength || recoL > maxLength ) return -1;

		Int_t check2 = GetDeltaThetaij( dEi, dEj, dthij, seg_size, ind );

		if ( check2 != 0 ) return -1.0;

		Double_t logL = 1e+16;

		Double_t bf = -666.0; // Double_t errs = -666.0;

		Double_t start1 = 0.0; Double_t end1 = 750.0;

		Double_t start2 = 0.0; Int_t end2 = 0.0; // 800.0;

		for ( Int_t k = start1; k <= end1; k++ )
		{
			Double_t p_test = 0.001 + k * 0.01;

			for ( Int_t l = start2; l <= end2; l++ )
			{
				Double_t res_test = 2.0; // 0.001+l*1.0;

				Double_t fv = my_mcs_llhd( p_test, res_test );

				if ( fv < logL )
				{
					bf = p_test;

					logL = fv;

					// errs = res_test;

				}

			}

		}

		p_mcs_2 = bf; LLbf = logL;

		p = p_mcs_2;

		// Fill debug tree
		// if(debug){
		// 	std::cout<<"huzzah!"<<std::endl;
		// 	std::cout<<"ind size is "<<ind.size()<<std::endl;
		// 			for (auto const& david : ind)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// 	std::cout<<"segx size is "<<segx.size()<<std::endl;
		// 	for (auto const& david : segx)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// 	std::cout<<"points_per_segment size is "<<points_per_segment.size()<<std::endl;
		// 	for (auto const& david : points_per_segment)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;

		_llbf = LLbf;

			_full_track_len = (trk.back().Position().Vect() - trk.front().Position().Vect()).Mag();
			_full_range_E = _myspline->GetMuMomentum(_full_track_len) / 1000. + 0.106;
			_full_MCS_E = p;

			for(int david = 0; david < ind.size(); david+=2){
				_seg_end_x = segx.at(size_t(david/2.)+1);
				_seg_end_y = segy.at(size_t(david/2.)+1);
				_seg_end_z = segz.at(size_t(david/2.)+1);
				// std::cout<<"_seg_end_x is "<<_seg_end_x<<std::endl;
				// if(david < ind.size() ){

					_n_traj_points = points_per_segment.at(david/2);
	
					// std::cout<<"_n_traj_points is "<<_n_traj_points<<std::endl;
					TVector3 segend = TVector3(_seg_end_x,_seg_end_y,_seg_end_z);
					TVector3 segstart   = TVector3(segx.at(size_t(david/2.)),
												 segy.at(size_t(david/2.)),
												 segz.at(size_t(david/2.)));

					//quick find the MCTrack trajectory point closest to segstart
					double mindist = 999999.;
					size_t this_idx = 0;
					for(size_t ipt = 0; ipt < trk.size(); ipt++){
						auto const &pt = trk.at(ipt);
						double dist = (pt.Position().Vect() - segstart).Mag();
						if ( dist < mindist ) {
							mindist = dist;
							this_idx = ipt;
						}
					}
					_true_segment_E = (trk.at(this_idx).E() - 106.)/1000.;
					double redlen = steps_size2 / 14.;
					_true_predicted_RMS = ( 13.6 / (_true_segment_E + 0.106) ) * ( 1.0 + 0.038 * TMath::Log( redlen ) ) * sqrt( redlen );

					_segment_E = _full_range_E - 0.106 - (segstart - trk.front().Position().Vect()).Mag() * kcal;
					_predicted_RMS = ( 13.6 / (_segment_E+0.106) ) * ( 1.0 + 0.038 * TMath::Log( redlen ) ) * sqrt( redlen );

					_segment_E_fromMCS = _full_MCS_E - 0.106 - (segstart - trk.front().Position().Vect()).Mag() * kcal;
					_predicted_RMS_fromMCS = ( 13.6 / (_segment_E_fromMCS+0.106) ) * ( 1.0 + 0.038 * TMath::Log( redlen ) ) * sqrt( redlen );

					_resid_dist = (segstart - trk.back().Position().Vect()).Mag();

// std::cout<<"For this segment, the "<<this_idx<<"th point on the mctrack is closest, at "<<mindist<<"cm away."<<std::endl;
					// std::cout<<"segend x is "<<segend.X()<<std::endl;
					// std::cout<<"segstart x is "<<segstart.X()<<std::endl;
					_seg_theta = (segend-segstart).Theta();
					_seg_phi   = (segend-segstart).Phi();
				// }
				// else{
				// 	_n_traj_points = 999999;
				// 	_seg_theta = -99999.;
				// 	_seg_phi   = -99999.;
				// }
				

				if (ind.at(david) == 1 && david < dthij.size() && (david-1) > 0){
					// std::cout<<"david is "<<david<<std::endl;
					// std::cout<<"size of dthij is "<<dthij.size()<<std::endl;
					_delta_theta_x = dthij.at(david);
					_delta_theta_y = dthij.at(david-1);
					_debug_tree->Fill();
				}
				if (ind.at(david) == 2 && (david+1) < dthij.size() && david > 0){
					_delta_theta_x = dthij.at(david+1);
					_delta_theta_y = dthij.at(david);
					_debug_tree->Fill();
				}
			}
		// }
		return p;

	}

	Double_t TrackMomentumCalculator::GetMomentumMultiScatterLLHD( const larlite::track &trk, 
		bool flip, 
	bool debug, 
	bool reweight,
	int run,
	int subrun,
	int eventid)
	{
	// 	std::cout<<"Start of GetMomentumMultiScatterLLHD!"<<std::endl;
	// 	std::cout<<"Track Start: ("<<trk.Vertex().X()<<","<<trk.Vertex().Y()<<","<<trk.Vertex().Z()<<"), ";
	// 	std::cout<<"end: ("<<trk.End().X()<<","<<trk.End().Y()<<","<<trk.End().Z()<<")"<<std::endl;
		_counter++;
		points_per_segment.clear();

		Double_t p = -1.0;

		std::vector<Float_t> recoX; std::vector<Float_t> recoY; std::vector<Float_t> recoZ;

		recoX.clear(); recoY.clear(); recoZ.clear();

		Int_t n_points = trk.NumberTrajectoryPoints();
		// std::cout<<"n_points is "<<n_points<<std::endl;
		if (!flip) {
			const TVector3 &start = trk.Vertex();

			for ( Int_t i = 0; i < n_points; i++ )
			{
				const TVector3 &pos = trk.LocationAtPoint(i);
				if( (pos - start).Mag() > max_len_to_analyze ) break;

				recoX.push_back( pos.X() ); recoY.push_back( pos.Y() ); recoZ.push_back( pos.Z() );

				// cout << " posX, Y, Z : " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << endl; getchar();

			}
		}
		else {
			const TVector3 &start = trk.End();
			for ( Int_t i = n_points - 1; i >= 0; i-- )
			{
				const TVector3 &pos = trk.LocationAtPoint(i);
				if( (pos - start).Mag() > max_len_to_analyze ) break;
				recoX.push_back( pos.X() ); recoY.push_back( pos.Y() ); recoZ.push_back( pos.Z() );
			}

		}
		Int_t my_steps = recoX.size();
		// std::cout<<" my_Steps = "<<my_steps<<std::endl;
		if ( my_steps < 2 ) {
			// std::cout<<" my_Steps = "<<my_steps<<std::endl;
			return -1.0;
		}

		Int_t check0 = GetRecoTracks( recoX, recoY, recoZ );
		// std::cout<<" check0 = "<<check0<<std::endl;
		if ( check0 != 0 ) {
			// std::cout<<" check0 = "<<check0<<std::endl;
			return -1.0;
		}

		seg_size = steps_size2;

		Int_t check1 = GetSegTracks2( recoX, recoY, recoZ );
		// std::cout<<"check1 = "<<check1<<std::endl;
		if ( check1 != 0 ) {
			// std::cout<<"check1 = "<<check1<<std::endl;
			return -1.0;
		}

		Int_t seg_steps = segx.size();
		// std::cout<<"seg_steps = "<<seg_steps<<std::endl;
		if ( seg_steps < 2 ) {
			// std::cout<<"seg_steps = "<<seg_steps<<std::endl;
			return -1;
		}

		Int_t seg_steps0 = seg_steps - 1;

		Double_t recoL = segL.at(seg_steps0);
		// std::cout<<"debug: recoL is "<<recoL<<std::endl;
		// std::cout<<"  segL is (";
		// for (auto const& david : segL)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;

		if ( recoL < minLength || recoL > maxLength ) {
			// std::cout<<"length problem. recoL is "<<recoL<<std::endl;
			// std::cout<<"note track len is "<<(trk.Vertex()-trk.End()).Mag()<<std::endl;
			return -1;
		}

		Int_t check2 = GetDeltaThetaij( dEi, dEj, dthij, seg_size, ind );

		// fill a smeared dthij vector so you don't have to re-draw random for every raster scan step
		smeared_dthij.clear();
		for(size_t david = 0; david < dthij.size(); david++)
			smeared_dthij.push_back( dthij.at(david) + gaus_smear->GetRandom() );


		// std::cout<<"end of GetDeltaThetaij. dthij is :"<<std::endl;
		// std::cout<<" (";
		// for(auto const david : dthij)
		// 	std::cout<<david<<",";
		// std::cout<<")"<<std::endl;
		// std::cout<<"and ind is :"<<std::endl;
		// std::cout<<" (";
		// for(auto const david : ind)
		// 	std::cout<<david<<",";
		// std::cout<<")"<<std::endl;

		// std::cout<<"check2 = "<<check2<<std::endl;
		if ( check2 != 0 ) {
			std::cout<<"check2 = "<<check2<<std::endl;
			return -1.0;
		}

		Double_t logL = 1e+16;

		Double_t bf = -666.0; // Double_t errs = -666.0;

		Double_t start1 = 0.0; Double_t end1 = 750.0;

		Double_t start2 = 0.0; Int_t end2 = 0.0; // 800.0;

		for ( Int_t k = start1; k <= end1; k++ )// += 100)//k++)
		{
			Double_t p_test = 0.001 + k * 0.01;
			// std::cout<<"this p_test is "<<p_test<<std::endl;

			for ( Int_t l = start2; l <= end2; l++ )
			{
				Double_t res_test = 2; // 0.001+l*1.0; // is this the resolution parameter?

				Double_t fv = my_mcs_llhd( p_test, res_test, reweight );

				if ( fv < logL )
				{
					bf = p_test;

					logL = fv;

					// errs = res_test;

				}

			}

		}
		// std::cout<<"best p was found to be "<<bf<<std::endl;
		p_mcs_2 = bf; LLbf = logL;

		p = p_mcs_2;

		// std::cout<<"MCS decided this track has momentum "<<p<<std::endl;

		// Fill debug tree
		if(debug){

			_run = run;
			_subrun = subrun;
			_eventid = eventid;

			_llbf = LLbf;

		// 	std::cout<<"huzzah!"<<std::endl;
		// 	std::cout<<"ind size is "<<ind.size()<<std::endl;
		// 			for (auto const& david : ind)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// 	std::cout<<"segx size is "<<segx.size()<<std::endl;
		// 	for (auto const& david : segx)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;
		// 	std::cout<<"points_per_segment size is "<<points_per_segment.size()<<std::endl;
		// 	for (auto const& david : points_per_segment)
		// 	std::cout<<david<<", ";
		// std::cout<<")"<<std::endl;

			_full_track_len = (trk.End() - trk.Vertex()).Mag();
			_full_MCS_E = p;
			_full_range_E = _myspline->GetMuMomentum(_full_track_len) / 1000. + 0.106;

			for(int david = 0; david < ind.size(); david+=2){
				_seg_end_x = segx.at(size_t(david/2.)+1);
				_seg_end_y = segy.at(size_t(david/2.)+1);
				_seg_end_z = segz.at(size_t(david/2.)+1);
				// std::cout<<"_seg_end_x is "<<_seg_end_x<<std::endl;
				// if(david < ind.size() ){
					_n_traj_points = points_per_segment.at(david/2);
					// std::cout<<"_n_traj_points is "<<_n_traj_points<<std::endl;
					TVector3 segend = TVector3(_seg_end_x,_seg_end_y,_seg_end_z);
					TVector3 segstart   = TVector3(segx.at(size_t(david/2.)),
												 segy.at(size_t(david/2.)),
												 segz.at(size_t(david/2.)));
					// std::cout<<"segend x is "<<segend.X()<<std::endl;
					// std::cout<<"segstart x is "<<segstart.X()<<std::endl;
					_seg_theta = (segend-segstart).Theta();
					_seg_phi   = (segend-segstart).Phi();


					double redlen = steps_size2 / 14.;
					_segment_E = _full_range_E - 0.106 - (segstart - trk.Vertex()).Mag() * kcal;
					_predicted_RMS = ( 13.6 / (_segment_E+0.106) ) * ( 1.0 + 0.038 * TMath::Log( redlen ) ) * sqrt( redlen );

					_segment_E_fromMCS = _full_MCS_E - 0.106 - (segstart - trk.Vertex()).Mag() * kcal;
					_predicted_RMS_fromMCS = ( 13.6 / (_segment_E_fromMCS+0.106) ) * ( 1.0 + 0.038 * TMath::Log( redlen ) ) * sqrt( redlen );

					_resid_dist = (segstart - trk.End()).Mag();

				// }
				// else{
				// 	_n_traj_points = 999999;
				// 	_seg_theta = -99999.;
				// 	_seg_phi   = -99999.;
				// }
				

				if (ind.at(david) == 1 && david < dthij.size() && (david-1) > 0){
					// std::cout<<"david is "<<david<<std::endl;
					// std::cout<<"size of dthij is "<<dthij.size()<<std::endl;
					_delta_theta_x = dthij.at(david);
					_delta_theta_y = dthij.at(david-1);
					_debug_tree->Fill();
				}
				if (ind.at(david) == 2 && (david+1) < dthij.size() && david > 0){
					_delta_theta_x = dthij.at(david+1);
					_delta_theta_y = dthij.at(david);
					_debug_tree->Fill();
				}
			}
		}

		return p;

	}

	Double_t TrackMomentumCalculator::GetMomentumMultiScatterChi2( const larlite::mctrack &trk ){
		  Double_t p = -1.0;

   std::vector<Float_t> recoX; std::vector<Float_t> recoY; std::vector<Float_t> recoZ;

   recoX.clear(); recoY.clear(); recoZ.clear();

   Int_t n_points = trk.size();

   for ( Int_t i=0; i<n_points; i++ )
     {
       const TVector3 &pos = trk[i].Position().Vect();

       recoX.push_back( pos.X() ); recoY.push_back( pos.Y() ); recoZ.push_back( pos.Z() );

       // cout << " posX, Y, Z : " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << endl;

     }

   Int_t my_steps = recoX.size();

   if ( my_steps<2 ) return -1.0;

   Int_t check0 = GetRecoTracks( recoX, recoY, recoZ );

   if ( check0!=0 ) return -1.0;

   seg_size = steps_size;

   Int_t check1 = GetSegTracks2( recoX, recoY, recoZ );

   if ( check1!=0 ) return -1.0;

   Int_t seg_steps = segx.size();

   if ( seg_steps<2 ) return -1;

   Int_t seg_steps0 = seg_steps-1;

   Double_t recoL = segL.at(seg_steps0);

   if ( seg_steps<2 || recoL<minLength || recoL>maxLength ) return -1;

   Double_t mean = 666.0; Double_t rms = 666.0; Double_t rmse = 666.0;

   nmeas = 0; Double_t max1=-999.0; Double_t min1=+999.0;

   for ( Int_t j=0; j<n_steps; j++ )
     {
       Double_t trial = steps.at( j );

       GetDeltaThetaRMS( mean, rms, rmse, trial );

       // cout << mean << ",  " << rms << ", " << rmse << ", " << trial << endl;

       xmeas[ nmeas ] = trial;

       ymeas[ nmeas ] = rms;

       eymeas[ nmeas ] = sqrt( pow( rmse, 2.0 ) + pow( 0.05*rms, 2.0 ) ); // <--- conservative syst. error to fix chi^{2} behaviour !!!

       // ymeas[ nmeas ] = 10.0; eymeas[ nmeas ] = 1.0; // <--- for debugging !

       if ( min1>ymeas[ nmeas ] ) min1=ymeas[ nmeas ];

       if ( max1<ymeas[ nmeas ] ) max1=ymeas[ nmeas ];

       nmeas++;

     }

   gr_meas = new TGraphErrors( nmeas, xmeas, ymeas, 0, eymeas );

   gr_meas->SetTitle( "(#Delta#theta)_{rms} versus material thickness; Material thickness in cm; (#Delta#theta)_{rms} in mrad" );

   gr_meas->SetLineColor( kBlack ); gr_meas->SetMarkerColor( kBlack ); gr_meas->SetMarkerStyle( 20 ); gr_meas->SetMarkerSize( 1.2 );

   gr_meas->GetXaxis()->SetLimits( ( steps.at( 0 )-steps.at( 0 ) ), ( steps.at( n_steps-1 )+steps.at( 0 ) ) );

   gr_meas->SetMinimum( 0.0 );

   gr_meas->SetMaximum( 1.80*max1 );

   // c1->cd();

   // gr_meas->Draw( "APEZ" );

   // c1->Update();

   // c1->WaitPrimitive();

   ROOT::Minuit2::Minuit2Minimizer *mP = new ROOT::Minuit2::Minuit2Minimizer( );

   ROOT::Math::Functor FCA( &my_mcs_chi2, 2 );

   mP->SetFunction( FCA );


   mP->SetLimitedVariable( 1, "#delta#theta", 0.0, 1.0, 0.0, 45.0 );

   mP->SetMaxFunctionCalls( 1.E9 );

   mP->SetMaxIterations( 1.E9 );

   mP->SetTolerance( 0.01 );

   mP->SetStrategy( 2 );

   mP->SetErrorDef( 1.0 );

   bool mstatus = mP->Minimize();

   mP->Hesse();

   const double *pars = mP->X();

   const double *erpars = mP->Errors();

   Double_t deltap = ( recoL*kcal )/2.0;

   p_mcs = pars[0]+deltap;

   p_mcs_e = erpars[0];

   if ( mstatus ) p = p_mcs;

   else p = -1.0;

   chi2 = mP->MinValue();

   delete mP;

   return p;

	}

}
}
#endif
