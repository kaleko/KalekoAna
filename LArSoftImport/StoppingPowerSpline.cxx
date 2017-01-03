#ifndef LARLITE_STOPPINGPOWERSPLINE_CXX
#define LARLITE_STOPPINGPOWERSPLINE_CXX

#include "StoppingPowerSpline.h"

namespace larlite {

	StoppingPowerSpline::StoppingPowerSpline() {


		// // Initialize the spline with hardcoded (fit) variables.
		// // See StoppingPowerSpline header file to see where these came from.
		// Float_t muon_momentum_GEV[49] = {0.0006,    0.0009,    0.0012,    0.0016,    0.0024,    0.0035,
		//                                  0.0048,    0.0061,    0.0075,    0.0092,    0.0126,    0.0173,    0.0228,    0.0322,
		//                                  0.0456,    0.0645,    0.1007,    0.1526,    0.2314,    0.3508,    0.5137,    0.7788,
		//                                  1.2655,    1.9185,    3.1153,    5.2407,    9.1212,   18.2491,   30.6789,   42.1547,
		//                                  57.5945,   84.3404,  115.4033,  163.2351,  222.6906,  274.1947,  349.5230,  461.2661,
		//                                  652.4494,  926.8821, 1270.1491, 1796.5898, 2824.0367, 3858.3852, 5650.1559,
		//                                  8273.9852, 12120.0000, 19100.0000, 29980.0000
		//                                 };

		// //This stopping power is in units of MeV cm^2 / g
		// Float_t stopping_power[49] = {142.067, 178.768, 202.739, 211.389, 198.548, 157.786, 115.341, 87.911,
		//                               68.418, 51.069, 33.629, 21.238, 14.281,  8.471,  5.131,  3.240,  2.107,  1.640,  1.447,  1.417,
		//                               1.447,  1.508,  1.573,  1.638,  1.744,  1.819,  1.937,  2.062,  2.195,  2.291,  2.439,  2.653,
		//                               2.945,  3.338,  3.863,  4.289,  4.964,  5.991,  7.698, 10.117, 13.358, 18.274, 27.178, 37.964,
		//                               54.151, 78.869, 114.871, 182.200, 282.532
		//                              };

		// Initialize the spline with hardcoded (fit) variables.
		// See StoppingPowerSpline header file to see where these came from.
		Float_t muon_momentum_GEV[80] = {
			0.000581,
			0.000822,
			0.001000,
			0.001000,
			0.002000,
			0.003000,
			0.005000,
			0.006000,
			0.007000,
			0.009000,
			0.010000,
			0.011000,
			0.012000,
			0.013000,
			0.013000,
			0.014000,
			0.015000,
			0.016000,
			0.018000,
			0.020000,
			0.021000,
			0.024000,
			0.027000,
			0.030000,
			0.034000,
			0.038000,
			0.043000,
			0.045000,
			0.048000,
			0.052000,
			0.055000,
			0.058000,
			0.061000,
			0.064000,
			0.068000,
			0.072000,
			0.078000,
			0.086000,
			0.095000,
			0.106000,
			0.118000,
			0.130000,
			0.143000,
			0.155000,
			0.160000,
			0.172000,
			0.189000,
			0.217000,
			0.234000,
			0.275000,
			0.329000,
			0.482000,
			0.730000,
			1.186000,
			1.798000,
			2.918000,
			4.908000,
			8.538000,
			17.074000,
			28.694000,
			39.418000,
			53.844000,
			78.827000,
			107.836000,
			152.494000,
			207.992000,
			256.059000,
			326.350000,
			430.601000,
			608.927000,
			864.841000,
			1184.870000,
			1675.559000,
			2632.961000,
			3596.543000,
			5265.307000,
			7708.372000,
			11285.001000,
			17783.637000,
			27903.365000
		};

		//This stopping power is in units of MeV cm^2 / g
		Float_t stopping_power[80] = {
			127.192,
			158.684,
			179.119,
			186.470,
			175.552,
			140.713,
			104.070,
			80.128,
			62.948,
			47.501,
			37.988,
			35.140,
			31.771,
			29.518,
			26.626,
			24.425,
			22.144,
			20.412,
			18.216,
			15.509,
			13.930,
			11.411,
			9.954,
			8.425,
			7.219,
			6.082,
			5.200,
			4.871,
			4.490,
			4.057,
			3.780,
			3.590,
			3.341,
			3.137,
			2.960,
			2.833,
			2.642,
			2.383,
			2.208,
			2.027,
			1.908,
			1.804,
			1.734,
			1.676,
			1.643,
			1.622,
			1.576,
			1.537,
			1.534,
			1.503,
			1.506,
			1.537,
			1.600,
			1.666,
			1.733,
			1.840,
			1.916,
			2.035,
			2.162,
			2.296,
			2.392,
			2.541,
			2.755,
			3.047,
			3.437,
			3.957,
			4.375,
			5.037,
			6.036,
			7.684,
			9.996,
			13.063,
			17.662,
			25.881,
			35.705,
			50.257,
			72.179,
			103.662,
			161.615,
			246.543
		};


		stopping_power_vs_momentum = new TGraph(80, muon_momentum_GEV, stopping_power);
		// stopping_power_vs_momentum_spline3 = new TSpline3("spline", stopping_power_vs_momentum);

		_argon_density = 1.396; //g/cm^3
		_mu_mass_GEV = 0.105658;
		_mu_mass_GEV_squared = _mu_mass_GEV * _mu_mass_GEV;

	}

	double StoppingPowerSpline::GetMudEdx(const double muon_momentum_GEV) {
		// return stopping_power_vs_momentum_spline3->Eval(muon_momentum_GEV) * _argon_density;
		return stopping_power_vs_momentum->Eval(muon_momentum_GEV) * _argon_density;
	}


	double StoppingPowerSpline::Getp(const double initial_muon_momentum_GEV, const double length_travelled) {

		size_t n_steps = 10;

		double step_size = length_travelled / n_steps;

		double current_p = initial_muon_momentum_GEV;
		double current_E = p_to_E(initial_muon_momentum_GEV);

		for (size_t i = 0; i < n_steps; ++i) {
			// std::cout << "At this point, length traversed is " << i*step_size
			//           << " of " << length_travelled << std::endl;
			// std::cout << " current_p is " << current_p << std::endl;
			// std::cout << " current_E is " << current_E << std::endl;
			// double this_dedx_GEV = stopping_power_vs_momentum_spline3->Eval(current_p) * _argon_density / 1000.;
			double this_dedx_GEV = stopping_power_vs_momentum->Eval(current_p) * _argon_density / 1000.;
			// std::cout << " this_dedx_GEV is " << this_dedx_GEV << std::endl;
			current_E -= (this_dedx_GEV * step_size);
			if ( current_E <= _mu_mass_GEV ) {
				// std::cout<<"WARNING: current_E less than mu mass. it is "<<current_E<<std::endl;
				return 0.;
			}
			current_p = E_to_p(current_E);
		}

		return current_p;
	}

}
#endif
