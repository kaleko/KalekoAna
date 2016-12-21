#ifndef LARLITE_STOPPINGPOWERSPLINE_CXX
#define LARLITE_STOPPINGPOWERSPLINE_CXX

#include "StoppingPowerSpline.h"

namespace larlite {

	StoppingPowerSpline::StoppingPowerSpline() {

		// Initialize the spline with hardcoded (fit) variables.
		// See StoppingPowerSpline header file to see where these came from.
		Float_t muon_momentum_GEV[49] = {0.0006,    0.0009,    0.0012,    0.0016,    0.0024,    0.0035,
		                                 0.0048,    0.0061,    0.0075,    0.0092,    0.0126,    0.0173,    0.0228,    0.0322,
		                                 0.0456,    0.0645,    0.1007,    0.1526,    0.2314,    0.3508,    0.5137,    0.7788,
		                                 1.2655,    1.9185,    3.1153,    5.2407,    9.1212,   18.2491,   30.6789,   42.1547,
		                                 57.5945,   84.3404,  115.4033,  163.2351,  222.6906,  274.1947,  349.5230,  461.2661,
		                                 652.4494,  926.8821, 1270.1491, 1796.5898, 2824.0367, 3858.3852, 5650.1559,
		                                 8273.9852, 12120.0000, 19100.0000, 29980.0000
		                                };

		//This stopping power is in units of MeV cm^2 / g
		Float_t stopping_power[49] = {142.067, 178.768, 202.739, 211.389, 198.548, 157.786, 115.341, 87.911,
		                              68.418, 51.069, 33.629, 21.238, 14.281,  8.471,  5.131,  3.240,  2.107,  1.640,  1.447,  1.417,
		                              1.447,  1.508,  1.573,  1.638,  1.744,  1.819,  1.937,  2.062,  2.195,  2.291,  2.439,  2.653,
		                              2.945,  3.338,  3.863,  4.289,  4.964,  5.991,  7.698, 10.117, 13.358, 18.274, 27.178, 37.964,
		                              54.151, 78.869, 114.871, 182.200, 282.532
		                             };


		stopping_power_vs_momentum = new TGraph(49, muon_momentum_GEV, stopping_power);
		stopping_power_vs_momentum_spline3 = new TSpline3("spline", stopping_power_vs_momentum);

	}

	double StoppingPowerSpline::GetMudEdx(const double muon_momentum_GEV) {
		double argon_density = 1.396; //g/cm^3
		return stopping_power_vs_momentum_spline3->Eval(muon_momentum_GEV) * argon_density;
	}


}
#endif
