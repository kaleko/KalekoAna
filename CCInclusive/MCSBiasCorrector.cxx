#ifndef MCSBIASCORRECTOR_CXX
#define MCSBIASCORRECTOR_CXX

#include "MCSBiasCorrector.h"

namespace larlite {

	MCSBiasCorrector::MCSBiasCorrector() {

		Float_t lengths[25] =
		{	97.5,  112.5,  127.5,  142.5,  157.5,  172.5,  187.5,  202.5,  217.5,  232.5,
			247.5,  262.5,  277.5,  292.5,  307.5,  322.5,  337.5,  352.5,  367.5,  382.5,
			397.5,  412.5,  427.5,  442.5,  457.5
		};

		Float_t biases_MC[25] =
		{	-1.29, 0.18, 0.18, 0.16, 0.15, 0.14, 0.13, 0.13, 0.13, 0.12,
			0.12, 0.12, 0.11, 0.11, 0.1, 0.1, 0.1, 0.09, 0.08, 0.08, 0.08, 0.08, 0.07, 0.07
		};

		Float_t biases_data[25] =
		{	-1.27, 0.15, 0.14, 0.12, 0.12, 0.11, 0.1, 0.08, 0.09, 0.08,
			0.09, 0.09, 0.08, 0.08, 0.08, 0.1, 0.1, 0.11, 0.12, 0.1, 0.09, 0.09, 0.1, 0.11
		};

		biasgraph_MC = new TGraph(25, lengths, biases_MC);
		biasgraph_data = new TGraph(25, lengths, biases_data);

		biasspline_MC = new TSpline3("biasspline_MC", biasgraph_MC);
		biasspline_data = new TSpline3("biasspline_data", biasgraph_data);

	}

	double MCSBiasCorrector::GetBiasFactor(const double muon_range, bool data_true_MC_false) {

		return data_true_MC_false ? biasspline_data->Eval(muon_range) : biasspline_MC->Eval(muon_range);

	}


}
#endif
