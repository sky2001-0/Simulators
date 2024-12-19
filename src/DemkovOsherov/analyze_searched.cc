#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../../lib/include/src/globals.h"
#include "../../lib/include/src/Range.h"
#include "../../lib/RootSupport/RootSupport.h"
#include "../../lib/Simulators.h"

#include "Model.h"

using namespace simulators_support;
using namespace rs::utils;



void analyze_searched(const std::filesystem::path treefile_path)
{
	CCheck();
	SetStyle();
	IgnoreWarning();

	auto f = rs::file::Open(treefile_path);
 	auto th = rs::TreeHelper(std::move(rs::file::GetObj<TTree>("tree", f.get())));

	const Int_t n = th.GetEntries();

	auto g1 = rs::graph::Create(
		n, nullptr, "prob_adiabatic", "adiabatic", "max_prob"
	);
	rs::graph::SetLimitY(g1.get(), {0., 1.});
	auto g2 = rs::graph::Create(
		n, nullptr, "prob_rabi", "rabi (GHz)", "max_prob"
	);
	rs::graph::SetLimitY(g2.get(), {0., 1.});

	double adiabatic_max;
	double rabi_rate_max;
	double max_prob_max = 0.;

	for (auto entry : csp::Range<Long64_t>(th.GetEntries())) {
		th.GetEntry(entry);

		const double adiabatic = std::get<double>(th.cget("adiabatic/D"));
		const double rabi_rate = std::get<double>(th.cget("rabi_rate/D"));
		const double max_prob = std::get<double>(th.cget("max_prob/D"));
		const double max_time = std::get<double>(th.cget("max_time_ns/D")) * u::ns;

		g1->SetPoint(entry, adiabatic, max_prob);
		g2->SetPoint(entry, rabi_rate, max_prob);
		if (max_prob > max_prob_max) {
			max_prob_max = max_prob;
			adiabatic_max = adiabatic;
			rabi_rate_max = rabi_rate;
		}
	}

	std::cout << adiabatic_max << std::endl;
	std::cout << rabi_rate_max << std::endl;
	std::cout << max_prob_max << std::endl;

	const std::string treefile_name = treefile_path.stem();

	rs::draw::FastSaveToFile(
		g1.get(), result_dirpath / (treefile_name + g1->GetName() + ".png"), "AP"
	);
	{
    auto c = std::make_unique<TCanvas>(g2->GetName(), g2->GetTitle());
		c->SetLogx();
		g2->Scale(ps::split_02 / u::tau / u::GHz, "X");
    g2->Draw("AP");
		c->SaveAs(
      (result_dirpath / (treefile_name + g2->GetName() + ".png")).c_str()
    );
	}
}
