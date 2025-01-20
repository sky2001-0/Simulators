#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "TGraph2D.h"

#include "../../lib/include/src/globals.h"
#include "../../lib/include/src/Range.h"
#include "../../lib/RootSupport/RootSupport.h"
#include "../../lib/Simulators.h"

#include "Model.h"

using namespace simulators_support;
using namespace rs::utils;
using namespace DemkovOsherovModel;


void analyze_searched(const std::filesystem::path treefile_path)
{
	CCheck();
	SetStyle();
	IgnoreWarning();
	CreateCustomPalette();

	const std::string treefile_name = treefile_path.stem();

	const double dipole_corr = (
		[treefile_name] () -> double {
			const std::string tmp = treefile_name.substr(4);
			if (tmp == "2") {
				return 1.;
			} else if (tmp == "3_0") {
				return std::sqrt(1. / 3.);
			} else if (tmp == "3_1") {
				return std::sqrt(1. / 2.);
			} else {
				throw std::runtime_error("Invalid name.");
			}
		} ()
	);

	auto f = rs::file::Open(treefile_path);
 	auto th = rs::TreeHelper(std::move(rs::file::GetObj<TTree>("tree", f.get())));

	const Int_t n = th.GetEntries();

	auto g0 = std::make_unique<TGraph2D>(n);

	auto g1 = rs::graph::Create(
		n, nullptr, "prob_adiabatic", "adiabatic", "max_prob"
	);
	{
		rs::graph::SetLimitY(g1.get(), {0., 1.});
	}

	auto g2 = rs::graph::Create(
		n, nullptr, "prob_rabi", "rabi (GHz)", "max_prob"
	);
	{
		rs::graph::SetLimitY(g2.get(), {0., 1.});
	}

	auto g3 = rs::graph::Create(
		n, nullptr, "time_adiabatic", "adiabatic", "max_time (ns)"
	);

	auto g4 = rs::graph::Create(
		n, nullptr, "time_rabi", "rabi (GHz)", "max_time (ns)"
	);

	double adiabatic_max;
	double rabi_rate_max;
	double max_prob_max = 0.;
	double max_time_max;

	for (auto entry : csp::Range<Long64_t>(n)) {
		th.GetEntry(entry);

		const double adiabatic = std::get<double>(th.cget("adiabatic/D"));
		const double rabi_rate = std::get<double>(th.cget("rabi_rate/D"));
		const double max_prob = std::get<double>(th.cget("max_prob/D"));
		const double max_time = std::get<double>(th.cget("max_time_ns/D"));

		g0->SetPoint(entry, adiabatic, rabi_rate, max_prob);
		g1->SetPoint(entry, adiabatic, max_prob);
		g2->SetPoint(entry, rabi_rate, max_prob);
		g3->SetPoint(entry, adiabatic, max_time);
		g4->SetPoint(entry, rabi_rate, max_time);
		if (max_prob > max_prob_max) {
			max_prob_max = max_prob;
			adiabatic_max = adiabatic;
			rabi_rate_max = rabi_rate;
			max_time_max = max_time;
		}
	}

	std::cout
		<< "Max Probability : "
		<< max_prob_max
		<< std::endl
		<< "Adiabatic : "
		<< adiabatic_max * dipole_corr
		<< std::endl
		<< "Max time : "
		<< max_time_max
		<< " ns"
		<< std::endl
		<< "Chirp Rate : "
		<< (
			std::pow(rabi_rate_max * ps::split_02 / adiabatic_max, 2) / 2.
			/ u::tau / u::GHz * u::ns
		)
		<< " GHz / ns"
		<< std::endl
		<< "Rabi Freq : "
		<< rabi_rate_max * dipole_corr * ps::split_02 / u::tau / u::GHz
		<< " GHz"
		<< std::endl;

	{
		g0->Scale(dipole_corr, "X");
		g0->Scale(dipole_corr * ps::split_02 / u::tau / u::GHz, "Y");

		g0->SetTitle("max_prob");
		g0->SetMinimum(0.);
    g0->SetMaximum(1.);

		auto c = std::make_unique<TCanvas>("prob", "prob");
		c->SetLogy();

		g0->Draw("TRI2");
    g0->GetXaxis()->SetTitle("adiabatic");
    g0->GetYaxis()->SetTitle("Rabi Freq (GHz)");
    g0->GetZaxis()->SetTitle("max_prob");
		c->SaveAs(
      (result_dirpath / (treefile_name + "prob.png")).c_str()
    );
	}
	{
		g1->Scale(dipole_corr, "X");
		rs::draw::FastSaveToFile(
			g1.get(), result_dirpath / (treefile_name + g1->GetName() + ".png"), "AP"
		);
	}
	{
    auto c = std::make_unique<TCanvas>(g2->GetName(), g2->GetTitle());
		c->SetLogx();
		g2->Scale(dipole_corr * ps::split_02 / u::tau / u::GHz, "X");
    g2->Draw("AP");
		c->SaveAs(
      (result_dirpath / (treefile_name + g2->GetName() + ".png")).c_str()
    );
	}
	{
		g3->Scale(dipole_corr, "X");
		rs::draw::FastSaveToFile(
			g3.get(), result_dirpath / (treefile_name + g3->GetName() + ".png"), "AP"
		);
	}
	{
    auto c = std::make_unique<TCanvas>(g4->GetName(), g4->GetTitle());
		c->SetLogx();
		g4->Scale(dipole_corr * ps::split_02 / u::tau / u::GHz, "X");
    g4->Draw("AP");
		c->SaveAs(
      (result_dirpath / (treefile_name + g4->GetName() + ".png")).c_str()
    );
	}
}
