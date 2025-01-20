#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "TScatter.h"

#include "../../lib/include/src/globals.h"
#include "../../lib/include/src/Range.h"
#include "../../lib/RootSupport/RootSupport.h"
#include "../../lib/Simulators.h"

#include "Model.h"

using namespace simulators_support;
using namespace rs::utils;
using namespace DemkovOsherovModel;



void compare_searched()
{
	CCheck();
	SetStyle();
	IgnoreWarning();
	CreateCustomPalette();

	auto f0 = rs::file::Open("result/tree3_0.root");
 	auto th0 = rs::TreeHelper(std::move(
		rs::file::GetObj<TTree>("tree", f0.get())
	));

	auto f1 = rs::file::Open("result/tree3_1.root");
 	auto th1 = rs::TreeHelper(std::move(
		rs::file::GetObj<TTree>("tree", f1.get())
	));

	const Int_t n = th0.GetEntries();

  std::vector<double> x_arr(n);
  std::vector<double> y_arr(n);
  std::vector<double> c_arr(n);
  std::vector<double> size_arr(n);

	for (auto entry : csp::Range<Long64_t>(n)) {
		th0.GetEntry(entry);
		th1.GetEntry(entry);

		const double adiabatic0 = std::get<double>(th0.cget("adiabatic/D"));
		const double rabi_rate0 = std::get<double>(th0.cget("rabi_rate/D"));
		const double max_prob0 = std::get<double>(th0.cget("max_prob/D"));
		const double max_time0 = std::get<double>(th0.cget("max_time_ns/D"));
		const double max_prob1 = std::get<double>(th1.cget("max_prob/D"));
		const double max_time1 = std::get<double>(th1.cget("max_time_ns/D"));

		const double rabi = rabi_rate0 * ps::split_02;
		const double chirp_rate = std::pow(rabi / adiabatic0, 2) / 2.;
		const double gdd = 1. / chirp_rate;
		const double rabi1 = rabi / std::sqrt(2.);

		const double value = (
			max_prob0 * std::exp(ps::gamma0 * (max_time0 - max_time1) * u::ns)
			+ 2. * max_prob1
		) / 3.;

		x_arr[entry] = std::log10(gdd / u::fs / u::fs);
		y_arr[entry] = std::log10(rabi1 / u::tau / u::GHz);
		c_arr[entry] = value;
		size_arr[entry] = 1;

		if (
			value > 0.96
			&& max_prob0 > 0.95
			&& max_prob1 > 0.95
			&& gdd < 1e12 / u::tau / 555. * u::fs * u::fs
			&& rabi1 < 30 * u::tau * u::GHz
		) {
			std::cout
				<< gdd / u::fs / u::fs
				<< std::endl
				<< chirp_rate / u::tau / u::GHz * u::ps
				<< std::endl
				<< rabi1 / u::tau / u::GHz
				<< std::endl
				<< value
				<< std::endl
				<< max_prob0
				<< "\t"
				<< max_prob1
				<< std::endl
				<< max_time0
				<< "\t"
				<< max_time1
				<< std::endl
				<< u::tau * 555. * u::GHz / chirp_rate / u::ns
				<< std::endl
				<< std::endl;
		}
	}

	{
		auto g = std::make_unique<TScatter>(
			n, x_arr.data(), y_arr.data(), c_arr.data(), size_arr.data()
		);
		g->SetMaxMarkerSize(0.1);
		g->SetTitle(
			"value;log10( gdd (fs^2) );log10( rabi_freq 1 (2 pi GHz) );bar(w)"
		);
		auto c = std::make_unique<TCanvas>("value", "value");
		c->SetRightMargin(0.15);
		g->Draw("A");

		TLatex lx;
		lx.SetTextSize(0.04);
		lx.SetTextFont(42);
		lx.DrawLatex(15.8, 3.1, "#bar{w}");

		rs::draw::DrawLineHorizontal(c.get(), std::log10(6.16));
		rs::draw::DrawLineVertical(c.get(), std::log10(1e12 / u::tau / 555.));
		rs::draw::DrawLine(c.get(), {12.2, -0.8}, {5.1, 2.9});

		c->SaveAs((result_dirpath / "prob3_01.png").c_str());
	}
}
