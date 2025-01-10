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



void compare_searched()
{
	CCheck();
	SetStyle();
	IgnoreWarning();

	auto f0 = rs::file::Open("result/tree3_0.root");
 	auto th0 = rs::TreeHelper(std::move(
		rs::file::GetObj<TTree>("tree", f0.get())
	));

	auto f1 = rs::file::Open("result/tree3_1.root");
 	auto th1 = rs::TreeHelper(std::move(
		rs::file::GetObj<TTree>("tree", f1.get())
	));

	const Int_t n = th0.GetEntries();

	for (auto entry : csp::Range<Long64_t>(n)) {
		th0.GetEntry(entry);
		th1.GetEntry(entry);

		const double adiabatic0 = std::get<double>(th0.cget("adiabatic/D"));
		const double rabi_rate0 = std::get<double>(th0.cget("rabi_rate/D"));
		const double max_prob0 = std::get<double>(th0.cget("max_prob/D"));
		const double max_time0 = std::get<double>(th0.cget("max_time_ns/D"));

		const double adiabatic1 = std::get<double>(th1.cget("adiabatic/D"));
		const double rabi_rate1 = std::get<double>(th1.cget("rabi_rate/D"));
		const double max_prob1 = std::get<double>(th1.cget("max_prob/D"));
		const double max_time1 = std::get<double>(th1.cget("max_time_ns/D"));

		if (
			max_prob0 < 0.95
			|| max_prob1 < 0.95
			|| (
				max_prob0 * std::exp(ps::gamma0 * (max_time0 - max_time1) * u::ns)
				+ 2. * max_prob1
			 ) / 3. < 0.96
			|| rabi_rate0 * ps::split_02 / std::sqrt(2.) / u::tau / u::GHz > 50.
		) {
			continue;
		}

		const double rabi = rabi_rate0 * ps::split_02;
		const double dipole = (
			std::pow(2., 8.5) * u::e * u::hbar
			/ std::pow(3., 5) / u::alpha / u::m_e / u::c
		);
		std::cout
			<< "u (GHz/ps), F (kW/mm^2), Omega (GHz) = \t"
			<< std::pow(rabi / adiabatic0, 2) / 2.  / u::tau / u::GHz * u::ps
			<< "\t"
			<< (
				u::c * u::epsilon_0 * std::pow(u::hbar * rabi / dipole, 2)
				/ u::kW * u::mm * u::mm
			)
			<< "\t"
			<< rabi / std::sqrt(2.) / u::tau / u::GHz
			<< std::endl
			<< "Value : "
			<< (
				max_prob0 * std::exp(ps::gamma0 * (max_time0 - max_time1) * u::ns)
				+ 2. * max_prob1
			) / 3.
			<< "(p0, p1, t0 (ps), t1 (ps)) = "
			<< max_prob0
			<< "\t"
			<< max_prob1
			<< "\t"
			<< max_time0 * u::ns / u::ps
			<< "\t"
			<< max_time1 * u::ns / u::ps
			<< std::endl
			<< std::endl;
	}
}
