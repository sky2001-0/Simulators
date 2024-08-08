#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>

#include "TF1.h"

#include "../../lib/Simulators.hh"
#include "../../lib/include/globals.hh"
#include "../../lib/include/Support.hh"
#include "../../lib/RootSupport/RootSupport.hh"

#include "DO_Model.hh"

using namespace simulators_support;
using namespace rs::utils;

constexpr std::size_t size_of_model = 3;
using State = std::array<std::complex<double>, size_of_model>;

const int j = 1;

const csp::Range adiabatic_range = {0.9, 0.1, 8};
const csp::Range rabi_rate_range = {0.2, 0.1, 40};


struct Results
{
  const std::unique_ptr<TGraph> g1_;
  const std::unique_ptr<TGraph> g2_;
  int graph_index_;

public:
  Results()
  : g1_(std::move(rs::graph::Create<TGraph>(
    adiabatic_range.get_size() * rabi_rate_range.get_size(),
    nullptr,
    "the time until the probability of excited state 2 reaches (1 - exp(-#pi #kappa^{2}))^{2}",
    "2 #kappa #sqrt{1 + (2 #varepsilon)^{2}}",
    "#sqrt{2 v}  t_{re}"
  ))),
    g2_(std::move(rs::graph::Create<TGraph>(
    adiabatic_range.get_size() * rabi_rate_range.get_size(),
    nullptr,
    "the time until the probability of excited state 2 reaches (1 - exp(-#pi #kappa^{2}))^{2}",
    "index",
    "#sqrt{2 v} t_{re}"
  ))),
    graph_index_(0)
  {
  }

  ~Results()
  {
    std::filesystem::create_directory(result_dirpath);

    g1_->Fit("pol1", "Q");
    rs::draw::SetLimit(
      g1_,
      {
        0.,
        4. * adiabatic_range.get_back() * rabi_rate_range.get_back() * 1.3
      }
    );
    g1_->SetMinimum(0.);
    rs::draw::FastSaveToFile(g1_, result_dirpath / "propto.png", "AP", 'L');
    rs::draw::FastSaveToFile(g2_, result_dirpath / "time_n.png", "AP");
  }

  void Set(const double x, const double y)
  {
    g1_->SetPoint(graph_index_, x, y);
    g2_->SetPoint(graph_index_, graph_index_, y);
    graph_index_++;
  }
};




void search()
{
  CCheck();

  static_assert(size_of_model == 3, "size_of_model == 3");

  rs::draw::SetStyle();
  rs::utils::IgnoreWarning();

  Results results;

  auto create_model = [] (const double adiabatic, const double rabi_rate)
    -> DO_Simulator<size_of_model>
  {
    constexpr double split_up = (
      ps::ps_ortho2P2_binding_energy - ps::ps_ortho2P1_binding_energy
    ) / u::hbar;

    constexpr double split_down = (
      ps::ps_ortho2P1_binding_energy - ps::ps_ortho2P0_binding_energy
    ) / u::hbar;

    if (j == 0) {
      return {
        {split_up, -split_down},
        std::pow(rabi_rate * split_up / adiabatic, 2) / 2.,
        {rabi_rate * 2. * split_up, rabi_rate * split_up}
      };

    } else if (std::abs(j) == 1) {
      return {
        split_up,
        std::pow(rabi_rate * split_up / adiabatic, 2) / 2.,
        rabi_rate * split_up
      };

    } else {
      throw std::runtime_error("j must be -1, 0, or 1.");
    }
  };

  auto simulate = [&results, &create_model] (
    const double adiabatic,
    const double rabi_rate
  )
  {
    auto sim = create_model(adiabatic, rabi_rate);

    State state = sim.InitialState();
    double time_n = -sim.get_inf();

    for (auto i : csp::Range(sim.get_size())) {
      sim.RK4(state, time_n);
      if (std::norm(state[1]) > std::pow(sim.get_theorical_probability(1), 2)) {
        break;
      }
    }
    results.Set(
      4. * adiabatic * std::sqrt(1. + std::pow(rabi_rate, 2)),
      time_n - sim.get_first_transition_time_n()
    );
  };

  csp::ProgressBar<std::size_t> pb(
    adiabatic_range.get_size() * rabi_rate_range.get_size()
  );
  for (auto adiabatic : adiabatic_range) {
    for (auto rabi_rate : rabi_rate_range) {
      simulate(adiabatic, rabi_rate);
      ++pb;
    }
  }
}
