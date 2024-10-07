#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>

#include "../../lib/Simulators.hh"
#include "../../lib/include/globals.hh"
#include "../../lib/include/Support.hh"
#include "../../lib/RootSupport/RootSupport.hh"

#include "Model.hh"

using namespace simulators_support;
using namespace rs::utils;
using namespace DemkovOsherovModel;

constexpr std::size_t size_of_model = 3;
const std::size_t first_state = 1;
const int magnetic_qnumber = 1;
const bool with_decay = true;

const double adiabatic = 1.2;
const double rabi_rate = 25;
const double omega_hw = u::tau * 500. * u::GHz;



struct ResultData
{
private:
  const double norml_;
  const Simulator<size_of_model>& sim_;

  int index_;

  std::unique_ptr<TGraph> g_cor_;
  std::unique_ptr<TGraph> g_pure_;
  std::array<std::unique_ptr<TGraph>, size_of_model> g_arr_;

  std::unique_ptr<TGraph> GraphInitialize(
    const std::string name, const std::string title
  )
  {
    return std::move(rs::graph::Create<TGraph>(
      sim_.get_size(),
      name.c_str(),
      ("time evolution of " + title).c_str(),
      "time (ns)",
      "probability"
    ));
  }

  double to_time_ns(const double time)
  {
    return time / norml_ / u::ns;
  }

public:
  ResultData() = delete;

  ResultData(const Simulator<size_of_model>& sim)
  : sim_(sim),
    norml_(ps_ortho1S2P::Calc_norml(adiabatic, rabi_rate))
  {
    index_ = 0;

    g_cor_ = std::move(GraphInitialize("g_correctly", "Total Probability"));
    g_pure_ = std::move(GraphInitialize("g_purely", "Purity"));
    for (auto i : csp::Range<std::size_t>(1, size_of_model)) {
      g_arr_[i] = std::move(GraphInitialize(
        "g_excited_" + std::to_string(i),
        "the probability of excited state " + std::to_string(i)
      ));
    }
    g_arr_[0] = std::move(GraphInitialize(
      "g_ground",
      "the probability of ground state"
    ));
  }

  void Set(const Simulator<size_of_model>::Mat& dmat, const double time)
  {
    g_cor_->SetPoint(index_, to_time_ns(time), std::abs(dmat.trace()));
    if (!with_decay) {
      g_pure_->SetPoint(index_, to_time_ns(time), std::abs((dmat * dmat).trace()));
    }

    for (std::size_t i = 0; i < size_of_model; ++i) {
      g_arr_[i]->SetPoint(
        index_, to_time_ns(time), std::abs(dmat.cgetf(i, i))
      );
    }
    ++index_;
  }

  void Save()
  {
    {
      rs::graph::SetLimitY(g_cor_, {0.999, 1.001});
      rs::draw::FastSaveToFile(
        g_cor_, result_dirpath / (std::string(g_cor_->GetName()) + ".png"), "AL"
      );
    }

    if (!with_decay) {
      rs::graph::SetLimitY(g_pure_, {0.999, 1.001});
      rs::draw::FastSaveToFile(
        g_pure_,
        result_dirpath / (std::string(g_pure_->GetName()) + ".png"), "AL"
      );
    }

    for (auto i : csp::Range<std::size_t>(size_of_model)) {
      const std::unique_ptr<TGraph>& g = g_arr_[i];

      if (size_of_model == 2) {
        rs::graph::SetLimit(
          g_arr_[i],
          {
            to_time_ns(-30. * sim_.get_transition()),
            to_time_ns(30. * sim_.get_transition())
          },
          {0., 1.}
        );
      } else {
        rs::graph::SetLimitY(g, {0., 1.});
      }

      auto c = std::make_unique<TCanvas>(g->GetName(), g->GetTitle());
      g->Draw("AL");

      rs::draw::DrawLineHorizontal(c, sim_.get_theorical_probability(i), kRed);
      // rs::draw::DrawLineVertical(c, 8. * adiabatic / norml_ / u::ns, kRed);
      rs::draw::DrawLineVertical(
        c, 2. * omega_hw / std::pow(norml_, 2) / u::ns, kBlue
      );

      c->SaveAs(
        (result_dirpath / (std::string(g->GetName()) + ".png")).c_str()
      );
    };
  }
};



void draw()
{
  CCheck();
  SetStyle();
  IgnoreWarning();

  auto sim = ps_ortho1S2P::Create<size_of_model>(
    adiabatic, rabi_rate, magnetic_qnumber, with_decay
  );

  ResultData result(sim);
  auto dmat = sim.InitialState(first_state);
  double time_n = -sim.get_inf();

  csp::ProgressBar<std::size_t> pb(sim.get_size());

  for (auto i : csp::Range(sim.get_size())) {
    result.Set(dmat, time_n);
    sim.RK4(dmat, time_n);
    ++pb;
  }

  result.Save();
}
