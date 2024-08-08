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

#include "DO_Model.hh"

using namespace simulators_support;
using namespace rs::utils;

constexpr std::size_t size_of_model = 2;
using State = std::array<std::complex<double>, size_of_model>;

const int j = 1;

const double adiabatic = 1.4;
const double rabi_rate = 6.;

const double omega_hw = u::tau * 500. * u::GHz;



void SimulateOnePass(
  const DO_Simulator<size_of_model>& sim,
  State state = {}
);
void SimulateTwoPass(const DO_Simulator<size_of_model>& sim);



void draw()
{
  CCheck();
  SetStyle();
  IgnoreWarning();

  auto func = [] () -> DO_Simulator<size_of_model>
  {
    const double split_up = (
      ps::ps_ortho2P2_binding_energy - ps::ps_ortho2P1_binding_energy
    ) / u::hbar;

    const double split_down = (
      ps::ps_ortho2P1_binding_energy - ps::ps_ortho2P0_binding_energy
    ) / u::hbar;

    if (j == 0) {
      if (size_of_model != 3) {
        throw std::runtime_error("'j' = 0 needs 'size_of_model' != 3");
      }

      // return {
      //   {split_up, -split_down},
      //   std::pow(rabi_rate * split_up / adiabatic, 2) / 2.,
      //   {rabi_rate * 2. * split_up, rabi_rate * split_up}
      // };

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

  auto sim = func();

  SimulateOnePass(sim);
  // SimulateTwoPass(sim);
}




struct ResultData
{
  using State = State;

private:
  int graph_index_;

  std::unique_ptr<TGraph> g_cor_;
  std::array<std::unique_ptr<TGraph>, size_of_model> g_arr_;

  const DO_Simulator<size_of_model>& sim_;

public:
  ResultData() = delete;

  ResultData(const DO_Simulator<size_of_model>& sim)
  : sim_(sim)
  {
    graph_index_ = 0;

    const int n = sim.get_size();
    auto initialize = [n] (const std::string title) {
      return std::move(rs::graph::Create<TGraph>(
        n,
        nullptr,
        ("time evolution of " + title).c_str(),
        "time (ns)",
        "probability"
      ));
    };

    g_cor_ = std::move(initialize("Total Probability"));

    if (size_of_model == 2) {
        g_arr_[0] = std::move(initialize(
          "the probability of excited state"
        ));
    } else {
      for (auto i : csp::Range(size_of_model - 1)) {
        g_arr_[i] = std::move(initialize(
          "the probability of excited state " + std::to_string(i + 1)
        ));
      }
    }
    g_arr_[size_of_model - 1] = std::move(initialize(
      "the probability of ground state"
    ));
  }

  std::array<std::unique_ptr<TGraph>, size_of_model>& get_g_arr_()
  {
    return g_arr_;
  }

  void Set(const State& state, const double time)
  {
    g_cor_->SetPoint(graph_index_, time / u::ns, f::Norm(state));

    for (auto i : csp::Range(size_of_model)) {
      g_arr_[i]->SetPoint(graph_index_, time / u::ns, std::norm(state[i]));
    }
    ++graph_index_;
  }

  void Save()
  {
    std::filesystem::create_directory(result_dirpath);

    rs::graph::SetLimitY(g_cor_, {0.999, 1.001});
    rs::draw::FastSaveToFile(g_cor_, result_dirpath / "g_correctly.png", "AL");

    auto func = [this](
      const std::size_t i,
      const std::string name
    ) {
      const std::unique_ptr<TGraph>& g = this->g_arr_[i];

      if (size_of_model == 2) {
        rs::graph::SetLimit(
          g,
          {
            sim_.to_time(-30. * sim_.get_step() * sim_.precision_) / u::ns,
            sim_.to_time(30. * sim_.get_step() * sim_.precision_) / u::ns
          },
          {0., 1.}
        );
      } else {
        rs::graph::SetLimitY(g, {0., 1.});
      }

      const std::filesystem::path filepath
        = result_dirpath / ("g_" + name + ".png");
      auto c = std::make_unique<TCanvas>(g->GetName(), g->GetTitle());
      g->Draw("AL");

      rs::draw::DrawLineHorizontal(
        c,
        this->sim_.get_theorical_probability(i),
        kRed
      );

      rs::draw::DrawLineVertical(
        c,
        std::sqrt(2 / sim_.get_chirp_rate())
          * sim_.get_adiabatic_parameter()
          / u::ns,
        kRed
      );

      // rs::draw::DrawLineVertical(
      //   c,
      //   omega_hw / sim_.get_chirp_rate() / u::ns,
      //   kBlue
      // );

      c->SaveAs(filepath.c_str());
    };

    for (auto i : csp::Range(size_of_model - 1)) {
      func(
        i,
        ("excited_" + std::to_string(i)).c_str()
      );
    }
    func(size_of_model - 1, "ground");
  }
};



void SimulateOnePass(
  const DO_Simulator<size_of_model>& sim,
  State state = {}
)
{
  ResultData result(sim);

  if (f::Norm(state) <= 0.1) {
    state = sim.InitialState();
  }
  double time_n = -sim.get_inf();

  csp::ProgressBar<std::size_t> pb(sim.get_size());

  for (auto i : csp::Range(sim.get_size())) {
    result.Set(state, sim.to_time(time_n));
    sim.RK4(state, time_n);
    ++pb;
  }

  result.Save();
}


void SimulateTwoPass(const DO_Simulator<size_of_model>& sim)
{
  ResultData result(sim);

  auto state = sim.InitialState();
  double time_n = -sim.get_inf();

  csp::ProgressBar<std::size_t> pb(sim.get_size() * 2);

  for (auto i : csp::Range(sim.get_size())) {
    result.Set(state, sim.to_time(time_n - sim.get_inf()));
    sim.RK4(state, time_n);
    ++pb;
    if (std::norm(state[size_of_model - 1]) < 0.1) {
      break;
    }
  }

  const double time_n_rel = time_n;
  time_n = -sim.get_inf();
  for (auto i : csp::Range(sim.get_size())) {
    result.Set(state, sim.to_time(time_n + time_n_rel));
    sim.RK4(state, time_n);
    ++pb;
  }

  result.Save();

  // auto f = rs::file::Create(result_dirpath / "test.root");
  // result.get_g_arr_()[size_of_model - 1]->Write();
}
