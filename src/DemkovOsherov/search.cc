#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "TF1.h"
#include "TTree.h"

#include "../../lib/Simulators.hh"
#include "../../lib/include/globals.hh"
#include "../../lib/include/Support.hh"
#include "../../lib/RootSupport/RootSupport.hh"

#include "Model.hh"

using namespace simulators_support;
using namespace rs::utils;
using namespace DemkovOsherovModel;

constexpr std::size_t size_of_model = 2;
const std::size_t first_state = 0;
const int magnetic_qnumber = 1;
const bool with_decay = true;

const auto adiabatic_range = csp::Range<double>(0.5, 2.9, 36, false);
const std::vector<double> rabi_rate_range = {
  0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10., 20.
};



class TreeObj
{
private:
  const std::size_t size_;
  const std::unique_ptr<TTree> tree_;

  double adiabatic_;
  double rabi_rate_;

  double max_prob_;
  double max_time_;
  double correctly_;

public:
  TreeObj(const std::size_t size)
  : size_(size),
    tree_(std::make_unique<TTree>("tree", "title"))
  {
    tree_->Branch("adiabatic", &adiabatic_, "adiabatic/D");
    tree_->Branch("rabi_rate", &rabi_rate_, "rabi_rate/D");
    tree_->Branch("max_prob", &max_prob_, "max_prob/D");
    tree_->Branch("max_time_ns", &max_time_, "max_time_ns/D");
  }

  ~TreeObj()
  {
    auto file = rs::file::Create(
      result_dirpath / ("tree" + std::to_string(size_of_model) + ".root")
    );
    tree_->Write();
  }

  void Simulate(const double adiabatic, const double rabi_rate)
  {
    adiabatic_ = adiabatic;
    rabi_rate_ = rabi_rate;
    max_prob_ = 0.;

    auto sim = ps_ortho1S2P::Create<size_of_model>(
      adiabatic_, rabi_rate_, magnetic_qnumber, with_decay
    );
    auto dmat = sim.InitialState();
    double time_n = -sim.get_inf();

    for (std::size_t i = 0; i < sim.get_size(); ++i) {
      sim.RK4(dmat, time_n);
      if (std::abs(dmat.getf(1, 1)) > max_prob_) {
        max_prob_ = std::abs(dmat.getf(1, 1));
        max_time_ = time_n;
      }
    }
    max_time_ /= ps_ortho1S2P::Calc_norml(adiabatic, rabi_rate);
    max_time_ /= u::ns;

    correctly_ = std::abs(dmat.trace());

    tree_->Fill();
  }
};


void search()
{
  CCheck();
  rs::utils::SetStyle();
  rs::utils::IgnoreWarning();

  auto to = TreeObj(size_of_model);
  csp::ProgressBar<std::size_t> pb(
    adiabatic_range.size() * rabi_rate_range.size()
  );
  for (auto adiabatic : adiabatic_range) {
    for (auto rabi_rate : rabi_rate_range) {
      to.Simulate(adiabatic, rabi_rate);
      ++pb;
    }
  }
}
