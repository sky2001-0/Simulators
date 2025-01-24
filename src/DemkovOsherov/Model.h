#ifndef DEMKOV_OSHEROV_MODEL_H
#define DEMKOV_OSHEROV_MODEL_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <utility>

#include "../../lib/include/src/globals.h"
#include "../../lib/include/src/Matrix.h"
#include "../../lib/RootSupport/RootSupport.h"



namespace ps
{
  constexpr double split_02 = (
    ps::ortho2P2_binding_energy - ps::ortho2P0_binding_energy
  ) / u::hbar;

  constexpr double split_12 = (
    ps::ortho2P2_binding_energy - ps::ortho2P1_binding_energy
  ) / u::hbar;

  const double gamma0 = (
    std::pow(2., 19)
    / std::pow(3., 11)
    / std::pow(u::c, 4)
    * std::pow(ps::trans_angfreq, 3)
    * std::pow(u::hbar / u::m_e, 2)
    / u::alpha
  );

  const double dipole0 = (
    std::pow(2., 8.5) * u::e * u::hbar
    / std::pow(3., 5) / u::alpha / u::m_e / u::c
  );
}


namespace DemkovOsherovModel
{


// Calculation
// i \partial_t = t \ketbra{0}{0} + energies[j - 1] \ketbra{j}{j} + (rabis[j - 1] \ketbra{0}{j} + std::conj(rabis[j - 1]) \ketbra{j}{0}) / 2
template <std::size_t SizeOfModel>
class Simulator
{
  static_assert(SizeOfModel > 1, "SizeOfModel must be greater than 1.");


public:
  using Mat = csp::Matrix<std::complex<double>, SizeOfModel, SizeOfModel>;
  using Arr = std::array<double, SizeOfModel - 1>;


  static Simulator<SizeOfModel> PsOrtho1S2P(
    const double chirp_rate,
    const double rabi,
    const int magnetic_qnumber,
    const bool with_decay
  )
  {
    const double norml = std::sqrt(std::abs(chirp_rate));
    const double rabi_n = rabi / norml;
    const double gamma_n = with_decay ? ps::gamma0 / norml : 0.;

    double sign;
    if (chirp_rate > 0) {
      sign = 1.;
    } else if (chirp_rate < 0) {
      sign = -1.;
    } else {
      throw std::runtime_error("chirp_rate is zero.");
    }

    if constexpr (SizeOfModel == 2) {
      if (std::abs(magnetic_qnumber) == 1) {
        return {{0.}, {rabi_n}, {gamma_n}, sign};
      }

    } else if constexpr (SizeOfModel == 3) {
      if (magnetic_qnumber == 0) {
        return {
          {0., 2. * ps::split_02 / norml},
          {rabi_n / std::sqrt(3.), rabi_n * std::sqrt(2. / 3.)},
          {gamma_n, gamma_n},
          sign
        };

      } else if (std::abs(magnetic_qnumber) == 1) {
        return {
          {0., 2. * ps::split_12 / norml},
          {rabi_n / std::sqrt(2.), rabi_n / std::sqrt(2.)},
          {gamma_n, gamma_n},
          sign
        };
      }
    }

    throw std::runtime_error("not implemented.");
  }


  static constexpr double kHowInf = 1e2;
  static constexpr double kPrecision = 1e5;


private:
  // hmat_ = -u::i * H / u::hbar
  Mat hmat_;
  const Arr decays_;
  const double time_transition_;
  const double time_step_;
  const double sign_;


public:
  Simulator() = delete;

  // energies[j] < energies[j + 1]
  Simulator(
    const std::array<double, SizeOfModel - 1>& energies,
    const std::array<std::complex<double>, SizeOfModel - 1>& rabis,
    const std::array<double, SizeOfModel - 1>& decays,
    const double sign
  )
  : hmat_(std::move(CreateMat(energies, rabis))),
    decays_(decays),
    time_transition_(CalcTransitionTime(rabis)),
    time_step_(time_transition_ / kPrecision),
    sign_(sign)
  {
  }

  ~Simulator() = default;

  Simulator(const Simulator& rh) = delete;

  Simulator(Simulator&& rh) = delete;

  Simulator& operator=(const Simulator& rh) = delete;

  Simulator& operator=(Simulator&& rh) = delete;


  double cget_transition() const noexcept
  {
    return time_transition_;
  }

  double time_start() const noexcept
  {
    return -time_transition_ * kHowInf;
  }

  std::size_t size() const noexcept
  {
    return 2 * kHowInf * kPrecision;
  }


  double CalcTheoricalProbability(const std::size_t i) const
  {
    auto func = [this] (const std::size_t level) {
      return 2. * std::pow(std::abs(this->hmat_.cgetf(0, level)), 2);
    };

    if (i >= SizeOfModel) {
      throw std::out_of_range("put of range.");
    }

    else if (i == 0) {
      double tmp = 1.;
      for (std::size_t j = 1; j < SizeOfModel; ++j) {
        tmp *= std::exp(-u::pi * func(j));
      }
      return tmp;
    }

    else {
      double tmp = 1.;
      if (sign_ > 0.) {
        for (std::size_t j = 1; j < i; ++j) {
          tmp *= std::exp(-u::pi * func(j));
        }
        tmp *= 1. - std::exp(-u::pi * func(i));
        return tmp;
      } else {
        for (std::size_t j = SizeOfModel - 1; i < j; --j) {
          tmp *= std::exp(-u::pi * func(j));
        }
        tmp *= 1. - std::exp(-u::pi * func(i));
        return tmp;
      }
    }
  }

  void RK4(Mat& dmat, double& time) noexcept
  {
    const Mat k1 = Step(dmat, time) * (time_step_ / 2.);
    const Mat k2 = Step(dmat + k1, time + time_step_ / 2.) * (time_step_ / 2.);
    const Mat k3 = Step(dmat + k2, time + time_step_ / 2.) * time_step_;
    const Mat k4 = Step(dmat + k3, time + time_step_) * (time_step_ / 6.);

    dmat += k1 * (1. / 3.);
    dmat += k2 * (2. / 3.);
    dmat += k3 * (1. / 3.);
    dmat += k4;

    time += time_step_;
  }

private:
  Mat Step(const Mat& dmat, const double time) noexcept
  {
    hmat_.getf(0, 0) = -sign_ * u::i * time;
    Mat result = hmat_.commute(dmat);

    for (std::size_t i = 1; i < SizeOfModel; ++i) {
      const double& decay = decays_[i - 1];
      result.getf(i, 0) -= decay / 2. * dmat.cgetf(i, 0);
      result.getf(0, i) -= decay / 2. * dmat.cgetf(0, i);
      const auto tmp = decay * dmat.cgetf(i, i);
      result.getf(0, 0) += tmp;
      result.getf(i, i) -= tmp;
    }

    return result;
  }

  static Mat CreateMat(
    const std::array<double, SizeOfModel - 1>& energies,
    const std::array<std::complex<double>, SizeOfModel - 1>& rabis
  )
  {
    Mat mat = {};
    for (std::size_t i = 1; i < SizeOfModel; ++i) {
      mat.getf(0, i) = rabis[i - 1] / 2.;
      mat.getf(i, 0) = std::conj(rabis[i - 1]) / 2.;
      mat.getf(i, i) = energies[i - 1];
    }

    mat *= -u::i;
    return mat;
  }

  static double CalcTransitionTime(
    const std::array<std::complex<double>, SizeOfModel - 1>& rabis
  )
  {
    const double rabi_max = std::abs(
      *std::max_element(
        rabis.begin(), rabis.end(),
        [](const auto& elem1, const auto& elem2) {
          return std::abs(elem1) < std::abs(elem2);
        }
      )
    );

    return std::max(1., rabi_max);
  }

public:
  static Mat InitialState(const std::size_t i = 0) noexcept
  {
    Mat mat = {};
    mat.getf(i, i) = 1.;
    return mat;
  }
};


template<std::size_t SizeOfModel>
class ResultData
{
private:
  const double norml_;
  const Simulator<SizeOfModel>& sim_;

  const std::filesystem::path result_dirpath_;

  int index_;

  std::unique_ptr<TGraph> g_cor_;
  std::unique_ptr<TGraph> g_pure_;
  std::array<std::unique_ptr<TGraph>, SizeOfModel> g_arr_;

  const bool graph_cut_;


public:
  ResultData() = delete;

  ResultData(
    const Simulator<SizeOfModel>& sim,
    const double chirp_rate,
    const std::filesystem::path result_dirpath,
    const bool graph_cut,
    const int mag = 1
  )
  : sim_(sim),
    norml_(std::sqrt(chirp_rate)),
    result_dirpath_(result_dirpath),
    graph_cut_(graph_cut)
  {
    index_ = 0;

    g_cor_ = std::move(GraphInitialize(
      "g_correctly", "Total Probability", mag
    ));
    g_pure_ = std::move(GraphInitialize("g_purely", "Purity", mag));
    for (std::size_t i = 0; i < SizeOfModel; ++i) {
      g_arr_[i] = std::move(GraphInitialize(
        "g_excited_" + std::to_string(i),
        "the probability of excited state " + std::to_string(i),
        mag
      ));
    }
    g_arr_[0] = std::move(GraphInitialize(
      "g_ground",
      "the probability of ground state",
      mag
    ));
  }

  ~ResultData() = default;


  void Set(
    const typename Simulator<SizeOfModel>::Mat& dmat,
    const double time_n
  )
  {
    g_cor_->SetPoint(
      index_, to_time_ns(time_n), std::abs(dmat.trace())
    );
    g_pure_->SetPoint(
      index_, to_time_ns(time_n), std::abs((dmat * dmat).trace())
    );
    for (std::size_t i = 0; i < SizeOfModel; ++i) {
      g_arr_[i]->SetPoint(
        index_, to_time_ns(time_n), std::abs(dmat.cgetf(i, i))
      );
    }
    ++index_;
  }

  void Save()
  {
    {
      rs::graph::SetLimitY(g_cor_.get(), {0.999, 1.001});
      rs::draw::FastSaveToFile(
        g_cor_.get(),
        result_dirpath_ / (std::string(g_cor_->GetName()) + ".png"),
        "AL"
      );
    }

    {
      rs::graph::SetLimitY(g_pure_.get(), {0., 1.01});
      rs::draw::FastSaveToFile(
        g_pure_.get(),
        result_dirpath_ / (std::string(g_pure_->GetName()) + ".png"),
        "AL"
      );
    }

    for (std::size_t i = 0; i < SizeOfModel; ++i) {
      const std::unique_ptr<TGraph>& g = g_arr_[i];

      if (graph_cut_) {
        rs::graph::SetLimit(
          g_arr_[i].get(),
          {
            to_time_ns(-30. * sim_.cget_transition()),
            to_time_ns(30. * sim_.cget_transition())
          },
          {0., 1.}
        );
      }

      auto c = std::make_unique<TCanvas>(g->GetName(), g->GetTitle());
      g->Draw("AL");

      rs::draw::DrawLineHorizontal(
        c.get(), sim_.CalcTheoricalProbability(i), kRed
      );

      c->SaveAs(
        (result_dirpath_ / (std::string(g->GetName()) + ".png")).c_str()
      );
    };
  }


private:
  double to_time_ns(const double time)
  {
    return time / norml_ / u::ns;
  }

  std::unique_ptr<TGraph> GraphInitialize(
    const std::string name, const std::string title, const int mag
  )
  {
    return std::move(rs::graph::Create<TGraph>(
      sim_.size() * mag,
      name.c_str(),
      ("time evolution of " + title).c_str(),
      "time (ns)",
      "probability"
    ));
  }
};


void CreateCustomPalette()
{
	const int num_color = 99;
	double stops[num_color];
	double red[num_color];
	double green[num_color];
	double blue[num_color];

	for (int i = 0; i < num_color; ++i) {
			stops[i] = static_cast<double>(i) / (num_color - 1);
			const double c = std::pow(stops[i], 6);
			red[i] = std::pow(c, 2);
			green[i] = 4. * (c - std::pow(c, 2));
			blue[i] = std::pow(1. - c, 2);
	}

	TColor::CreateGradientColorTable(num_color, stops, red, green, blue, 255);
	gStyle->SetNumberContours(num_color);
}


}



#endif // DEMKOV_OSHEROV_MODEL_H
