#ifndef DO_MODEL_HH
#define DO_MODEL_HH

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

#include "../../lib/include/globals.hh"
#include "../../lib/include/Matrix.hh"



namespace DemkovOsherovModel
{
/*
  Calculation
    2 i \partial_t = 2 v t \ketbra{0}{0} + 2 E / \hbar \ketbra{j}{j} + \Omega \ketbra{0}{j} + \Omega^* \ketbra{j}{0}
*/
template <std::size_t size_>
class Simulator
{
  static_assert(size_ > 1, "size_ must be greater than 1.");

public:
  using Mat = csp::math::Matrix<std::complex<double>, size_, size_>;
  using Arr = std::array<double, size_ - 1>;

  static constexpr double how_inf_ = 1e2;
  static constexpr double precision_ = 1e5;

private:
  const std::complex<double> factor_ = 2. * u::i;
  Mat h_mat_; // = H / (2 i hbar)
  const Arr decays_;

  const double time_transition_;
  const double time_inf_;
  const double time_step_;

  std::complex<double> hget(const std::size_t row, const std::size_t col) const
  {
    return h_mat_.cgetf(row, col) * factor_;
  }

  double Calc_transition_time(
    const std::array<std::complex<double>, size_ - 1>& rabis
  )
  {
    return std::max(
      1.,
      std::abs(*std::max_element(
        rabis.begin(),
        rabis.end(),
        [](const auto& elem1, const auto& elem2) {
          return std::abs(elem1) < std::abs(elem2);
        }
      ))
    );
  }

  Mat Step(const Mat& d_mat, const double time) noexcept
  {
    h_mat_.getf(0, 0) = time / factor_;
    Mat result = h_mat_.commute(d_mat);

    for (std::size_t i = 1; i < size_; ++i) {
      const double& decay = decays_[i - 1];
      result.getf(i, 0) -= decay / 2. * d_mat.cgetf(i, 0);
      result.getf(0, i) -= decay / 2. * d_mat.cgetf(0, i);
      const std::complex<double> tmp = decay * d_mat.cgetf(i, i);
      result.getf(0, 0) += tmp;
      result.getf(i, i) -= tmp;
    }

    return result;
  }

public:
  Simulator() = delete;

  Simulator(
    const std::array<double, size_ - 1>& energies,
    const std::array<std::complex<double>, size_ - 1>& rabis,
    const std::array<double, size_ - 1>& decays
  )
  : h_mat_(0.),
    decays_(decays),
    time_transition_(Calc_transition_time(rabis)),
    time_inf_(time_transition_ * how_inf_),
    time_step_(time_transition_ / precision_)
  {
    for (std::size_t i = 1; i < size_; ++i) {
      h_mat_.getf(0, i) = rabis[i - 1] / factor_;
      h_mat_.getf(i, 0) = std::conj(rabis[i - 1]) / factor_;
      h_mat_.getf(i, i) = energies[i - 1] / factor_;
    }

    const double time_transition = std::max(
      1.,
      std::abs(*std::max_element(
        rabis.begin(), rabis.end(),
        [](const auto& elem1, const auto& elem2) {
          return std::abs(elem1) < std::abs(elem2);
        }
      ))
    );
  }

  ~Simulator() noexcept = default;

  Simulator(const Simulator& rh) = delete;

  Simulator(Simulator&& rh) = delete;

  Simulator& operator=(const Simulator& rh) = delete;

  Simulator& operator=(Simulator&& rh) = delete;


  double get_transition() const noexcept
  {
    return time_transition_;
  }

  double get_inf() const noexcept
  {
    return time_inf_;
  }

  std::size_t get_size() const noexcept
  {
    return 2 * how_inf_ * precision_;
  }

  double get_first_transition_time_n() const noexcept
  {
    return std::abs(hget(1, 1));
  }

  double get_theorical_probability(const std::size_t i) const
  {
    if (i >= size_) {
      throw std::out_of_range("put of range.");
    }

    else if (i == 0) {
      double tmp = 1.;
      for (std::size_t j = 1; j < size_; ++j) {
        tmp *= std::exp(-u::pi * std::pow(std::abs(hget(0, j)), 2));
      }
      return tmp;
    }

    else {
      double tmp = 1.;
      for (std::size_t j = 1; j < i; ++j) {
        tmp *= std::exp(-u::pi * std::pow(std::abs(hget(0, j)), 2));
      }
      tmp *= 1. - std::exp(-u::pi * std::pow(std::abs(hget(0, i)), 2));
      return tmp;
    }
  }


  void RK4(Mat& dmat, double& time) noexcept
  {
    const Mat k1 = Step(dmat, time) * (time_step_ / 2.);
    const Mat k2 = Step(dmat + k1, time + time_step_ / 2.) * (time_step_ / 2.);
    const Mat k3 = Step(dmat + k2, time + time_step_ / 2.) * time_step_;
    const Mat k4 = Step(dmat + k3, time + time_step_);

    for (std::size_t i = 0; i < size_ * size_; ++i) {
      dmat.getf(i) += (
        k1.cgetf(i) / 3.
        + k2.cgetf(i) * 2. / 3.
        + k3.cgetf(i) / 3.
        + k4.cgetf(i) / 6. * time_step_
      );
    }

    time += time_step_;
  }


  static Mat InitialState(const std::size_t i = 0) noexcept
  {
    Mat initial(0.);
    initial.getf(i, i) = 1.;
    return initial;
  }
};


namespace ps_ortho1S2P
{
  constexpr double split_02 = (
    ps::ortho2P2_binding_energy - ps::ortho2P0_binding_energy
  ) / u::hbar;

  constexpr double split_12 = (
    ps::ortho2P2_binding_energy - ps::ortho2P1_binding_energy
  ) / u::hbar;

  const double gamma0 = (
    std::pow(2., 17)
    / std::pow(3., 10)
    / std::pow(u::c, 4)
    * std::pow(ps::trans_angfreq, 3)
    * std::pow(u::hbar / u::m_e, 2)
    / u::alpha
  );


  double Calc_norml(const double adiabatic, const double rabi_rate) noexcept
  {
    return rabi_rate * split_02 / adiabatic;
  }


  template <std::size_t size>
  Simulator<size> Create(
    const double adiabatic,
    const double rabi_rate,
    const int magnetic_qnumber,
    const bool with_decay
  )
  {
    const double norml = Calc_norml(adiabatic, rabi_rate);
    const double gamma = with_decay ? gamma0 / norml : 0.;

    if constexpr (size == 2) {
      if (std::abs(magnetic_qnumber) == 1) {
        return {{0.}, {adiabatic}, {gamma}};
      }

    } else if constexpr (size == 3) {
      if (magnetic_qnumber == 0) {
        return {
          {0., 2. * split_02 / norml},
          {adiabatic * 2. / 3., adiabatic * 4. / 3.},
          {4. / 3. * gamma, gamma}
        };

      } else if (std::abs(magnetic_qnumber) == 1) {
        return {
          {0., 2. * split_12 / norml},
          {adiabatic, adiabatic},
          {3. / 2. * gamma, 3. / 2. * gamma}
        };
      }
    }

    throw std::runtime_error("not implemented.");
  }
}

}



#endif // DO_MODEL_HH
