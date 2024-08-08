#ifndef DO_MODEL_HH
#define DO_MODEL_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numeric>
#include <type_traits>

#include "../../lib/include/Support.hh"



namespace f
{
  template <std::size_t size>
  std::array<std::complex<double>, size> Add(
    const std::array<std::complex<double>, size>& initial_state,
    const std::array<std::complex<double>, size>& add_state,
    const double factor
  )
  {
    std::array<std::complex<double>, size> result;
    std::transform(
      initial_state.begin(),
      initial_state.end(),
      add_state.begin(),
      result.begin(),
      [factor](const std::complex<double>& a, const std::complex<double>& b) {
        return a + b * factor;
      }
    );
    return result;
  }


  template <typename DividableType, std::size_t size>
  std::array<DividableType, size> Divide(
    const std::array<DividableType, size>& arr,
    const double factor
  )
  {
    static_assert(
      std::is_constructible_v<
        DividableType,
        decltype(std::declval<DividableType>() / std::declval<double>())
      >,
      "Template parameter 'DividableType' must be dividable by double."
    );

    std::array<DividableType, size> result;
    std::transform(
      arr.begin(),
      arr.end(),
      result.begin(),
      [factor](const DividableType& a) {
        return a / factor;
      }
    );
    return result;
  }


  template <std::size_t size>
  double Norm(const std::array<std::complex<double>, size>& state)
  {
    return std::accumulate(
      state.begin(), state.end(), 0.,
      [](const double acc, const std::complex<double>& elem) {
        return acc + std::norm(elem);
      }
    );
  }
}



template <std::size_t size_>
class DO_Simulator
{
  static_assert(size_ > 1, "size_ must be bigger than 1.");
  using State = std::array<std::complex<double>, size_>;
  using DiagOfMat = std::array<double, size_ - 1>;
  using OffOfMat = std::array<std::complex<double>, size_ - 1>;

public:
  static constexpr double how_inf_ = 70.;
  static constexpr double precision_ = 1e5;

  static State InitialState() noexcept
  {
    State initial = {0.};
    initial[size_ - 1] = 1.;
    return initial;
  }

private:
  const double normalizer_;
  const DiagOfMat diag_;
  const OffOfMat off_diag_;

  const double transition_time_n_;
  const double time_n_inf_;
  const double time_n_step_;
  const std::complex<double> diff_factor_;
  const std::size_t step_size_;

  static DiagOfMat ExpandDiagOfMat(const double omega_split) noexcept
  {
    DiagOfMat diag;
    for (auto i : csp::Range(size_ - 1)) {
      diag[i] = (size_ / 2. - i - 1.) * omega_split;
    }
    return diag;
  }

  static OffOfMat ExpandOffOfMat(
    const std::complex<double> rabi_angfreq
  ) noexcept
  {
    OffOfMat off_diag;
    off_diag.fill(rabi_angfreq);
    return off_diag;
  }

  double CalcTransitionTimeN() const noexcept
  {
    return std::sqrt(u::tau) * std::max(
      1.,
      std::abs(*std::max_element(
        off_diag_.begin(),
        off_diag_.end(),
        [](const auto& a, const auto& b) {
          return std::abs(a) < std::abs(b);
        }
      )) / u::pi
    );
  }

public:
  DO_Simulator() = delete;

  DO_Simulator(
    const DiagOfMat&& omega_split,
    const double chirp_rate,
    const OffOfMat&& rabi_angfreq
  )
  : normalizer_(std::sqrt(2. * std::abs(chirp_rate))),
    diag_(f::Divide(omega_split, normalizer_)),
    off_diag_(f::Divide(rabi_angfreq, normalizer_)),
    transition_time_n_(CalcTransitionTimeN()),
    time_n_inf_(transition_time_n_ * how_inf_),
    time_n_step_(transition_time_n_ / precision_),
    diff_factor_(-0.5 * u::i * time_n_step_),
    step_size_(std::ceil(2. * time_n_inf_ / time_n_step_))
  {
    if (chirp_rate == 0.) {
      throw std::invalid_argument("chirp rate must not be 0.");
    }
  }

  DO_Simulator(
    const double omega_split,
    const double chirp_rate,
    const std::complex<double> rabi_angfreq
  )
  : DO_Simulator<size_>(
    ExpandDiagOfMat(omega_split),
    chirp_rate,
    ExpandOffOfMat(rabi_angfreq)
  )
  {
  }

  ~DO_Simulator() noexcept = default;

  double get_adiabatic_parameter(const std::size_t i = 0) const
  {
    return std::abs(off_diag_.at(i));
  }

  double get_omega_split(const std::size_t i = 0) const
  {
    return diag_.at(i) * normalizer_;
  }

  double get_inf() const noexcept
  {
    return time_n_inf_;
  }

  int get_size() const noexcept
  {
    return step_size_;
  }

  double get_chirp_rate() const noexcept
  {
    return std::pow(normalizer_, 2) / 2.;
  }

  double get_step() const noexcept
  {
    return time_n_step_;
  }

  double get_first_transition_time_n() const noexcept
  {
    return 2. * diag_[size_ - 2];
  }

  double get_theorical_probability(const std::size_t i) const
  {
    if (i >= size_) {
      throw std::out_of_range("put of range.");
    }

    else if (i == size_ - 1) {
      double tmp = 1.;
      for (auto j : csp::Range(size_ - 1)) {
        tmp *= std::exp(
          -u::pi * std::pow(std::abs(off_diag_[j]), 2)
        );
      }
      return tmp;
    }

    else {
      double tmp = 1.;
      for (auto j : csp::Range(i + 1, size_ - 1)) {
        tmp *= std::exp(
          -u::pi * std::pow(std::abs(off_diag_[j]), 2)
        );
      }
      tmp *= (
        1.
        - std::exp(
          -u::pi * std::pow(std::abs(off_diag_[i]), 2)
        )
      );
      return tmp;
    }
  }

private:
  State Diff(const State& state, const double time_n) const noexcept
  {
    State diff_state;
    std::complex<double> tmp = 0.;
    for (auto i : csp::Range(size_ - 1)) {
      diff_state[i] = diff_factor_ * (
        diag_[i] * state[i] + std::conj(off_diag_[i]) * state[size_ - 1]
      );
      tmp += off_diag_[i] * state[i];
    }
    diff_state[size_ - 1] = diff_factor_ * (
      tmp + time_n * state[size_ - 1]
    );
    return diff_state;
  }

public:
  double to_time(const double time_n) const noexcept
  {
    return time_n / normalizer_;
  }


  void RK4(State& state, double& time_n) const noexcept
  {
    const State k1 = Diff(state, time_n);
    const State k2 = Diff(
      f::Add(state, k1, 1. / 2.), time_n + time_n_step_ / 2.
    );
    const State k3 = Diff(
      f::Add(state, k2, 1. / 2.), time_n + time_n_step_ / 2.
    );
    const State k4 = Diff(
      f::Add(state, k3, 1.), time_n + time_n_step_
    );

    for (auto i : csp::Range(size_)) {
      state[i] += (k1[i] + (k2[i] + k3[i]) * 2. + k4[i]) / 6.;
    }

    time_n += time_n_step_;
  }
};



#endif // DO_MODEL_HH
