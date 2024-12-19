#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <cmath>
#include <stdexcept>
#include <string>

#include "../../lib/include/src/globals.h"
#include "../../lib/include/src/Support.h"



namespace u
{
  constexpr double hydrogen_reduced_mass = u::m_e / (1. + u::m_e / u::m_p);
}


namespace ps
{
  constexpr double reduced_mass = u::m_e / 2.;
}


namespace Ionization
{
const std::string l_to_str(const int l)
{
  switch (l) {
    case 0 : {
      return "S";
    }
    case 1 : {
      return "P";
    }
    case 2 : {
      return "D";
    }
    case 3 : {
      return "F";
    }
    default : {
      throw std::runtime_error("Not implemented.");
    }
  }
}


class Calculator
{
private:
  const double a0_;
  const double n_;
  const double l_;


public:
  Calculator() = delete;

  explicit Calculator(const double reduced_mass, const int n, const int l)
  : a0_(u::hbar / u::alpha / reduced_mass / u::c),
    n_(n),
    l_(l)
  {
    if (n <= 0 || l < 0 || n <= l) {
      throw std::domain_error("no proper n or l.");
    }
  }

  ~Calculator() = default;

  Calculator(const Calculator& rh) = delete;

  Calculator(Calculator&& rh) = delete;

  Calculator& operator=(const Calculator& rh) = delete;

  Calculator& operator=(Calculator&& rh) = delete;


  csp::Range<double> GetRange()
  {
    const double binding_wavelength = (
      2. * u::tau * std::pow(n_, 2) / u::alpha * a0_
    );

    return {binding_wavelength / 6., binding_wavelength, 100, false};
  }

  double Calc(const double light_lambda)
  {
    const double coeff = (
      16. * std::pow(u::pi, 2) / 3.
      * std::pow(a0_, 3) / light_lambda / (2 * l_ + 1.)
    );
    const double factor = to_factor(light_lambda);

    if (l_ == 0) {
      if (n_ == 1) {
        return coeff * std::pow(func_m1_to_m0(factor), 2);
      }
      if (n_ == 2) {
        return coeff * std::pow(func_m2_to_m1(factor), 2);
      }
    } else {
      if (n_ == l_ + 1) {
        return coeff * (
          (l_ + 1) * std::pow(func_m1_to_m0(factor), 2)
          + l_ * std::pow(func_m1_to_m2(factor), 2)
        );
      }
      if (n_ == l_ + 2) {
        return coeff * (
          (l_ + 1) * std::pow(func_m2_to_m1(factor), 2)
          + l_ * std::pow(func_m2_to_m3(factor), 2)
        );
      }
    }
    throw std::runtime_error("not implemented for n, l");
  }


private:
  double to_factor(const double light_lambda)
  {
    /*
      factor := bohr_radius * wave_number_of_scattering_particle
    */

    const double rest = (
      2. * u::tau * a0_ / u::alpha / light_lambda - 1. / std::pow(n_, 2)
    );

    if (rest < 0.) {
      std::cout << light_lambda << std::endl;
      throw std::runtime_error("cannot be ionized.");
    }

    return std::sqrt(rest);
  }

  double func_m1_to_m0(const double factor)
  {
    double tmp = 1.;
    for (int s = 1; s <= n_; ++s) {
      tmp *= 1. + std::pow(s * factor, 2);
    }

    return (
      std::pow(n_, 2)
      * std::sqrt(u::pi / 2. / csp::math::Factorial(std::round(2. * n_ - 1.)))
      * 4. * std::pow(4. * n_, n_)
      * std::sqrt(tmp)
      / std::sqrt(1. - std::exp(-u::tau / factor))
      * std::exp(-2. / factor * std::atan(n_ * factor))
      / std::pow(1. + std::pow(n_ * factor, 2), n_ + 2)
    );
  }

  double func_m2_to_m1(const double factor)
  {
    return (
      0.5
      * std::sqrt(2. * n_ - 1.)
      * std::sqrt(1. + std::pow(n_ * factor, 2))
      * func_m1_to_m0(factor)
    );
  }

  double func_m1_to_m2(const double factor)
  {
    return (
      0.5 / n_
      * std::sqrt(1. + std::pow(n_ * factor, 2))
      / std::sqrt(1. + std::pow((n_ - 1.) * factor, 2))
      * func_m1_to_m0(factor)
    );
  }

  double func_m2_to_m3(const double factor)
  {
    return (
      (4. + (n_ - 1.) * (1. + std::pow(n_ * factor, 2))) / 2. / n_
      * std::sqrt(2. * n_ - 1.)
      / std::sqrt(1. + std::pow((n_ - 2.) * factor, 2))
      * func_m1_to_m2(factor)
    );
  }
};

}



#endif // CALCULATOR_H
