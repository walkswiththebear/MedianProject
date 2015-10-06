//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_EXPONENTIAL_DISTRIBUTION_PIVOT_09_27_2015_HPP
#define TMB_EXPONENTIAL_DISTRIBUTION_PIVOT_09_27_2015_HPP

#include <boost/math/special_functions/erf.hpp>
#include "distribution_specific_pivot_base.hpp"

/*
 * Pivot calculator for exponentially distributed data.
 */

namespace median_project
{
namespace read_only_numerical_quick_median_detail
{

/**
 * Pivot calculator for exponential distributions.
 */
class exponential_distribution_pivot : public distribution_specific_pivot_base
{
  private:
    virtual double cdf(double x) const override
    {
        double lambda = 1.0 / m_total_sequence_mean;
        return 1.0 - std::exp(-(lambda * x));
    }

    virtual double quantile(double x) const override
    {
        double lambda = 1.0 / m_total_sequence_mean;
        return  - std::log(1.0 - x) / lambda;
    }
};
} // end namespace read_only_numerical_quick_median_detail
} // end namespace median_project

#endif // TMB_EXPONENTIAL_DISTRIBUTION_PIVOT_09_27_2015_HPP
