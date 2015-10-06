//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_NORMAL_DISTRIBUTION_PIVOT_09_27_2015_HPP
#define TMB_NORMAL_DISTRIBUTION_PIVOT_09_27_2015_HPP

#include <boost/math/special_functions/erf.hpp>
#include "distribution_specific_pivot_base.hpp"

/*
 * Pivot calculator for normally distributed data.
 */

namespace median_project
{
namespace read_only_numerical_quick_median_detail
{

/**
 * Pivot calculator for normal distributions.
 */
class normal_distribution_pivot : public distribution_specific_pivot_base
{
  private:
    virtual double cdf(double x) const override
    {
        return .5 * (1.0 + boost::math::erf((x - m_total_sequence_mean) / (sqrt(2.0) * m_total_sequence_std_dev)));
    }

    virtual double quantile(double x) const override
    {
        return m_total_sequence_mean + m_total_sequence_std_dev * sqrt(2.0) * boost::math::erf_inv(2.0 * x - 1.0);
    }
};
} // end namespace read_only_numerical_quick_median_detail
} // end namespace median_project

#endif // TMB_NORMAL_DISTRIBUTION_PIVOT_09_27_2015_HPP
