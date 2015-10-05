//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_DISTRIBUTION_SPECIFIC_PIVOT_BASE_10_03_2015_HPP
#define TMB_DISTRIBUTION_SPECIFIC_PIVOT_BASE_10_03_2015_HPP

#include <boost/math/special_functions/erf.hpp>
#include "pivot_functor_base.hpp"

/*
 * Base class for pivot calculators that are specific to a distribution.
 */

namespace median_project
{
namespace read_only_numerical_quick_median_detail
{

/**
 * Base class for pivot calculators that are specific to a distribution.
 */
class distribution_specific_pivot_base : public pivot_functor_base
{
  public:
    virtual double operator()(double median_lower_bound,
                              double median_upper_bound,
                              int num_elements_less_than_median_lower_bound,
                              int num_elements_greater_than_median_upper_bound) const override
    {
        double size_of_interval_of_median_candidates =
            static_cast<double>(m_total_sequence_length - num_elements_less_than_median_lower_bound -
                                num_elements_greater_than_median_upper_bound);
                                
        double num_desired_elements_less_than_or_equal_to_pivot =
            static_cast<double>(m_total_sequence_length) / 2.0 -
            static_cast<double>(num_elements_less_than_median_lower_bound);
            
        double lambda = num_desired_elements_less_than_or_equal_to_pivot /
                        size_of_interval_of_median_candidates;

        double pivot = quantile((1.0 - lambda) * cdf(median_lower_bound) + lambda * cdf(median_upper_bound));
        return std::min(median_upper_bound, std::max(pivot, median_lower_bound));
    }

  private:
    virtual double cdf(double x) const = 0;
    virtual double quantile(double x) const = 0;
};
} // end namespace read_only_numerical_quick_median_detail
} // end namespace median_project

#endif // TMB_DISTRIBUTION_SPECIFIC_PIVOT_BASE_10_03_2015_HPP
