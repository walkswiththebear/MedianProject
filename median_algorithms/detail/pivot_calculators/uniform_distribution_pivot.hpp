//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_UNIFORM_DISTRIBUTION_PIVOT_09_27_2015_HPP
#define TMB_UNIFORM_DISTRIBUTION_PIVOT_09_27_2015_HPP

#include "pivot_functor_base.hpp"

/*
 * Pivot calculator for uniformly distributed data.
 */

namespace median_project
{
namespace read_only_numerical_quick_median_detail
{

/**
 * Pivot calculator for uniform distributions: use the midpoint of sequence min
 * and max.
 */
class uniform_distribution_pivot : public pivot_functor_base
{
  public:

    double operator()(double median_lower_bound, double median_upper_bound) const override
    {
        return median_lower_bound / 2.0 + median_upper_bound / 2.0;
    }
};
} // end namespace read_only_numerical_quick_median_detail
} // end namespace median_project

#endif // TMB_UNIFORM_DISTRIBUTION_PIVOT_09_27_2015_HPP
