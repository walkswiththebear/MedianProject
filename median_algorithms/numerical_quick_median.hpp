//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef NUMERICAL_QUICK_MEDIAN_10_06_2015_HPP
#define NUMERICAL_QUICK_MEDIAN_10_06_2015_HPP

/*
 * Median algorithm based on quickselect, using numerical pivoting.
 *
 * The algorithm works for numerical data only. This means that for a variable x
 * whose type is the value type of the sequence, the function static_cast<double>(x)
 * must define and order-preserving embedding from the sequence values into the reals.
 *
 * Arguments are begin and end iterators to the data. These must be random access
 * iterators. The return value is a double. For an uneven number of elements, it is
 * the median converted to a double. For an even number of elements, it is the mean
 * of the left and right endpoint of the median interval.
 */

#include "detail/numerical_quick_median_detail.hpp"
#include "detail/pivot_calculators/standard_numerical_pivot.hpp"
#include "detail/pivot_calculators/uniform_distribution_pivot.hpp"
#include "detail/pivot_calculators/normal_distribution_pivot.hpp"
#include "detail/pivot_calculators/exponential_distribution_pivot.hpp"
#include "detail/no_op_median_performance_stats.hpp"

namespace median_project
{

/**
 * Function numerical_quick_median
 * ===============================
 *
 * Median algorithm for numerical data. Use this algorithm when there is no
 * information or conjecture regarding the distribution of the data. This algorithm will
 * experience the worst case performance (O(n^2)) for the sequence \{2^(-i)\}_{i\in N).
 *
 */
template <typename RandomAccessIterator>
double numerical_quick_median(RandomAccessIterator begin, RandomAccessIterator end)
{
    no_op_median_performance_stats performance_stats;
    read_only_numerical_quick_median_detail::standard_numerical_pivot pivot_calculator;
    return numerical_quick_median_detail::numerical_quick_median_internal(
        begin, end, pivot_calculator, performance_stats);
}

/**
 * Function numerical_quick_median_for_uniform_distributions
 * =========================================================
 *
 * Median algorithm for uniformly distributed numerical data.
 *
 */
template <typename RandomAccessIterator>
double numerical_quick_median_for_uniform_distributions(RandomAccessIterator begin, RandomAccessIterator end)
{
    no_op_median_performance_stats performance_stats;
    read_only_numerical_quick_median_detail::uniform_distribution_pivot pivot_calculator;
    return numerical_quick_median_detail::numerical_quick_median_internal(
        begin, end, pivot_calculator, performance_stats);
}

/**
 * Function numerical_quick_median_for_normal_distributions
 *=========================================================
 *
 * Median algorithm for normally distributed numerical data.
 *
 */
template <typename RandomAccessIterator>
double numerical_quick_median_for_normal_distrubtions(RandomAccessIterator begin, RandomAccessIterator end)
{
    no_op_median_performance_stats performance_stats;
    read_only_numerical_quick_median_detail::normal_distribution_pivot pivot_calculator;
    return numerical_quick_median_detail::numerical_quick_median_internal(
        begin, end, pivot_calculator, performance_stats);
}

/**
 * Function numerical_quick_median_for_exponential_distributions
 * =============================================================
 *
 * Median algorithm for exponentially distributed numerical data.
 *
 */
template <typename RandomAccessIterator>
double numerical_quick_median_for_exponential_distrubtions(RandomAccessIterator begin, RandomAccessIterator end)
{
    no_op_median_performance_stats performance_stats;
    read_only_numerical_quick_median_detail::exponential_distribution_pivot pivot_calculator;
    return numerical_quick_median_detail::numerical_quick_median_internal(
        begin, end, pivot_calculator, performance_stats);
}

} // end namespace median_project

#endif // NUMERICAL_QUICK_MEDIAN_10_06_2015_HPP
