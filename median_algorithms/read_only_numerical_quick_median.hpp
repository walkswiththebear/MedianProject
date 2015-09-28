//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_READ_ONLY_NUMERICAL_QUICK_MEDIAN_08_24_2015_HPP
#define TMB_READ_ONLY_NUMERICAL_QUICK_MEDIAN_08_24_2015_HPP

/*
 * Non-modifying median algorithm in average N * log(N) time. This is appropriate
 * for situations where the data must not be modified, and there are not enough
 * resources to make a copy of the data.
 *
 * The algorithm works for numerical data only. This means that for a variable x
 * whose type is the value type of the sequence, the function static_cast<double>(x)
 * must define and order-preserving embedding from the sequence values into the reals.
 *
 * The difference between this algorithm and the corresponding general algorithm
 * read_only_quick_median is that this algorithm calculates the "pivot",
 * that is, the next element to be considered as a median candidate, from the current
 * lower and upper bound for the median. This means that whether or not the worst case
 * performance is hit does not depend on the ordering of the elements in the sequence,
 * but on their numerical distribution. This can be advantageous, mostly so when it is
 * known that the distribution is uniform. It can also be treacherous: the worst case
 * or near worst case happens when the distribution is rather one-sided, e.g., it is 
 * exponential. This means that hitting the worst case performance is not per se
 * unlikely, as it is with the random pivot choice. Therefore, using this algorithm
 * requires discretion in the sense that one must have knowledge or make conjectures
 * about the distribution of the elements of the sequence.
 *
 *
 * Arguments are begin and end iterators to the data. These do not  have to be better
 * than forward iterators. The return value is a double. For an uneven
 * number of elements, it is the median converted to a double. For an even number of
 * elements, it is the mean of the left and right endpoint of the median interval.
 *
 * A brief description of the algorithm can be found at lines 195 - 219 of the implementation
 * file "read_only_numerical_quick_median_detail.hpp".
 */

#include "detail/read_only_numerical_quick_median_detail.hpp"
#include "detail/pivot_calculators.hpp"
#include "detail/no_op_median_performance_stats.hpp"

namespace median_project
{

/**
 * Function read_only_numerical_quick_median
 * =========================================
 *
 * Non-modifying median algorithm for numerical data in average N * log(N) time.
 * Iterators need to be forward or better.
 *
 */
template <typename Iterator>
double read_only_numerical_quick_median(Iterator begin, Iterator end)
{
    no_op_median_performance_stats performance_stats;
    read_only_numerical_quick_median_detail::uniform_distribution_pivot pivot_calculator;
    return read_only_numerical_quick_median_detail::read_only_numerical_quick_median_internal(
        begin, end, pivot_calculator, performance_stats);
}
} // end namespace median_project

#endif // TMB_READ_ONLY_NUMERICAL_QUICK_MEDIAN_08_24_2015_HPP
