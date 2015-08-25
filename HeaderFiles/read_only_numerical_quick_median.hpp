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
 * Arguments are begin and end iterators to the data. These do not  have to be better
 * than forward iterators. The return value is a double. For an uneven
 * number of elements, it is the median converted to a double. For an even number of
 * elements, it is the mean of the left and right endpoint of the median interval.
 *
 * A brief description of the algorithm can be found at lines 245 - 270 of the implementation
 * file "read_only_numerical_quick_median_detail.hpp". There will soon be an article
 * on my personal web site with more details about background, context, existing work, etc.
 */

#include "read_only_numerical_quick_median_detail.hpp"
namespace median_project
{

/*
 * Function read_only_numerical_quick_median
 * ==============================================
 *
 * Non-modifying median algorithm for numerical data in average N * log(N) time.
 * Iterators need to be forward or better.
 *
 */
template <typename Iterator>
std::pair<Iterator, Iterator> read_only_numerical_quick_median(Iterator begin, Iterator end)
{
    no_op_median_performance_stats::no_op_performance_stats performance_stats;
    return read_only_numerical_quick_median_detail::read_only_numerical_quick_median_internal(
        begin, end, performance_stats);
}
} // end namespace median_project

#endif // TMB_READ_ONLY_NUMERICAL_QUICK_MEDIAN_08_24_2015_HPP
