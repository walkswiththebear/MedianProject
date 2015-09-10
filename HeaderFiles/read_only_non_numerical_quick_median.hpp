//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_READ_ONLY_NON_NUMERICAL_QUICK_MEDIAN_07_21_2015_HPP
#define TMB_READ_ONLY_NON_NUMERICAL_QUICK_MEDIAN_07_21_2015_HPP

/*
 * Non-modifying median algorithm in average N * log(N) time. This is appropriate
 * for situations where the data must not be modified, and there are not enough
 * resources to make a copy of the data.
 *
 * The algorithm works for numerical and non-numerical data, that is, elements
 * do not have to support arithmetic operations such as calculating a mean (midpoint)
 * of two elements. For numerical data, there is a better algorithm, namely,
 * read_only_numerical_quick_median.
 *
 * Arguments are begin and end iterators to the data. These do not have to be better
 * than forward iterators. The return value is a pair of iterators. For an uneven
 * number of elements, both iterators point to the median. For an even number of
 * elements, they point to the left and right endpoint of the median interval.
 *
 * There are two versions of this algorithm. The first is for all kinds of data,
 * including sorted data. Performance improves if the iterators are bidirectional
 * rather than just forward, and even more if they are random access.
 *
 * The second version is for the special situation where the data is essentially the
 * result of a random walk, that is, there is no reason to assume that sorted data
 * is particularly likely to occur. In this case, the second version tends to perform
 * noticeably better than the first. This is true in particular for iterators that are
 * not random access. The performance of the second version does not vary significantly
 * for different iterator categories.
 *
 * A brief description of the algorithm can be found at lines 372 - 489 of the implementation
 * file "read_only_non_numerical_quick_median_detail.hpp".
 */

#include "read_only_non_numerical_quick_median_detail.hpp"
#include "no_op_median_performance_stats.hpp"

namespace median_project
{

/**
 * Function read_only_non_numerical_quick_median
 * ==============================================
 *
 * Non-modifying median algorithm for non-numerical and numerical data in average
 * N * log(N) time. For numerical data, consider using read_only_numerical_quick_median.
 *
 * For best results, use random access iterators. Bidirectional iterators will work
 * well, but performance may suffer. Forward iterators will cause noticeable
 * performance degradation. In the special case where the data is essentially the
 * result of a random walk (sorted data is not particularly likely), the performance
 * hit caused by non-random-access iterators can be avoided by using the alternate
 * version read_only_non_numerical_quick_median_random_data.
 */
template <typename Iterator>
std::pair<Iterator, Iterator> read_only_non_numerical_quick_median(Iterator begin, Iterator end)
{
    no_op_median_performance_stats performance_stats;
    return read_only_non_numerical_quick_median_detail::read_only_non_numerical_quick_median_internal(
        begin, end, read_only_non_numerical_quick_median_detail::standard_pivoting_strategy(), performance_stats);
}

/**
 * Function read_only_non_numerical_quick_median_random_data
 * =========================================================
 *
 * Non-modifying median algorithm for non-numerical and numerical data in average 
 * N * log(N) time. This version of the algorithm improves performance if the data is 
 * essentially the result of a random walk, that is, it has no particular likelihood 
 * of being sorted. This is true in particular when the iterators used are no better 
 * than forward iterators. For numerical data, consider using read_only_numerical_quick_median.
 *
 * IMPORTANT: For sorted data, the performance of this version of the algorithm will degrade
 * to N^2. Use this algorithm only if your data is essentially random, that is, has little or
 * no tendency to be sorted. If sorted data is likely, use read_only_non_numerical_quick_median
 * instead, regardless of the iterator category.
 */
template <typename Iterator>
std::pair<Iterator, Iterator> read_only_non_numerical_quick_median_random_data(Iterator begin, Iterator end)
{
    no_op_median_performance_stats performance_stats;
    return read_only_non_numerical_quick_median_detail::read_only_non_numerical_quick_median_internal(
        begin,
        end,
        read_only_non_numerical_quick_median_detail::pivoting_strategy_for_random_data(),
        performance_stats);
}
} // end namespace median_project

#endif // TMB_READ_ONLY_NON_NUMERICAL_QUICK_MEDIAN_07_21_2015_HPP
