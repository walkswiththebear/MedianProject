//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

// Revision History
// ================
//
// 24 Aug 2015 (Thomas Becker) Created

#ifndef TMB_NO_OP_MEDIAN_PERFORMANCE_STATS_08_24_2015_HPP
#define TMB_NO_OP_MEDIAN_PERFORMANCE_STATS_08_24_2015_HPP

/*
 * No-op performance stats for median algorithms in production.
 */

#include <tuple>
#include <assert.h>
#include "read_only_non_numerical_quick_median_detail.hpp"
#include "read_only_numerical_quick_median_detail.hpp"

namespace tmb_algorithms
{
namespace no_op_median_performance_stats
{

/*
* No-op performance stats class, to be used as the default in the algorithm.
*/
class no_op_performance_stats
{
    template <typename Iterator, typename PivotingsStrategy, typename PerformanceStats>
    friend typename std::pair<Iterator, Iterator>
    read_only_non_numerical_quick_median_detail::read_only_non_numerical_quick_median_internal(
        Iterator begin,
        Iterator end,
        PivotingsStrategy pivoting_strategy,
        PerformanceStats &performance_stats);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<Iterator, int, int> read_only_non_numerical_quick_median_detail::trim_sequence_left(
        Iterator active_sequence_begin,
        bool median_lower_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_lower_bound,
        bool median_upper_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_upper_bound,
        PerformanceStats &performance_stats);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<Iterator, int, int> read_only_non_numerical_quick_median_detail::trim_sequence_right(
        Iterator active_sequence_end,
        bool median_lower_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_lower_bound,
        bool median_upper_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_upper_bound,
        PerformanceStats &performance_stats,
        std::forward_iterator_tag);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<Iterator, int, int> read_only_non_numerical_quick_median_detail::trim_sequence_right(
        Iterator active_sequence_end,
        bool median_lower_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_lower_bound,
        bool median_upper_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_upper_bound,
        PerformanceStats &performance_stats,
        std::bidirectional_iterator_tag);

    friend class read_only_non_numerical_quick_median_detail::standard_pivoting_strategy;

    friend class read_only_non_numerical_quick_median_detail::pivoting_strategy_for_random_data;

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<int, int, int> read_only_non_numerical_quick_median_detail::count_elements(
        Iterator begin,
        Iterator end,
        typename std::iterator_traits<Iterator>::value_type pivot,
        PerformanceStats &performance_stats);

    template <typename Iterator, typename PerformanceStats>
    friend double read_only_numerical_quick_median_detail::read_only_numerical_quick_median_internal(
        Iterator begin,
        Iterator end,
        PerformanceStats &performance_stats);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<double, double, int>
    read_only_numerical_quick_median_detail::get_initial_sequence_data(Iterator begin,
                                                                       Iterator end,
                                                                       PerformanceStats& performance_stats);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<Iterator, int, int>
    read_only_numerical_quick_median_detail::trim_sequence_left(Iterator active_sequence_begin,
                                                                double median_lower_bound,
                                                                double median_upper_bound,
                                                                PerformanceStats &performance_stats);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<Iterator, int, int>
    read_only_numerical_quick_median_detail::trim_sequence_right(Iterator active_sequence_end,
                                                                 double median_lower_bound,
                                                                 double median_upper_bound,
                                                                 PerformanceStats &performance_stats,
                                                                 std::forward_iterator_tag);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<Iterator, int, int>
    read_only_numerical_quick_median_detail::trim_sequence_right(Iterator active_sequence_end,
                                                                 double median_lower_bound,
                                                                 double median_upper_bound,
                                                                 PerformanceStats &performance_stats,
                                                                 std::bidirectional_iterator_tag);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<int, int, int, double, double>
    read_only_numerical_quick_median_detail::count_elements(Iterator begin,
                                                            Iterator end,
                                                            double pivot,
                                                            PerformanceStats &performance_stats);

  private:
    void set_sequence_length(int len)
    {
    }
    void increment_pivot_count()
    {
    }
    void add_comparisons(int comps)
    {
    }
    void update_averages()
    {
    }
};
} // end namespace no_op_median_performance_stats
} // end namespace tmb_algorithms

#endif // TMB_NO_OP_MEDIAN_PERFORMANCE_STATS_08_24_2015_HPP
