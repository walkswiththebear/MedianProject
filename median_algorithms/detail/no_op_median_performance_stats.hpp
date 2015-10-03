//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_NO_OP_MEDIAN_PERFORMANCE_STATS_08_24_2015_HPP
#define TMB_NO_OP_MEDIAN_PERFORMANCE_STATS_08_24_2015_HPP

/*
 * No-op performance stats for median algorithms in production.
 */

#include <tuple>
#include <assert.h>

/*
 * NOTE: The methods that increment the counts in the performance stats class are
 * private. Therefore, the median helper functions that call these methods must
 * be friends of the performance stats class. Hence the lengthy forward and friend
 * declarations here. This is of course not very important. It was more an
 * experiment in friends and templates.
 */

namespace median_project
{

namespace read_only_quick_median_detail
{
template <typename Iterator, typename PivotingsStrategy, typename PerformanceStats>
typename std::pair<Iterator, Iterator> read_only_quick_median_internal(Iterator begin,
                                                                       Iterator end,
                                                                       PivotingsStrategy pivoting_strategy,
                                                                       PerformanceStats &performance_stats);

template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int>
trim_sequence_left(Iterator active_sequence_begin,
                   bool median_lower_bound_found,
                   typename std::iterator_traits<Iterator>::value_type median_lower_bound,
                   bool median_upper_bound_found,
                   typename std::iterator_traits<Iterator>::value_type median_upper_bound,
                   PerformanceStats &performance_stats);

template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int>
trim_sequence_right(Iterator active_sequence_end,
                    bool median_lower_bound_found,
                    typename std::iterator_traits<Iterator>::value_type median_lower_bound,
                    bool median_upper_bound_found,
                    typename std::iterator_traits<Iterator>::value_type median_upper_bound,
                    PerformanceStats &performance_stats,
                    std::forward_iterator_tag);

template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int>
trim_sequence_right(Iterator active_sequence_end,
                    bool median_lower_bound_found,
                    typename std::iterator_traits<Iterator>::value_type median_lower_bound,
                    bool median_upper_bound_found,
                    typename std::iterator_traits<Iterator>::value_type median_upper_bound,
                    PerformanceStats &performance_stats,
                    std::bidirectional_iterator_tag);

class standard_pivoting_strategy;

class pivoting_strategy_for_random_data;

template <typename Iterator, typename PerformanceStats>
std::tuple<int, int, int> count_elements(Iterator begin,
                                         Iterator end,
                                         typename std::iterator_traits<Iterator>::value_type pivot,
                                         PerformanceStats &performance_stats);
}

namespace read_only_numerical_quick_median_detail
{
template <typename Iterator, typename PivotCalculator, typename PerformanceStats>
double read_only_numerical_quick_median_internal(Iterator begin,
                                                 Iterator end,
                                                 PivotCalculator &pivot_calculator,
                                                 PerformanceStats &performance_stats);

template <typename Iterator, typename PerformanceStats>
std::tuple<int, double, double, double>
get_initial_sequence_data(Iterator begin, Iterator end, PerformanceStats &performance_stats);

template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int> trim_sequence_left(Iterator active_sequence_begin,
                                                  double median_lower_bound,
                                                  double median_upper_bound,
                                                  PerformanceStats &performance_stats);

template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int> trim_sequence_right(Iterator active_sequence_end,
                                                   double median_lower_bound,
                                                   double median_upper_bound,
                                                   PerformanceStats &performance_stats,
                                                   std::forward_iterator_tag);

template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int> trim_sequence_right(Iterator active_sequence_end,
                                                   double median_lower_bound,
                                                   double median_upper_bound,
                                                   PerformanceStats &performance_stats,
                                                   std::bidirectional_iterator_tag);

class pivot_functor_base;

template <typename Iterator, typename PerformanceStats>
std::tuple<int, int, int, double, double>
count_elements(Iterator begin, Iterator end, double pivot, PerformanceStats &performance_stats);
}

/**
* No-op performance stats class for median, to be used as the default in the algorithm.
*/
class no_op_median_performance_stats
{
    // Friends
    // =======

    template <typename Iterator, typename PivotingsStrategy, typename PerformanceStats>
    friend typename std::pair<Iterator, Iterator>
    read_only_quick_median_detail::read_only_quick_median_internal(Iterator begin,
                                                                   Iterator end,
                                                                   PivotingsStrategy pivoting_strategy,
                                                                   PerformanceStats &performance_stats);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<Iterator, int, int> read_only_quick_median_detail::trim_sequence_left(
        Iterator active_sequence_begin,
        bool median_lower_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_lower_bound,
        bool median_upper_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_upper_bound,
        PerformanceStats &performance_stats);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<Iterator, int, int> read_only_quick_median_detail::trim_sequence_right(
        Iterator active_sequence_end,
        bool median_lower_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_lower_bound,
        bool median_upper_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_upper_bound,
        PerformanceStats &performance_stats,
        std::forward_iterator_tag);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<Iterator, int, int> read_only_quick_median_detail::trim_sequence_right(
        Iterator active_sequence_end,
        bool median_lower_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_lower_bound,
        bool median_upper_bound_found,
        typename std::iterator_traits<Iterator>::value_type median_upper_bound,
        PerformanceStats &performance_stats,
        std::bidirectional_iterator_tag);

    friend class read_only_quick_median_detail::standard_pivoting_strategy;

    friend class read_only_quick_median_detail::pivoting_strategy_for_random_data;

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<int, int, int>
    read_only_quick_median_detail::count_elements(Iterator begin,
                                                  Iterator end,
                                                  typename std::iterator_traits<Iterator>::value_type pivot,
                                                  PerformanceStats &performance_stats);

    template <typename Iterator, typename PivotCalculator, typename PerformanceStats>
    friend double read_only_numerical_quick_median_detail::read_only_numerical_quick_median_internal(
        Iterator begin,
        Iterator end,
        PivotCalculator &pivot_calculator,
        PerformanceStats &performance_stats);

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<int, double, double, double>
    read_only_numerical_quick_median_detail::get_initial_sequence_data(Iterator begin,
                                                                       Iterator end,
                                                                       PerformanceStats &performance_stats);

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

    friend class read_only_numerical_quick_median_detail::pivot_functor_base;

    template <typename Iterator, typename PerformanceStats>
    friend std::tuple<int, int, int, double, double>
    read_only_numerical_quick_median_detail::count_elements(Iterator begin,
                                                            Iterator end,
                                                            double pivot,
                                                            PerformanceStats &performance_stats);

    // No-op methods
    // =============

  public:
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
} // end namespace median_project

#endif // TMB_NO_OP_MEDIAN_PERFORMANCE_STATS_08_24_2015_HPP
