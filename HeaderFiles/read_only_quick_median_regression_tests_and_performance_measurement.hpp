// Regression tests with performance measurenment for OptimizedBruteForceMedian algorithm.
//

#pragma once
#include <cmath>
#include "read_only_non_numerical_quick_median.hpp"
#include "read_only_numerical_quick_median.hpp"
using namespace tmb_algorithms;

/*
 * Regression test class for read-only median algorithms. The only public
 * method is run_tests(). See its documentation below for details.
 *
 * TO DO: Make this more elegant and understandable, possibly by using a
 * general test framework.
 */
class read_only_quick_median_regression_tests_and_performance_measurement
{

  public:
    /*
     * Tests the following algorithms:
     *
     * read_only_non_numerical_quick_median
     * read_only_non_numerical_quick_median_random_data
     * read_only_numerical_quick_median

     * For each of these algorithms, the following tests are run:
     *
     * 1) A few white box tests. These calculates the median of hand-made
     *    data that targets certain "critical" cases of the implementation.
     *
     * 2) A configurable number of black box tests. These calculate the median
     *    of data that is generated with random number generator, both with and
     *    without duplicates. The variables m_monte_carlo_count and
     *    m_log10_of_size_of_largest_data_set configure the test. Please search
     *    the code for these variables for details.
     */
    static void run_tests();

  private:
    /**
      * Class for keeping track of performance
      */
    class performance_stats
    {
        template <typename Iterator, typename PivotingsStrategy, typename PerformanceStats>
        friend typename std::pair<Iterator, Iterator>
        read_only_non_numerical_quick_median_detail::read_only_non_numerical_quick_median_internal(
            Iterator begin,
            Iterator end,
            PivotingsStrategy pivoting_strategy,
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
        friend std::tuple<Iterator, int, int> read_only_non_numerical_quick_median_detail::trim_sequence_left(
            Iterator active_sequence_begin,
            bool median_lower_bound_found,
            typename std::iterator_traits<Iterator>::value_type median_lower_bound,
            bool median_upper_bound_found,
            typename std::iterator_traits<Iterator>::value_type median_upper_bound,
            PerformanceStats &performance_stats);

        template <typename Iterator, typename PerformanceStats>
        friend std::tuple<Iterator, int, int>
        read_only_numerical_quick_median_detail::trim_sequence_left(Iterator active_sequence_begin,
                                                                    double median_lower_bound,
                                                                    double median_upper_bound,
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
        friend std::tuple<Iterator, int, int>
        read_only_numerical_quick_median_detail::trim_sequence_right(Iterator active_sequence_end,
                                                                     double median_lower_bound,
                                                                     double median_upper_bound,
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

        template <typename Iterator, typename PerformanceStats>
        friend std::tuple<Iterator, int, int>
        read_only_numerical_quick_median_detail::trim_sequence_right(Iterator active_sequence_end,
                                                                     double median_lower_bound,
                                                                     double median_upper_bound,
                                                                     PerformanceStats &performance_stats,
                                                                     std::bidirectional_iterator_tag);

        friend class read_only_non_numerical_quick_median_detail::standard_pivoting_strategy;

        friend class read_only_non_numerical_quick_median_detail::pivoting_strategy_for_random_data;

        template <typename Iterator, typename PerformanceStats>
        friend std::tuple<int, int, int> read_only_non_numerical_quick_median_detail::count_elements(
            Iterator begin,
            Iterator end,
            typename std::iterator_traits<Iterator>::value_type compare_to_element,
            PerformanceStats &performance_stats);

        template <typename Iterator, typename PerformanceStats>
        friend std::tuple<int, int, int, double, double>
        read_only_numerical_quick_median_detail::count_elements(Iterator begin,
                                                                Iterator end,
                                                                double compare_to_element,
                                                                PerformanceStats &performance_stats);

      public:
        performance_stats()
        {
        }
        double get_pivot_count_average_to_logN_ratio()
        {
            return m_pivot_count_average / std::log2(static_cast<double>(m_sequence_length));
        }
        double get_comparison_count_average_to_NlogN_ratio()
        {
            return m_comparison_count_average /
                   (static_cast<double>(m_sequence_length) * std::log2(static_cast<double>(m_sequence_length)));
        }

      private:
        void set_sequence_length(int len)
        {
            m_sequence_length = len;
        }
        void increment_pivot_count()
        {
            ++m_pivot_count;
        }
        void add_comparisons(int comps)
        {
            m_comparison_count += comps;
        }
        void update_averages()
        {

            m_pivot_count_average =
                m_pivot_count_average * (static_cast<double>(m_set_averages_call_count) /
                                         (static_cast<double>(m_set_averages_call_count) + 1.0)) +
                static_cast<double>(m_pivot_count) * (1.0 / (static_cast<double>(m_set_averages_call_count) + 1.0));

            m_comparison_count_average =
                m_comparison_count_average * (static_cast<double>(m_set_averages_call_count) /
                                              (static_cast<double>(m_set_averages_call_count) + 1.0)) +
                static_cast<double>(m_comparison_count) *
                    (1.0 / (static_cast<double>(m_set_averages_call_count) + 1.0));

            m_pivot_count = 0;
            m_comparison_count = 0;
            ++m_set_averages_call_count;
        }

        int m_sequence_length = 0;
        int m_pivot_count = 0;
        int m_comparison_count = 0;
        double m_pivot_count_average = 0;
        double m_comparison_count_average = 0;
        int m_set_averages_call_count = 0;
    };

    /*
     * This function is called by the tests in lieu of the algorithm to be tested.
     * The function branches to the currently to be tested algorithm by testing the
     * member variable m_which_algorithm.
     */
    template <typename Iterator, typename PerformanceStats>
    static typename std::iterator_traits<Iterator>::value_type
    tested_algorithm(Iterator begin, Iterator end, PerformanceStats &performance_stats);

    // Internal helper functions
    //
    static void test_algorithm();
    static void run_white_box_tests();
    template <typename Iterator> static void run_a_few_shuffles(Iterator begin, Iterator end);
    static void run_random_tests();
    static void run_fixed_length_test(int len);
    static void run_no_duplicates_fixed_length_test(int len);
    static void run_fixed_length_test_with_duplicates(int len, int lower_bound, int upper_bound);
    template <typename Iterator, typename PerformanceStats>
    static void verify_median(Iterator begin, Iterator end, PerformanceStats &performance_stats);
    static void check_true(bool cond);

    // Static member variables
    static int m_check_true_count;
    static int m_monte_carlo_count;
    static int m_log10_of_size_of_largest_data_set;
    static int m_which_algorithm;
};
