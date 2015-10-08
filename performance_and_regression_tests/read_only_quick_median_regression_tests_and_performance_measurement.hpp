// Regression tests with performance measurenment for OptimizedBruteForceMedian algorithm.
//

#pragma once
#include <cmath>
#include <random>
#include "median_algorithms/read_only_quick_median.hpp"
#include "median_algorithms/read_only_numerical_quick_median.hpp"
#include "median_algorithms/numerical_quick_median.hpp"
using namespace median_project;

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
     * read_only_quick_median
     * read_only_quick_median_random_data
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
      *
      * NOTE: The methods that increment the counts are private. Therefore, the median
      * helper functions that call these methods must be friends of this class. Hence
      * the lengthy friend declarations here. This is of course not very important. It
      * was more an experiment in friends and templates.
      */

    class performance_stats
    {

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

        /*
         * NOTE: Strictly speaking, these methods should be private. But the friends thing
         * was really getting out of hand. Really.
         */

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
    static void test_top_level_algorithms();
    static void test_numerical_median_for_distributions();
    static void test_numerical_median_for_uniform_distribution(size_t num_elems, int num_reps, std::mt19937 &generator);
    static void test_numerical_median_for_normal_distribution(size_t num_elems, int num_reps, std::mt19937 &generator);
    static void
    test_numerical_median_for_exponential_distribution(size_t num_elems, int num_reps, std::mt19937 &generator);
    static void
    test_modifying_numerical_median_for_distributions();

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
