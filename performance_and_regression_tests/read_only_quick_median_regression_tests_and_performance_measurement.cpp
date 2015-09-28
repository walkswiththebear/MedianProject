// Regression tests with performance measurenment for OptimizedBruteForceMedian algorithm.
//

#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <random>
#include "read_only_quick_median_regression_tests_and_performance_measurement.hpp"
#include "performance_and_regression_tests/include/any_iterator/any_iterator.hpp"

void read_only_quick_median_regression_tests_and_performance_measurement::run_tests()
{
    std::cout << "Testing top level algorithms\n";
    std::cout << "============================\n\n";
    test_top_level_algorithms();

    std::cout << "Testing read_only_quick_median for (mostly) random access iterators\n";
    std::cout << "===================================================================\n\n";
    m_which_algorithm = 0;
    test_algorithm();

    std::cout << "\n";
    std::cout << "Testing read_only_quick_median for bidirectional iterators\n";
    std::cout << "==========================================================\n\n";
    m_which_algorithm = 1;
    test_algorithm();

    std::cout << "\n";
    std::cout << "Testing read_only_quick_median for forward iterators\n";
    std::cout << "====================================================\n\n";
    m_which_algorithm = 2;
    test_algorithm();

    std::cout << "\n";
    std::cout << "Testing read_only_quick_median_random_data\n";
    std::cout << "==========================================\n\n";
    m_which_algorithm = 3;
    test_algorithm();

    std::cout << "\n";
    std::cout << "Testing read_only_quick_median_random_data for bidirectional iterators\n";
    std::cout << "======================================================================\n\n";
    m_which_algorithm = 4;
    test_algorithm();

    std::cout << "\n";
    std::cout << "Testing read_only_quick_median_random_data for forward iterators\n";
    std::cout << "================================================================\n\n";
    m_which_algorithm = 5;
    test_algorithm();
    
    std::cout << "\n";
    std::cout << "Testing read_only_numerical_quick_median\n";
    std::cout << "========================================\n\n";
    m_which_algorithm = 6;
    test_algorithm();

    std::cout << m_check_true_count << " tests performed.\n";
}

void read_only_quick_median_regression_tests_and_performance_measurement::test_top_level_algorithms()
{
    std::vector<double> vec = std::vector<double>{ 7., 6., 5., 4., 3., 2., 1. };
    
    std::cout << "Testing read_only_quick_median...";
    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> median_pair = 
        read_only_quick_median(vec.cbegin(), vec.cend());
    check_true(*median_pair.first == 4.0);
    check_true(*median_pair.second == 4.0);
    std::cout << "done.\n";
    
    std::cout << "Testing read_only_quick_median_random_data...";
    median_pair =
        read_only_quick_median_random_data(vec.cbegin(), vec.cend());
    check_true(*median_pair.first == 4.0);
    check_true(*median_pair.second == 4.0);
    std::cout << "done.\n";

    std::cout << "Testing read_only_numerical_quick_median...";
    check_true(read_only_numerical_quick_median(vec.cbegin(), vec.cend()) == 4.0);
    std::cout << "done.\n\n";
}

template <typename Iterator, typename PerformanceStats>
typename std::iterator_traits<Iterator>::value_type
read_only_quick_median_regression_tests_and_performance_measurement::tested_algorithm(
    Iterator begin,
    Iterator end,
    PerformanceStats &performance_stats)
{
    typedef typename std::iterator_traits<Iterator>::value_type value_type;
    typedef typename IteratorTypeErasure::any_iterator<value_type, std::bidirectional_iterator_tag, value_type const &>
    bidirectional_it;
    typedef typename IteratorTypeErasure::any_iterator<value_type, std::forward_iterator_tag, value_type const &>
    forward_it;
    if (m_which_algorithm == 0)
    {
        std::pair<Iterator, Iterator> median_pos_pair =
            read_only_quick_median_detail::read_only_quick_median_internal(
                begin,
                end,
                read_only_quick_median_detail::standard_pivoting_strategy(),
                performance_stats);

        return (*std::get<0>(median_pos_pair) + *std::get<1>(median_pos_pair)) / 2.0;
    }
    else if (m_which_algorithm == 1)
    {
        std::pair<bidirectional_it, bidirectional_it> median_pos_pair =
            read_only_quick_median_detail::read_only_quick_median_internal(
                bidirectional_it(begin),
                bidirectional_it(end),
                read_only_quick_median_detail::standard_pivoting_strategy(),
                performance_stats);

        return (*std::get<0>(median_pos_pair) + *std::get<1>(median_pos_pair)) / 2.0;
    }
    else if (m_which_algorithm == 2)
    {
        std::pair<forward_it, forward_it> median_pos_pair =
            read_only_quick_median_detail::read_only_quick_median_internal(
                forward_it(begin),
                forward_it(end),
                read_only_quick_median_detail::standard_pivoting_strategy(),
                performance_stats);

        return (*std::get<0>(median_pos_pair) + *std::get<1>(median_pos_pair)) / 2.0;
    }
    else if (m_which_algorithm == 3)
    {
        std::pair<Iterator, Iterator> median_pos_pair =
            read_only_quick_median_detail::read_only_quick_median_internal(
                begin,
                end,
                read_only_quick_median_detail::pivoting_strategy_for_random_data(),
                performance_stats);

        return (*std::get<0>(median_pos_pair) + *std::get<1>(median_pos_pair)) / 2.0;
    }
    else if (m_which_algorithm == 4)
    {
        std::pair<bidirectional_it, bidirectional_it> median_pos_pair =
            read_only_quick_median_detail::read_only_quick_median_internal(
                bidirectional_it(begin),
                bidirectional_it(end),
                read_only_quick_median_detail::pivoting_strategy_for_random_data(),
                performance_stats);

        return (*std::get<0>(median_pos_pair) + *std::get<1>(median_pos_pair)) / 2.0;
    }
    else if (m_which_algorithm == 5)
    {
        std::pair<forward_it, forward_it> median_pos_pair =
            read_only_quick_median_detail::read_only_quick_median_internal(
                forward_it(begin),
                forward_it(end),
                read_only_quick_median_detail::pivoting_strategy_for_random_data(),
                performance_stats);

        return (*std::get<0>(median_pos_pair) + *std::get<1>(median_pos_pair)) / 2.0;
    }
    else if (m_which_algorithm == 6)
    {
        read_only_numerical_quick_median_detail::uniform_distribution_pivot pivot_calculator;
        return read_only_numerical_quick_median_detail::read_only_numerical_quick_median_internal(
            begin, end, pivot_calculator, performance_stats);
    }
    else
    {
        check_true(false);
        return static_cast<typename std::iterator_traits<Iterator>::value_type>(0);
    }
}

void read_only_quick_median_regression_tests_and_performance_measurement::test_algorithm()
{
    std::cout << "Running median white box tests...";
    run_white_box_tests();
    std::cout << "done.\n";
    std::cout << "Running median black box tests...\n\n";
    run_random_tests();
    std::cout << "done.\n";
}

void read_only_quick_median_regression_tests_and_performance_measurement::run_white_box_tests()
{
    no_op_median_performance_stats stats;

    // One element.
    std::list<double> lis{ 1. };
    verify_median(lis.begin(), lis.end(), stats);
    // Be redundant to test verify_median.
    check_true(tested_algorithm(lis.begin(), lis.end(), stats) == 1.0);

    // One element an uneven number of times
    lis = std::list<double>{ 1., 1., 1., 1., 1., 1., 1. };
    verify_median(lis.begin(), lis.end(), stats);
    // Be redundant to test verify_median.
    check_true(tested_algorithm(lis.begin(), lis.end(), stats) == 1.0);

    // One element an even number of times
    lis.push_back(1.0);
    verify_median(lis.begin(), lis.end(), stats);
    // Be redundant to test verify_median.
    check_true(tested_algorithm(lis.begin(), lis.end(), stats) == 1.0);

    // Leave it in there...
    std::vector<double> vec = std::vector<double>{ 2., 3., 6., 7., 5., 4., 1 };
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Descending sequence of uneven length
    vec = std::vector<double>{ 7., 6., 5., 4., 3., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    std::vector<double> tmp_vec(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Ascending sequence of uneven length
    std::reverse(vec.begin(), vec.end());
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Ascending sequence of even length
    vec.push_back(8.);
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.5);

    // Descending sequence of even length
    std::reverse(vec.begin(), vec.end());
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.5);

    // Uneven number of elements, "almost descending"
    vec.push_back(9.);
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 5.0);

    // Uneven number of elements, "almost ascending"
    std::reverse(vec.begin(), vec.end());
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 5.0);

    // Add a duplicate at the beginning, for the heck of it.
    vec.insert(vec.begin(), 9.0);
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 5.5);

    // Add another duplicate at the beginning, for the heck of it.
    vec.insert(vec.begin(), 9.0);
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 6.);

    // Add yet another duplicate at the beginning, for the heck of it.
    vec.insert(vec.begin(), 9.0);
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 6.5);

    // Even number of elements, "almost ascending"
    vec = std::vector<double>{ 8., 1., 2., 3., 4., 5., 6., 7. };
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.5);

    // Even number of elements, "almost descending"
    std::reverse(vec.begin(), vec.end());
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.5);

    // Add a duplicate at the end, for the heck of it.
    vec.push_back(8.);
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 5.0);

    // Add another duplicate at the end, for the heck of it.
    vec.push_back(8.);
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 5.5);

    // Add yet another duplicate at the end, for the heck of it.
    vec.push_back(8.);
    verify_median(vec.begin(), vec.end(), stats);
    tmp_vec = std::vector<double>(vec.begin(), vec.end());
    run_a_few_shuffles(tmp_vec.begin(), tmp_vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 6.0);

    // Even number of elements, median between two duplicates.
    vec = std::vector<double>{ 7., 6., 5., 4., 4., 3., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Uneven number of elements, median inside a set of duplicates
    vec = std::vector<double>{ 7., 6., 5., 4., 4., 4., 3., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Even number of elements, median in the middle of four duplicates.
    vec = std::vector<double>{ 7., 6., 5., 4., 4., 4., 4., 3., 2., 1. };
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());

    // Leave it in there...
    vec = std::vector<double>{ 3., 4., 4., 1., 5., 6., 4., 7., 3., 2., 4. };
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Uneven number of elements, median inside a set of duplicates (again).
    vec = std::vector<double>{ 7., 6., 5., 4., 4., 4., 4., 3., 3., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Uneven number of elements, median inside a set of duplicates (again,
    // with more duplicates).
    vec = std::vector<double>{ 7., 6., 5., 5., 4., 4., 4., 3., 3., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Even number of elements, median between third and fourth of four duplicates.
    vec = std::vector<double>{ 7., 6., 5., 4., 4., 4., 4., 3., 3., 3., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Uneven number of elements, median first of four duplicates, preceded by a non-duplicate.
    vec = std::vector<double>{ 9., 8., 7., 6., 5., 4., 4., 4., 4., 3., 3. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Even number of elements, median between a non-duplicate and a set of duplicates.
    vec = std::vector<double>{ 9., 8., 7., 6., 5., 4., 4., 4., 4., 3. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.5);

    // Uneven number of elements, median first of four duplicates, preceded by four duplicates.
    vec = std::vector<double>{ 5., 4., 4., 4., 4., 3., 3., 3., 3., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 3.0);

    // Uneven number of elements, median last of four duplicates, followed by a non-duplicate.
    vec = std::vector<double>{ 7., 6., 5., 4., 4., 4., 4., 3., 2., 2., 2., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Even number of elements, median between a set of duplicates and a non-duplicate.
    vec = std::vector<double>{ 7., 6., 5., 4., 4., 4., 4., 3., 2., 2., 2., 2., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 3.5);

    // Uneven number of elements, median last of four duplicates, followed by four duplicates.
    vec = std::vector<double>{ 7., 6., 5., 4., 4., 4., 4., 3., 3., 3., 3., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 4.0);

    // Even number of elements, median between two sets of duplicates.
    vec = std::vector<double>{ 7., 6., 5., 4., 4., 4., 4., 3., 3., 3., 3., 3., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 3.5);

    // Uneven number of elements, median a non-duplicate followed by four duplicates
    vec = std::vector<double>{ 9., 8., 7., 6., 5., 4., 4., 4., 4. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 5.0);

    // Uneven number of elements, median a non-duplicate preceded by a set of duplicates.
    vec = std::vector<double>{ 6., 5., 4., 4., 4., 4., 3., 2., 2., 2., 2., 2., 1. };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == 3.0);

    // This one didn't work with the first version of the numerical median algorithm.
    vec = std::vector<double>{ -4., -5., -4., -3., -1., 2., 3., 4., 4., 3. };
    verify_median(vec.begin(), vec.end(), stats);
    // Be redundant to test verify_median.
    check_true(tested_algorithm(vec.begin(), vec.end(), stats) == .5);

    // Just for the heck of it...
    vec.push_back(0.);
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
    //
    vec.push_back(0.);
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());

    // To test a case that is kind of special in torben.c
    vec = std::vector<double>{ .5, 2., 3., 3.5 };
    verify_median(vec.begin(), vec.end(), stats);
    run_a_few_shuffles(vec.begin(), vec.end());
}

template <typename Iterator>
void read_only_quick_median_regression_tests_and_performance_measurement::run_a_few_shuffles(Iterator begin,
                                                                                             Iterator end)
{
    no_op_median_performance_stats stats;

    for (int i = 0; i < 5; ++i)
    {
        std::shuffle(begin, end, std::default_random_engine());
        verify_median(begin, end, stats);
    }
}

void read_only_quick_median_regression_tests_and_performance_measurement::run_random_tests()
{

    int fixed_length = 1;
    for (int i = 1; i <= m_log10_of_size_of_largest_data_set; ++i)
    {
        fixed_length *= 10;
        run_fixed_length_test(fixed_length);
        run_fixed_length_test(fixed_length + 1);
    }
}

void read_only_quick_median_regression_tests_and_performance_measurement::run_fixed_length_test(int len)
{

    run_no_duplicates_fixed_length_test(len);
    run_fixed_length_test_with_duplicates(len, -2 * len, 2 * len);
    run_fixed_length_test_with_duplicates(len, -len / 2, len / 2);
    run_fixed_length_test_with_duplicates(len, -len / 5, len / 5);
    run_fixed_length_test_with_duplicates(len, 0, len / 10);
}

void read_only_quick_median_regression_tests_and_performance_measurement::run_no_duplicates_fixed_length_test(int len)
{

    // Fill vector with integers 0 through len - 1
    std::vector<double> sequence(len);
    double count = 0.0;
    for (double &elem : sequence)
    {
        elem = count;
        count += 1.0;
    }
    check_true(sequence[len - 1] == static_cast<double>(len) - 1);

    // Do it once for the sorted sequence, but let it not distort the
    // average complexity results.
    no_op_median_performance_stats unused_stats;
    verify_median(sequence.cbegin(), sequence.cend(), unused_stats);
    run_a_few_shuffles(sequence.begin(), sequence.end());

    // Fill vector with a "skewed uniform" distribution, shuffle,
    // check median, and keep track of performance stats
    //
    std::cout << "  " << len << " elements, no duplicates, " << m_monte_carlo_count << "  times.\n" << std::flush;
    performance_stats stats = performance_stats();
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_int_distribution<> uniform_distribution(1, 1000);
    for (int i = 0; i < m_monte_carlo_count; ++i)
    {
        std::vector<double> sequence(len);
        double previous_elem = 0.0;
        double count = 0.0;
        for (double &elem : sequence)
        {
            elem = previous_elem + count + static_cast<double>(uniform_distribution(generator));
            previous_elem = elem;
            ++count;
        }

        std::random_shuffle(sequence.begin(), sequence.end());
        verify_median(sequence.cbegin(), sequence.cend(), stats);
    }

    std::cout << "    n=" << len << ", (average pivot count)/log(n)=" << stats.get_pivot_count_average_to_logN_ratio()
              << ", (average num comparisons)/nlog(n)=" << stats.get_comparison_count_average_to_NlogN_ratio() << "\n\n"
              << std::flush;
}

void read_only_quick_median_regression_tests_and_performance_measurement::run_fixed_length_test_with_duplicates(
    int len,
    int lowerBound,
    int upperBound)
{

    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_int_distribution<> uniform_distribution(lowerBound, upperBound);

    std::vector<double> sequence(len);
    performance_stats stats = performance_stats();
    std::cout << "  " << len << " elements between " << lowerBound << " and " << upperBound << ", "
              << m_monte_carlo_count << "  times.\n" << std::flush;
    for (int i = 0; i < m_monte_carlo_count; ++i)
    {

        // Fill vector with random integers
        for (double &elem : sequence)
        {
            elem = uniform_distribution(generator);
        }

        verify_median(sequence.cbegin(), sequence.cend(), stats);
    }

    std::cout << "    n=" << len << ", (average pivot count)/log(n)=" << stats.get_pivot_count_average_to_logN_ratio()
              << ", (average num comparisons)/nlog(n)=" << stats.get_comparison_count_average_to_NlogN_ratio() << "\n\n"
              << std::flush;
}

template <typename Iterator, typename PerformanceStats>
void
read_only_quick_median_regression_tests_and_performance_measurement::verify_median(Iterator begin,
                                                                                   Iterator end,
                                                                                   PerformanceStats &performance_stats)
{
    std::vector<double> sequence_copy(begin, end);
    int len = sequence_copy.size();

    // Calculate the median with our algorithm.
    double my_median = tested_algorithm(begin, end, performance_stats);

    // Use std::nth_element to check our result.
    //
    double their_median = 0.0;
    std::vector<double>::iterator upper_mid_pos = sequence_copy.begin() + len / 2;
    std::nth_element(sequence_copy.begin(), upper_mid_pos, sequence_copy.end());
    double their_median_upper = *upper_mid_pos;
    if (len % 2 == 0)
    {
        std::vector<double>::iterator lower_mid_pos = upper_mid_pos - 1;
        std::nth_element(sequence_copy.begin(), lower_mid_pos, sequence_copy.end());
        their_median = (their_median_upper + *lower_mid_pos) / 2.0;
    }
    else
    {
        their_median = their_median_upper;
    }

    check_true(my_median == their_median);
}

void read_only_quick_median_regression_tests_and_performance_measurement::check_true(bool cond)
{
    ++m_check_true_count;
    if (!cond)
    {
        throw std::runtime_error("Test failed");
    }
}

int read_only_quick_median_regression_tests_and_performance_measurement::m_check_true_count = 0;
int read_only_quick_median_regression_tests_and_performance_measurement::m_monte_carlo_count =
    10000; // recommended 10000
int read_only_quick_median_regression_tests_and_performance_measurement::m_log10_of_size_of_largest_data_set =
    4; // recommended 5
int read_only_quick_median_regression_tests_and_performance_measurement::m_which_algorithm = -1;
