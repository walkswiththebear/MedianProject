//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_READ_ONLY_NUMERICAL_QUICK_MEDIAN_DETAIL_08_24_2015_HPP
#define TMB_READ_ONLY_NUMERICAL_QUICK_MEDIAN_DETAIL_08_24_2015_HPP

/*
 * Implementation details for read-only numerical quick median algorithm.
 */

#include <algorithm>
#include <tuple>
#include <limits>
#include <assert.h>

namespace median_project
{
namespace read_only_numerical_quick_median_detail
{

/**
 * Calculates the mean of two doubles. Guards against overflow. If the calculated mean
 * lies outside of the interval (very small numbers!), it is set to the min or max of
 * the two numbers, as the case may be.
 */
inline double mean(double x, double y)
{
    double mean = x / 2.0 + y / 2.0;
    mean = std::min(std::max(x, y), std::max(std::min(x, y), mean));
    return mean;
}

/**
 * Trims the sequence at beginning.This is a minor optimization that is just too obvious to
 * be passed up.
 */
template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int> trim_sequence_left(Iterator active_sequence_begin,
                                                  double median_lower_bound,
                                                  double median_upper_bound,
                                                  PerformanceStats &performance_stats)
{
    int num_discarded_elements_less_than_median_lower_bound = 0;
    int num_discarded_elements_greater_than_median_upper_bound = 0;
    while (true)
    {
        if (static_cast<double>(*active_sequence_begin) < median_lower_bound)
        {
            performance_stats.add_comparisons(1);
            ++num_discarded_elements_less_than_median_lower_bound;
            ++active_sequence_begin;
        }
        else if (static_cast<double>(*active_sequence_begin) > median_upper_bound)
        {
            performance_stats.add_comparisons(2);
            ++num_discarded_elements_greater_than_median_upper_bound;
            ++active_sequence_begin;
        }
        else
        {
            performance_stats.add_comparisons(2);
            break;
        }
    }

    return std::tuple<Iterator, int, int>(active_sequence_begin,
                                          num_discarded_elements_less_than_median_lower_bound,
                                          num_discarded_elements_greater_than_median_upper_bound);
}

/**
 * Trims sequence at end. This is a minor optimization that can be applied
 * only with bidirectional iterators or better.
 */
template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int> trim_sequence_right(Iterator active_sequence_end,
                                                   double median_lower_bound,
                                                   double median_upper_bound,
                                                   PerformanceStats &performance_stats,
                                                   std::forward_iterator_tag)
{
    return std::tuple<Iterator, int, int>(active_sequence_end, 0, 0);
}
//
template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int> trim_sequence_right(Iterator active_sequence_end,
                                                   double median_lower_bound,
                                                   double median_upper_bound,
                                                   PerformanceStats &performance_stats,
                                                   std::bidirectional_iterator_tag)
{
    Iterator active_sequence_last = active_sequence_end;
    --active_sequence_last;

    int num_discarded_elements_less_than_median_lower_bound = 0;
    int num_discarded_elements_greater_than_median_upper_bound = 0;
    while (true)
    {
        if (static_cast<double>(*active_sequence_last) < median_lower_bound)
        {
            performance_stats.add_comparisons(1);
            ++num_discarded_elements_less_than_median_lower_bound;
            --active_sequence_last;
        }
        else if (static_cast<double>(*active_sequence_last) > median_upper_bound)
        {
            performance_stats.add_comparisons(2);
            ++num_discarded_elements_greater_than_median_upper_bound;
            --active_sequence_last;
        }
        else
        {
            performance_stats.add_comparisons(2);
            break;
        }
    }

    ++active_sequence_last;
    return std::tuple<Iterator, int, int>(active_sequence_last,
                                          num_discarded_elements_less_than_median_lower_bound,
                                          num_discarded_elements_greater_than_median_upper_bound);
}

/**
* Counts the number of elements less than, greater than, and equal to a given element in the sequence
* [begin, end). Also keeps track of the maximum of elements less than the pivot and the minimum of
* elements greater than the pivot.
*/
template <typename Iterator, typename PerformanceStats>
std::tuple<int, int, int, double, double>
count_elements(Iterator begin, Iterator end, double pivot, PerformanceStats &performance_stats)
{
    int less_than_count = 0, equal_to_count = 0, greater_than_count = 0;
    double max_of_elements_less_than_pivot = -std::numeric_limits<double>::max();
    double min_of_elements_greater_than_pivot = std::numeric_limits<double>::max();
    for (Iterator run = begin; run != end; ++run)
    {
        if (static_cast<double>(*run) < pivot)
        {
            max_of_elements_less_than_pivot = std::max(max_of_elements_less_than_pivot, static_cast<double>(*run));
            ++less_than_count;
        }
        else if (static_cast<double>(*run) > pivot)
        {
            min_of_elements_greater_than_pivot =
                std::min(min_of_elements_greater_than_pivot, static_cast<double>(*run));
            ++greater_than_count;
        }
        else
        {
            ++equal_to_count;
        }
    }

    performance_stats.increment_pivot_count();
    performance_stats.add_comparisons(2 * less_than_count + 2 * equal_to_count + 3 * greater_than_count);
    return std::tuple<int, int, int, double, double>(less_than_count,
                                                     equal_to_count,
                                                     greater_than_count,
                                                     max_of_elements_less_than_pivot,
                                                     min_of_elements_greater_than_pivot);
}

/**
* The function read_only_numerical_quick_median forwards to this "internal" function.
*/
template <typename Iterator, typename PivotCalculator, typename PerformanceStats>
double read_only_numerical_quick_median_internal(Iterator begin,
                                                 Iterator end,
                                                 PivotCalculator &pivot_calculator,
                                                 PerformanceStats &performance_stats)
{
    if (begin == end)
    {
        throw std::runtime_error("Cannot calculate the median of an empty set.");
    }

    /*
     * The algorithm is an optimized version of the brute force algorithm, taking its
     * cues from quickselect, which in turn is based on quicksort.
     *
     * Quicksort works by selecting a pivot, putting it into its final position in
     * the sort order, then recursively applying itself to the two subsequences on
     * either side of the pivot. Quickselect proceeds in the same way, except that
     * it does not have to apply itself to both subsequences: it is obvious which
     * subsequence contains the n-th element. (This optimization is what causes
     * quickselect to have average complexity O(N) rather than O(N logN).) One way
     * of looking at that is to say when quickselect deals with a pivot, it may not
     * have found the n-th element, but it has found an upper or lower bound, as the
     * case may be, for the n-th element. Obviously, a read-only version of quickselect
     * can do the same thing. The difference is that it cannot pass to geometrically
     * shrinking subsequences, because it cannot do any sorting. Instead, it has to
     * keep track of the current upper and lower bounds. It then makes use of that
     * information by skipping all elements that are not within these bounds when it
     * comes to selecting a pivot.
     *
     * The description above is identical to the one that we gave for the non-numerical
     * case. What's different here (that is, in the numerical case), is this:
     * instead of selecting a pivot from the available elements between the current
     * lower and upper bound, we *calculate* a pivot, namely, as the mean of the current
     * lower and upper bound. This reflects the assumption that the median won't be too
     * far from the mean.
     */

    // If the number of elements is even, the median is an interval of
    // which both ends must be found.
    //
    bool median_interval_left_endpoint_found = false;
    double median_interval_left_endpoint = 0.;
    bool median_interval_right_endpoint_found = false;
    double median_interval_right_endpoint = 0.;

    // As a minor optimization, we trim the sequence at the beginning and
    // the end, discarding elements that are not candidates for a pivot.
    //
    int total_length_of_sequence = 0;
    Iterator active_sequence_begin = begin;
    Iterator active_sequence_end = end;
    int num_discarded_elements_less_than_median_lower_bound = 0;
    int num_discarded_elements_greater_than_median_upper_bound = 0;

    pivot_calculator.initialize(begin, end, performance_stats);
    total_length_of_sequence = pivot_calculator.get_total_sequence_length();
    double median_lower_bound = pivot_calculator.get_total_sequence_min();
    double median_upper_bound = pivot_calculator.get_total_sequence_max();
    performance_stats.set_sequence_length(total_length_of_sequence);

    // To determine the numerical pivot, we need the number of elements below and above
    // the current median lower bound and upper bound, respectively. By keeping that count
    // across loop iterations, we can update it efficiently, that is, without additional
    // comparisons.
    //
    int num_elements_less_than_median_lower_bound = 0;
    int num_elements_greater_than_median_upper_bound = 0;

    // Main loop for selecting and processing pivots.
    //
    double median = 0.0;
    while (true)
    {
        /*
         * Trim the sequence on both sides. This is a minor optimization with
         * little effect on the running time, but it is too obvious to not do it.
         * Trimming at the end can only be done for bidirectional iterators or
         * better, hence the two separate function calls.
         */

        // Trim the sequence at the beginning.
        //
        std::tuple<Iterator, int, int> trim_left_result =
            trim_sequence_left(active_sequence_begin, median_lower_bound, median_upper_bound, performance_stats);
        //
        // Trim the sequence at the end.
        //
        std::tuple<Iterator, int, int> trim_right_result =
            trim_sequence_right(active_sequence_end,
                                median_lower_bound,
                                median_upper_bound,
                                performance_stats,
                                typename std::iterator_traits<Iterator>::iterator_category());
        //
        active_sequence_begin = std::get<0>(trim_left_result);
        active_sequence_end = std::get<0>(trim_right_result);
        num_discarded_elements_less_than_median_lower_bound +=
            (std::get<1>(trim_left_result) + std::get<1>(trim_right_result));
        num_discarded_elements_greater_than_median_upper_bound +=
            (std::get<2>(trim_left_result) + std::get<2>(trim_right_result));

        /*
         * Calculate the pivot and count the number of elements less than, equal to, and greater
         * than the pivot in the subsequence [run, end).
         */

        double pivot = pivot_calculator(median_lower_bound,
                                        median_upper_bound,
                                        num_elements_less_than_median_lower_bound,
                                        num_elements_greater_than_median_upper_bound);

        std::tuple<int, int, int, double, double> element_counts =
            count_elements(active_sequence_begin, active_sequence_end, pivot, performance_stats);

        int num_elements_less_than_pivot =
            std::get<0>(element_counts) + num_discarded_elements_less_than_median_lower_bound;
        int num_elements_equal_to_pivot = std::get<1>(element_counts);
        int num_elements_greater_than_pivot =
            std::get<2>(element_counts) + num_discarded_elements_greater_than_median_upper_bound;
        double max_of_less_than_pivot = std::get<3>(element_counts);
        double min_of_greater_than_pivot = std::get<4>(element_counts);

        /*
         * Check what we found: new lower bound, new upper bound, median position, left or right endpoint of median
         * interval.
         */

        // Too many elements above pivot: new lower bound found.
        //
        if (num_elements_greater_than_pivot > (total_length_of_sequence / 2) + 1)
        {
            median_lower_bound = min_of_greater_than_pivot;
            num_elements_less_than_median_lower_bound = num_elements_less_than_pivot + num_elements_equal_to_pivot;
        }
        //
        // Too many elements below pivot: new upper bound found.
        //
        else if (num_elements_less_than_pivot > (total_length_of_sequence / 2) + 1)
        {
            median_upper_bound = max_of_less_than_pivot;
            num_elements_greater_than_median_upper_bound =
                num_elements_greater_than_pivot + num_elements_equal_to_pivot;
        }
        //
        // Median position or median interval endpoint(s) found. Here, the cases of even and odd number of
        // elements diverge.
        //
        else
        {
            // Sequence length is an odd number: median found.
            //
            if (total_length_of_sequence % 2 == 1)
            {
                if (num_elements_greater_than_pivot == (total_length_of_sequence / 2) + 1)
                {
                    median = min_of_greater_than_pivot;
                }
                else if (num_elements_less_than_pivot == (total_length_of_sequence / 2) + 1)
                {
                    median = max_of_less_than_pivot;
                }
                else
                {
                    median = pivot;
                }
                break;
            }
            //
            // Sequence length is an even number: one or both median interval endpoints found.
            //
            // In the case distinctions below, p stands for the pivot, and e_n stands for the
            // element at 1-based position n. So e_{n/2} is the left endpoint of the median
            // interval.
            //
            else
            {
                // e_{n/2 - 1} != e_{n/2} and p in [e_{n/2 - 1}, e_{n/2})
                //
                if (num_elements_greater_than_pivot == (total_length_of_sequence / 2) + 1)
                {
                    median_interval_left_endpoint_found = true;
                    median_interval_left_endpoint = min_of_greater_than_pivot;

                    // We may have to continue in search of the other median interval endpoint, so also record this
                    // as a lower bound.
                    median_lower_bound = median_interval_left_endpoint;
                    num_elements_less_than_median_lower_bound =
                        num_elements_less_than_pivot + num_elements_equal_to_pivot;
                }
                //
                // e_{n/2} != e_{n/2 + 1} and p == e_{n/2}
                //
                else if (num_elements_greater_than_pivot == total_length_of_sequence / 2 &&
                         num_elements_equal_to_pivot > 0)
                {
                    median_interval_left_endpoint_found = true;
                    median_interval_left_endpoint = pivot;

                    median_interval_right_endpoint_found = true;
                    median_interval_right_endpoint = min_of_greater_than_pivot;
                }
                //
                // e_{n/2} != e_{n/2 + 1} and p in (e_{n/2}, e_{n/2 + 1})
                //
                else if (num_elements_greater_than_pivot == total_length_of_sequence / 2)
                {
                    median_interval_left_endpoint_found = true;
                    median_interval_left_endpoint = max_of_less_than_pivot;

                    median_interval_right_endpoint_found = true;
                    median_interval_right_endpoint = min_of_greater_than_pivot;
                }
                //
                // e_{n/2} != e_{n/2 + 1} and p == e_{n/2 + 1}
                //
                else if (num_elements_less_than_pivot == total_length_of_sequence / 2 &&
                         num_elements_equal_to_pivot > 0)
                {
                    median_interval_left_endpoint_found = true;
                    median_interval_left_endpoint = max_of_less_than_pivot;

                    median_interval_right_endpoint_found = true;
                    median_interval_right_endpoint = pivot;
                }
                //
                // e_{n/2 + 1} != e_{n/2 + 2} and p in (e_{n/2 + 1}, e_{n/2 + 2}]
                //
                else if (num_elements_less_than_pivot == (total_length_of_sequence / 2) + 1)
                {
                    median_interval_right_endpoint_found = true;
                    median_interval_right_endpoint = max_of_less_than_pivot;

                    // We may have to continue in search of the other median interval endpoint, so also record this
                    // as an upper bound.
                    median_upper_bound = median_interval_right_endpoint;
                    num_elements_greater_than_median_upper_bound =
                        num_elements_greater_than_pivot + num_elements_equal_to_pivot;
                }
                //
                // Left and right endpoint of median interval must have been the same, and that was the pivot.
                //
                else
                {
                    assert(num_elements_equal_to_pivot > 0);
                    median = pivot;
                    break;
                }

                if (median_interval_left_endpoint_found && median_interval_right_endpoint_found)
                {
                    median = mean(median_interval_left_endpoint, median_interval_right_endpoint);
                    break;
                }
            }
        }
    }

    performance_stats.update_averages();
    return median;
}
} // end namespace read_only_numerical_quick_median_detail
} // end namespace median_project

#endif // TMB_READ_ONLY_NUMERICAL_QUICK_MEDIAN_DETAIL_08_24_2015_HPP
