//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_NUMERICAL_QUICK_MEDIAN_DETAIL_10_06_2015_HPP
#define TMB_NUMERICAL_QUICK_MEDIAN_DETAIL_10_06_2015_HPP

/*
 * Implementation detail for quick median algorithm with numerical (distribution-specific) pivoting.
 */

#include <algorithm>
#include <tuple>
#include <limits>
#include <assert.h>

namespace median_project
{
namespace numerical_quick_median_detail
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
* The function numerical_quick_median forwards to this "internal" function.
*/
template <typename RandomAccessIterator, typename PivotCalculator, typename PerformanceStats>
double numerical_quick_median_internal(RandomAccessIterator begin,
                                       RandomAccessIterator end,
                                       PivotCalculator &pivot_calculator,
                                       PerformanceStats &performance_stats)
{
    if (begin == end)
    {
        throw std::runtime_error("Cannot calculate the median of an empty set.");
    }

    pivot_calculator.initialize(begin, end, performance_stats);
    int total_length_of_sequence = pivot_calculator.get_total_sequence_length();
    double current_interval_min = pivot_calculator.get_total_sequence_min();
    double current_interval_max = pivot_calculator.get_total_sequence_max();
    performance_stats.set_sequence_length(total_length_of_sequence);

    RandomAccessIterator current_interval_begin = begin;
    RandomAccessIterator current_interval_end = end;
    int num_elements_below_current_interval = 0;
    int num_elements_above_current_interval = 0;
    double median = 0.0;
    while (true)
    {

        // Select a pivot as close to the median as you can guess by the data distribution
        double pivot = pivot_calculator(current_interval_min,
                                        current_interval_max,
                                        num_elements_below_current_interval,
                                        num_elements_above_current_interval);

        // Swap elements so that everything less than the pivot is to the left of everything that is
        // greater than or equal to the pivot.
        //
        RandomAccessIterator up_runner = current_interval_begin;
        RandomAccessIterator down_runner = (current_interval_end - 1);
        int num_elements_equal_to_pivot = 0;
        double max_of_less_than_pivot = -std::numeric_limits<double>::max();
        double min_of_greater_than_pivot = std::numeric_limits<double>::max();
        // Sweep up until a candidate for swapping is encountered, or the runners meet.
        //
        while (up_runner != down_runner)
        {
            while (static_cast<double>(*up_runner) < pivot && up_runner != down_runner)
            {
                if (static_cast<double>(*up_runner) == pivot)
                {
                    ++num_elements_equal_to_pivot;
                }
                else
                {
                    max_of_less_than_pivot = std::max(max_of_less_than_pivot, *up_runner);
                }
                ++up_runner;
            }

            // Sweep down until a candidate for swapping is encountered, or the runners meet.
            //
            while (static_cast<double>(*down_runner) >= pivot && down_runner != up_runner)
            {
                if (static_cast<double>(*down_runner) == pivot)
                {
                    ++num_elements_equal_to_pivot;
                }
                else
                {
                    min_of_greater_than_pivot = std::min(min_of_greater_than_pivot, *down_runner);
                }
                --down_runner;
            }

            // Swap if the runners haven't met yet.
            //
            if (up_runner != down_runner)
            {
                std::swap(*up_runner, *down_runner);
            }
            //
            // If the two runners are now equal, then this is the last run through the loop,
            // the element currently pointed to by both runners is greater than or equal to the
            // pivot (see exercise below), and that element has not been processed yet for number of
            // elements equal to pivot and min of elements greater than pivot.
            //
            else
            {
                if (*down_runner == pivot)
                {
                    ++num_elements_equal_to_pivot;
                }
                else
                {
                    min_of_greater_than_pivot = std::min(min_of_greater_than_pivot, *down_runner);
                }
            }
        }
        
        /*
         * Exercise: Prove that at this point, both up_runner and down-runner point at the first
         * position of the right-hand subinterval, that is, the one that contains all elements
         * that are greater than or equal to the pivot.
         */

        // Count elements
        //
        int num_elements_in_left_subinterval = up_runner - current_interval_begin;
        int num_elements_less_than_pivot = num_elements_in_left_subinterval + num_elements_below_current_interval;

        int num_elements_in_right_subinterval = current_interval_end - down_runner;
        int num_elements_greater_than_pivot_in_right_subinterval =
            num_elements_in_right_subinterval - num_elements_equal_to_pivot;
        int num_elements_greater_than_pivot =
            num_elements_greater_than_pivot_in_right_subinterval + num_elements_above_current_interval;

        /*
         * Check what we found: median position, left or right endpoint of median
         * interval, or a subinterval to search in.
         */

        // Too many elements below pivot: continue search in left-hand subinterval.
        //
        if (num_elements_less_than_pivot > (total_length_of_sequence + 1) / 2)
        {
            current_interval_end = up_runner;
            current_interval_max = max_of_less_than_pivot;
            num_elements_above_current_interval += num_elements_in_right_subinterval;
        }
        //
        // Too many elements above pivot: continue search in right-hand subinterval.
        //
        else if (num_elements_greater_than_pivot > (total_length_of_sequence + 1) / 2)
        {
            current_interval_begin = up_runner;
            current_interval_min = num_elements_equal_to_pivot > 0 ? pivot : min_of_greater_than_pivot;
            num_elements_below_current_interval += num_elements_in_left_subinterval;
        }
        //
        // Median position or median interval endpoint(s) found. Here, the cases of even and odd number of
        // elements diverge.
        //
        else
        {
            // Sequence length is an odd number: median is the greatest element of the left-hand interval
            // or the least element of the right hand interval.
            //
            if (total_length_of_sequence % 2 == 1)
            {
                if (num_elements_less_than_pivot == (total_length_of_sequence + 1) / 2)
                {
                    median = max_of_less_than_pivot;
                }
                else if (num_elements_greater_than_pivot == (total_length_of_sequence + 1) / 2)
                {
                    median = min_of_greater_than_pivot;
                }
                else
                {
                    median = pivot;
                }
                break;
            }
            //
            // Sequence length is an even number: median interval endpoints found.
            //
            // NOTE: Here, (total_length_of_sequence + 1) / 2 == total_length_of_sequence / 2
            //
            else
            {
                double median_interval_left_endpoint = 0.0;
                double median_interval_right_endpoint = 0.0;
                if (num_elements_less_than_pivot == total_length_of_sequence / 2)
                {
                    median_interval_left_endpoint = max_of_less_than_pivot;
                    median_interval_right_endpoint =
                        num_elements_equal_to_pivot > 0 ? pivot : min_of_greater_than_pivot;
                }
                else if (num_elements_greater_than_pivot == total_length_of_sequence / 2)
                {
                    median_interval_left_endpoint = num_elements_equal_to_pivot > 0 ? pivot : max_of_less_than_pivot;
                    median_interval_right_endpoint = min_of_greater_than_pivot;
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

                median = mean(median_interval_left_endpoint, median_interval_right_endpoint);
                break;

            } // end even/odd distinction

        } // end median found yes/no?

    } // end while-loop

    performance_stats.update_averages();
    return median;
}
} // end namespace numerical_quick_median_detail
} // end namespace median_project

#endif // TMB_NUMERICAL_QUICK_MEDIAN_DETAIL_10_06_2015_HPP
