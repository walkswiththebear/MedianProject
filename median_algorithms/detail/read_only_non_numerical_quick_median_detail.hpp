//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_READ_ONLY_NON_NUMERICAL_QUICK_MEDIAN_DETAIL_07_21_2015_HPP
#define TMB_READ_ONLY_NON_NUMERICAL_QUICK_MEDIAN_DETAIL_07_21_2015_HPP

/*
 * Implementation details for read-only non-numerical median algorithms.
 */

#include <tuple>
#include <random>
#include <assert.h>

namespace median_project
{
namespace read_only_non_numerical_quick_median_detail
{

/**
 * Function that trims a sequence at beginning.
 */
template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int>
trim_sequence_left(Iterator active_sequence_begin,
                   bool median_lower_bound_found,
                   typename std::iterator_traits<Iterator>::value_type median_lower_bound,
                   bool median_upper_bound_found,
                   typename std::iterator_traits<Iterator>::value_type median_upper_bound,
                   PerformanceStats &performance_stats)
{
    int num_discarded_elements_less_than_or_equal_to_median_lower_bound = 0;
    int num_discarded_elements_greater_than_or_equal_to_median_upper_bound = 0;
    while (true)
    {
        if (median_lower_bound_found && *active_sequence_begin <= median_lower_bound)
        {
            performance_stats.add_comparisons(1);
            ++num_discarded_elements_less_than_or_equal_to_median_lower_bound;
            ++active_sequence_begin;
        }
        else if (median_upper_bound_found && *active_sequence_begin >= median_upper_bound)
        {
            performance_stats.add_comparisons(2);
            ++num_discarded_elements_greater_than_or_equal_to_median_upper_bound;
            ++active_sequence_begin;
        }
        else
        {
            performance_stats.add_comparisons(2);
            break;
        }
    }

    return std::tuple<Iterator, int, int>(active_sequence_begin,
                                          num_discarded_elements_less_than_or_equal_to_median_lower_bound,
                                          num_discarded_elements_greater_than_or_equal_to_median_upper_bound);
}

/**
 * Function that trims a sequence at end. This is a minor optimization that can be applied
 * only with bidirectional iterators or better.
 */
template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int>
trim_sequence_right(Iterator active_sequence_end,
                    bool median_lower_bound_found,
                    typename std::iterator_traits<Iterator>::value_type median_lower_bound,
                    bool median_upper_bound_found,
                    typename std::iterator_traits<Iterator>::value_type median_upper_bound,
                    PerformanceStats &performance_stats,
                    std::forward_iterator_tag)
{
    return std::tuple<Iterator, int, int>(active_sequence_end, 0, 0);
}
//
template <typename Iterator, typename PerformanceStats>
std::tuple<Iterator, int, int>
trim_sequence_right(Iterator active_sequence_end,
                    bool median_lower_bound_found,
                    typename std::iterator_traits<Iterator>::value_type median_lower_bound,
                    bool median_upper_bound_found,
                    typename std::iterator_traits<Iterator>::value_type median_upper_bound,
                    PerformanceStats &performance_stats,
                    std::bidirectional_iterator_tag)
{
    Iterator active_sequence_last = active_sequence_end;
    --active_sequence_last;

    int num_discarded_elements_less_than_or_equal_to_median_lower_bound = 0;
    int num_discarded_elements_greater_than_or_equal_to_median_upper_bound = 0;
    while (true)
    {
        if (median_lower_bound_found && *active_sequence_last <= median_lower_bound)
        {
            performance_stats.add_comparisons(1);
            ++num_discarded_elements_less_than_or_equal_to_median_lower_bound;
            --active_sequence_last;
        }
        else if (median_upper_bound_found && *active_sequence_last >= median_upper_bound)
        {
            performance_stats.add_comparisons(2);
            ++num_discarded_elements_greater_than_or_equal_to_median_upper_bound;
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
                                          num_discarded_elements_less_than_or_equal_to_median_lower_bound,
                                          num_discarded_elements_greater_than_or_equal_to_median_upper_bound);
}

/**
* Standard pivoting strategy: pick a random element in the active
* sequence, then from there, look for the nearest element that lies
* strictly between the current minimum and maximum. There is a 
* version for bidirectional iterators or better, and one for forward 
* iterators.
* 
* NOTE: The first call to the pivoting strategy will be made with
* an incorrect value of length_of_sequence, namely, 0. We just return
* the first element of the sequence in that case.
*/
class standard_pivoting_strategy
{
public:

    /*
     * Constructor from begin and end iterator of sequence. The iterators
     * are needed only to determine the length of the sequence. If the 
     * length of the sequence is known, and the iterators are not random
     * access, it is better to use the constructor from the length of the
     * sequence.
     */
    template <typename Iterator>
    standard_pivoting_strategy(Iterator begin, Iterator end) :
        standard_pivoting_strategy(std::distance(begin, end)){}

    standard_pivoting_strategy(int total_length_of_sequence) :
        m_generator(std::random_device()()),
        m_pivot_position_generator(0, total_length_of_sequence - 1){
    }

    /*
     * Version for bidirectional iterator or better.
     *
     * NOTE: The first call to the pivoting strategy will be made with
     * an incorrect value of length_of_sequence, namely, 0. We just return
     * the first element of the sequence in that case.
     */
    template <typename Iterator, typename PerformanceStats>
    Iterator operator()(Iterator begin,
                        Iterator end,
                        bool lower_bound_found,
                        typename std::iterator_traits<Iterator>::value_type lower_bound,
                        bool upper_bound_found,
                        typename std::iterator_traits<Iterator>::value_type upper_bound,
                        int length_of_sequence,
                        PerformanceStats &performance_stats,
                        std::bidirectional_iterator_tag) 
    {

        // Attention: read NOTE above.
        //
        if (length_of_sequence == 0 || length_of_sequence == 1)
        {
            return begin;
        }

        std::reverse_iterator<Iterator> reverse_end(begin);
        int pivot_search_start_index = m_pivot_position_generator(m_generator) % length_of_sequence;
        Iterator pivot_position_locator_forward = begin;
        std::advance(pivot_position_locator_forward, pivot_search_start_index);
        std::reverse_iterator<Iterator> pivot_position_locator_backward(pivot_position_locator_forward);

        // Find the pivot candidate (that is, an element that's strictly between
        // the bounds) that's closest to the midpoint.
        //
        while (true)
        {
            if (pivot_position_locator_forward != end)
            {
                if ((lower_bound_found && *pivot_position_locator_forward <= lower_bound) ||
                    (upper_bound_found && *pivot_position_locator_forward >= upper_bound))
                {
                    performance_stats.add_comparisons(1);
                    ++pivot_position_locator_forward;
                }
                else
                {
                    performance_stats.add_comparisons(2);
                    return pivot_position_locator_forward;
                }
            }

            if (pivot_position_locator_backward != reverse_end)
            {
                if ((lower_bound_found && *pivot_position_locator_backward <= lower_bound) ||
                    (upper_bound_found && *pivot_position_locator_backward >= upper_bound))
                {
                    performance_stats.add_comparisons(1);
                    --pivot_position_locator_backward;
                }
                else
                {
                    performance_stats.add_comparisons(2);
                    return --pivot_position_locator_backward.base();
                }
            }
        }
    }

    /*
     * Version for forward iterators.
     *
     * NOTE: The first call to the pivoting strategy will be made with
     * an incorrect value of length_of_sequence, namely, 0. We just return
     * the first element of the sequence in that case.
     */
    template <typename Iterator, typename PerformanceStats>
    Iterator operator()(Iterator begin,
                        Iterator end,
                        bool lower_bound_found,
                        typename std::iterator_traits<Iterator>::value_type lower_bound,
                        bool upper_bound_found,
                        typename std::iterator_traits<Iterator>::value_type upper_bound,
                        int length_of_sequence,
                        PerformanceStats &performance_stats,
                        std::forward_iterator_tag)
    {
        Iterator pivot_candiate_left = begin;
        int num_steps_past_last_pivot_candidate = -1;
        int num_steps_past_desired_pivot_index = -1;
        int current_index = 0;
        bool on_or_past_desired_pivot_index = false;

        // Attention: read "NOTE" above.
        //  
        if (length_of_sequence == 0)
        {
            return begin;
        }

        // Here, using a random index seems to perform noticeably worse than
        // always starting the search in the middle of the sequence.
        int desired_pivot_index = length_of_sequence / 2;

        // Find the pivot candidate (that is, an element that's strictly between
        // the bounds) that's closest to the desired index.
        //
        for (Iterator run = begin; run != end; ++run, ++current_index)
        {
            on_or_past_desired_pivot_index = current_index >= desired_pivot_index;

            if ((lower_bound_found && *run <= lower_bound) || (upper_bound_found && *run >= upper_bound))
            {
                performance_stats.add_comparisons(1);
                if (on_or_past_desired_pivot_index)
                {
                    ++num_steps_past_desired_pivot_index;
                }
                else if (num_steps_past_last_pivot_candidate >= 0)
                {
                    ++num_steps_past_last_pivot_candidate;
                }
            }
            else
            {
                performance_stats.add_comparisons(1);

                if (on_or_past_desired_pivot_index)
                {
                    if (num_steps_past_last_pivot_candidate >= 0)
                    {
                        ++num_steps_past_last_pivot_candidate;
                    }

                    if (num_steps_past_last_pivot_candidate == -1 ||
                        num_steps_past_desired_pivot_index <= num_steps_past_last_pivot_candidate)
                    {
                        return run;
                    }
                    else
                    {
                        return pivot_candiate_left;
                    }
                }
                else
                {
                    pivot_candiate_left = run;
                    num_steps_past_last_pivot_candidate = 0;
                }
            }
        }

        // The loop above returns only after it has moved past the midpoint of the
        // sequence. If no pivot candidate was found there, there must have been one
        // before the midpoint.
        //
        return pivot_candiate_left;
    }

    private:
        std::mt19937 m_generator;
        std::uniform_int_distribution<> m_pivot_position_generator;
};

/**
* Pivoting strategy for random (non-sorted) data: use the first
* suitable candidate.
*/
class pivoting_strategy_for_random_data
{
  public:
    template <typename Iterator, typename PerformanceStats>
    Iterator operator()(Iterator begin,
                        Iterator end,
                        bool lower_bound_found,
                        typename std::iterator_traits<Iterator>::value_type lower_bound,
                        bool upper_bound_found,
                        typename std::iterator_traits<Iterator>::value_type upper_bound,
                        int length_of_sequence,
                        PerformanceStats &performance_stats,
                        std::forward_iterator_tag) const
    {
        // Find the first pivot candidate (that is, an element that's strictly between
        // the bounds).
        //
        for (Iterator run = begin; run != end; ++run)
        {
            performance_stats.add_comparisons(1);

            if ((!lower_bound_found || *run > lower_bound) && (!upper_bound_found || *run < upper_bound))
            {
                return run;
            }
        }

        assert(false);
        return begin;
    }
};

/**
* Counts the number of elements less than, greater than, and equal to a given element in the sequence
* [begin, end).
*/
template <typename Iterator, typename PerformanceStats>
std::tuple<int, int, int> count_elements(Iterator begin,
                                         Iterator end,
                                         typename std::iterator_traits<Iterator>::value_type pivot,
                                         PerformanceStats &performance_stats)
{
    int less_than_count = 0, equal_to_count = 0, greater_than_count = 0;
    for (Iterator run = begin; run != end; ++run)
    {
        if (*run < pivot)
        {
            ++less_than_count;
        }
        else if (*run > pivot)
        {
            ++greater_than_count;
        }
        else
        {
            ++equal_to_count;
        }
    }

    performance_stats.increment_pivot_count();
    performance_stats.add_comparisons(less_than_count + 2 * equal_to_count + 2 * greater_than_count);
    return std::tuple<int, int, int>(less_than_count, equal_to_count, greater_than_count);
}

/**
* The functions read_only_non_numerical_quick_median and
* read_only_non_numerical_quick_median_random_data forward to this "internal"
* function.
*/
template <typename Iterator, typename PivotingStrategy, typename PerformanceStats>
std::pair<Iterator, Iterator> read_only_non_numerical_quick_median_internal(Iterator begin,
                                                                            Iterator end,
                                                                            PivotingStrategy pivoting_strategy,
                                                                            PerformanceStats &performance_stats)
{
    typedef typename std::iterator_traits<Iterator>::value_type value_type;

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
    */

    bool median_lower_bound_found = false;
    value_type median_lower_bound = value_type();
    int num_discarded_elements_less_than_or_equal_to_median_lower_bound = 0;
    //
    bool median_upper_bound_found = false;
    value_type median_upper_bound = value_type();
    int num_discarded_elements_greater_than_or_equal_to_median_upper_bound = 0;

    // If the number of elements is even, the median is an interval of
    // which both ends must be found.
    //
    bool median_interval_left_endpoint_found = false;
    Iterator median_interval_left_endpoint_pos = begin;
    bool median_interval_right_endpoint_found = false;
    Iterator median_interval_right_endpoint_pos = end;

    // Loop for selecting and processing pivots. Each pivoting either
    // finds the median, in which case we're done, or it gives us a new
    // lower or upper bound, as the case may be, for the median.
    //
    int total_length_of_sequence = 0;
    Iterator active_sequence_begin = begin;
    Iterator active_sequence_end = end;
    std::pair<Iterator, Iterator> median_pos_pair(begin, end);
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
        std::tuple<Iterator, int, int> trim_left_result = trim_sequence_left(active_sequence_begin,
                                                                             median_lower_bound_found,
                                                                             median_lower_bound,
                                                                             median_upper_bound_found,
                                                                             median_upper_bound,
                                                                             performance_stats);

        // Trim the sequence at the end.
        //
        std::tuple<Iterator, int, int> trim_right_result =
            trim_sequence_right(active_sequence_end,
                                median_lower_bound_found,
                                median_lower_bound,
                                median_upper_bound_found,
                                median_upper_bound,
                                performance_stats,
                                typename std::iterator_traits<Iterator>::iterator_category());

        active_sequence_begin = std::get<0>(trim_left_result);
        active_sequence_end = std::get<0>(trim_right_result);
        num_discarded_elements_less_than_or_equal_to_median_lower_bound +=
            (std::get<1>(trim_left_result) + std::get<1>(trim_right_result));
        num_discarded_elements_greater_than_or_equal_to_median_upper_bound +=
            (std::get<2>(trim_left_result) + std::get<2>(trim_right_result));

        /*
         * Select a pivot and count the number of elements less than, equal to, and greater
         * than pivot in the subsequence [run, end).
         *
         * NOTE: The argument total_length_of_sequence is actually 0 during the
         * first iteration. That's fine. It just means that the first pivot is
         * the first element of the sequence.
         */

        Iterator pivot_pos = pivoting_strategy(active_sequence_begin,
                                               active_sequence_end,
                                               median_lower_bound_found,
                                               median_lower_bound,
                                               median_upper_bound_found,
                                               median_upper_bound,
                                               total_length_of_sequence -
                                                   num_discarded_elements_less_than_or_equal_to_median_lower_bound -
                                                   num_discarded_elements_greater_than_or_equal_to_median_upper_bound,
                                               performance_stats,
                                               typename std::iterator_traits<Iterator>::iterator_category());

        std::tuple<int, int, int> element_counts =
            count_elements(active_sequence_begin, active_sequence_end, *pivot_pos, performance_stats);

        int num_elements_less_than_pivot =
            std::get<0>(element_counts) + num_discarded_elements_less_than_or_equal_to_median_lower_bound;
        int num_elements_equal_to_pivot = std::get<1>(element_counts);
        int num_elements_greater_than_pivot =
            std::get<2>(element_counts) + num_discarded_elements_greater_than_or_equal_to_median_upper_bound;
        assert(num_elements_equal_to_pivot >= 1);

        // Store the element count. (This may happen several times.)
        //
        if (active_sequence_begin == begin && active_sequence_end == end)
        {
            total_length_of_sequence =
                num_elements_less_than_pivot + num_elements_equal_to_pivot + num_elements_greater_than_pivot;
            performance_stats.set_sequence_length(total_length_of_sequence);
        }
        else
        {
            assert(total_length_of_sequence ==
                   num_elements_less_than_pivot + num_elements_equal_to_pivot + num_elements_greater_than_pivot);
        }

        /*
         * Check what we found: new lower bound, new upper bound, median position, left or right endpoint of
         * median
         * interval.
         */

        // Too many elements above pivot: new lower bound found.
        //
        if (num_elements_greater_than_pivot > total_length_of_sequence / 2)
        {
            median_lower_bound_found = true;
            median_lower_bound = *pivot_pos;
        }
        //
        // Too many elements below pivot: new upper bound found.
        //
        else if (num_elements_less_than_pivot > total_length_of_sequence / 2)
        {
            median_upper_bound_found = true;
            median_upper_bound = *pivot_pos;
        }
        //
        // Median position or median interval endpoint(s) found. Here, the cases of even and odd number of
        // elements diverge.
        //
        else
        {
            // Sequence length is an even number: median found.
            //
            if (total_length_of_sequence % 2 == 1)
            {
                median_pos_pair = std::pair<Iterator, Iterator>(pivot_pos, pivot_pos);
                break;
            }
            //
            // Sequence length is an odd number: one or both median interval endpoints found.
            //
            else
            {

                // Less than half are below: median interval left endpoint found.
                //
                if (num_elements_less_than_pivot < total_length_of_sequence / 2)
                {
                    median_interval_left_endpoint_pos = pivot_pos;
                    median_interval_left_endpoint_found = true;

                    // We may have to continue in search of the other median interval endpoint, so also record
                    // this as a lower bound.
                    median_lower_bound_found = true;
                    median_lower_bound = *pivot_pos;
                }

                // Less than half are above: median interval right endpoint found.
                //
                // NOTE: This is not an else to the above if. We may have found both
                // median interval endpoints (duplicates!).
                //
                if (num_elements_greater_than_pivot < total_length_of_sequence / 2)
                {
                    median_interval_right_endpoint_pos = pivot_pos;
                    median_interval_right_endpoint_found = true;

                    // We may have to continue in search of the other median interval endpoint, so also record
                    // this as an upper bound.
                    median_upper_bound_found = true;
                    median_upper_bound = *pivot_pos;
                }

                if (median_interval_left_endpoint_found && median_interval_right_endpoint_found)
                {
                    median_pos_pair = std::pair<Iterator, Iterator>(median_interval_left_endpoint_pos,
                                                                    median_interval_right_endpoint_pos);
                    break;
                }
            }
        }
    }

    performance_stats.update_averages();
    return median_pos_pair;
}
} // end namespace read_only_non_numerical_quick_median_detail
} // end namespace median_project

#endif // TMB_READ_ONLY_NON_NUMERICAL_QUICK_MEDIAN_DETAIL_07_21_2015_HPP
