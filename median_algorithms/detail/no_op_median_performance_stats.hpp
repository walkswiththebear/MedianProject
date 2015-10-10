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

namespace median_project
{

/**
* No-op performance stats class for median, to be used as the default in the algorithms.
*/
class no_op_median_performance_stats
{

    // No-op methods
    // =============

    /*
     * NOTE: Strictly speaking, these methods should be private. But the friends thing
     * was really getting out of hand. Really.
     */

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
