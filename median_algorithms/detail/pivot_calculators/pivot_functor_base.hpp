//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_PIVOT_FUNCTOR_BASE_09_28_2015_HPP
#define TMB_PIVOT_FUNCTOR_BASE_09_28_2015_HPP

/*
 * Base class for all pivot functors
 */

namespace median_project
{
namespace read_only_numerical_quick_median_detail
{

/**
 * Abstract base class for all pivot functors: provides default data members for
 * length, min, max, mean, and standard deviation, initialization, and getters for
 * these, and abstract methods.
 */
class pivot_functor_base
{
  public:
    /**
     * Default constructor. Initializes member variables to NaN.
     */
    pivot_functor_base()
        : m_total_sequence_length(std::numeric_limits<int>::quiet_NaN())
        , m_total_sequence_min(std::numeric_limits<double>::quiet_NaN())
        , m_total_sequence_max(std::numeric_limits<double>::quiet_NaN())
        , m_total_sequence_mean(std::numeric_limits<double>::quiet_NaN())
        , m_total_sequence_std_dev(std::numeric_limits<double>::quiet_NaN())
    {
    }

    /**
     * Calculates basic sequence data and stores it in member variables.
     */
    template <typename Iterator, typename PerformanceStats>
    void initialize(Iterator begin, Iterator end, PerformanceStats &performance_stats)
    {
        assert(begin != end);
        int length = 0;
        double min_value = std::numeric_limits<double>::max();
        double max_value = -std::numeric_limits<double>::max();
        double mean = 0.0;
        double variance = 0.0;
        for (Iterator run = begin; run != end; ++run)
        {
            ++length;
            min_value = std::min(min_value, static_cast<double>(*run));
            max_value = std::max(max_value, static_cast<double>(*run));

            double delta = *run - mean;
            mean = mean + delta / static_cast<double>(length);
            variance = variance + delta * (*run - mean);
        }
        variance /= static_cast<double>(length);

        performance_stats.add_comparisons(2 * length);

        m_total_sequence_length = length;
        m_total_sequence_min = min_value;
        m_total_sequence_max = max_value;
        m_total_sequence_mean = mean;
        m_total_sequence_std_dev = sqrt(variance);
    }

    /**
     * Getters for member variables that are used outside of the pivot functor.
     */
    double get_total_sequence_length()
    {
        return m_total_sequence_length;
    }
    //
    double get_total_sequence_min()
    {
        return m_total_sequence_min;
    }
    //
    double get_total_sequence_max()
    {
        return m_total_sequence_max;
    }

    /**
     * Function call operator, to be overridden by derived classes. Calculates the next
     * pivot from the max and min of the remaining interval.
     */
    virtual double operator()(double median_lower_bound, double median_upper_bound) const = 0;

  protected:
    int m_total_sequence_length;
    double m_total_sequence_min;
    double m_total_sequence_max;
    double m_total_sequence_mean;
    double m_total_sequence_std_dev;
};
} // end namespace read_only_numerical_quick_median_detail
} // end namespace median_project

#endif // TMB_PIVOT_FUNCTOR_BASE_09_28_2015_HPP
