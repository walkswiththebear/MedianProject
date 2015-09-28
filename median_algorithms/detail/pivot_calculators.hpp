//  (C) Copyright Thomas Becker 2015. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef TMB_NO_OP_PIVOT_CALCULATORS_09_27_2015_HPP
#define TMB_NO_OP_PIVOT_CALCULATORS_09_27_2015_HPP

/*
 * Standard pivot calculators for numerical mean algorithms.
 */

namespace median_project
{
    namespace read_only_numerical_quick_median_detail
    {

        /**
         * Pivot calculator for uniform distributions: use the midpoint of sequence min
         * and max.
         */
        class uniform_distribution_pivot
        {
        public:
            void initialize(int total_sequence_length, double sequence_min,
                double sequence_max, double sequence_mean){};
            double operator() (double median_lower_bound, double median_upper_bound) const {
                return median_lower_bound / 2.0 + median_upper_bound / 2.0;
            }
        };
    } // end namespace read_only_numerical_quick_median_detail
} // end namespace median_project

#endif // TMB_NO_OP_PIVOT_CALCULATORS_09_27_2015_HPP

