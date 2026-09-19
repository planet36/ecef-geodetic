// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

/// running stats class
/**
* \file
* \author John D. Cook
* \author Steven Ward
* \sa https://www.johndcook.com/blog/skewness_kurtosis/
*
* XXX: Do not compile with -ffinite-math-only (included with -ffast-math (included with -Ofast)).
* This affects the behavior of functions std::fmin, std::fmax.
*/

#pragma once

#include <cmath>
#include <concepts>
#include <iterator>
#include <limits>

template <std::floating_point T = double>
class running_stats
{
private:
    T M1 = 0;
    T M2 = 0;
    T M3 = 0;
    T M4 = 0;
    T _sum = 0;
    T _min = std::numeric_limits<T>::quiet_NaN();
    T _max = std::numeric_limits<T>::quiet_NaN();
    T _sum_abs = 0;
    T _min_abs = std::numeric_limits<T>::quiet_NaN();
    T _max_abs = std::numeric_limits<T>::quiet_NaN();
    long long n = 0;

public:
    /// Reset all statistics to their initial state
    void clear()
    {
        M1 = 0;
        M2 = 0;
        M3 = 0;
        M4 = 0;
        _sum = 0;
        _min = std::numeric_limits<T>::quiet_NaN();
        _max = std::numeric_limits<T>::quiet_NaN();
        _sum_abs = 0;
        _min_abs = std::numeric_limits<T>::quiet_NaN();
        _max_abs = std::numeric_limits<T>::quiet_NaN();
        n = 0;
    }

    /// Add a single value to the running statistics
    void push(const T x)
    {
        const auto n1 = n;
        n++;
        const auto delta = x - M1;
        const auto delta_n = delta / n;
        const auto delta_n2 = delta_n * delta_n;
        const auto term1 = delta * delta_n * n1;
        M1 += delta_n;
        M4 += term1 * delta_n2 * (n * n - 3 * n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
        M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
        M2 += term1;

        _sum += x;
        _min = std::fmin(_min, x);
        _max = std::fmax(_max, x);
        const auto abs_x = std::abs(x);
        _sum_abs += abs_x;
        _min_abs = std::fmin(_min_abs, abs_x);
        _max_abs = std::fmax(_max_abs, abs_x);
    }

    /// Add a range of values to the running statistics
    template <std::forward_iterator It>
    void push(It first, It last)
    {
        for (; first != last; ++first)
        {
            push(*first);
        }
    }

    /// get the number of values pushed
    [[nodiscard]] constexpr auto num_data_values() const { return n; }

    /// get the mean, or NaN if no values were pushed
    [[nodiscard]] constexpr auto mean() const
    {
        return (n > 0) ? M1 : std::numeric_limits<T>::quiet_NaN();
    }

    /// get the sample variance, or NaN if fewer than 2 values were pushed
    [[nodiscard]] auto variance() const
    {
        return (n > 1) ? M2 / (n - 1) : std::numeric_limits<T>::quiet_NaN();
    }

    /// get the sample standard deviation, or NaN if fewer than 2 values were pushed
    [[nodiscard]] auto standard_deviation() const { return std::sqrt(variance()); }

    /// get the skewness, or NaN if fewer than 2 values were pushed
    [[nodiscard]] auto skewness() const { return std::sqrt(n) * M3 / std::pow(M2, 1.5); }

    /// get the kurtosis, or NaN if fewer than 2 values were pushed
    [[nodiscard]] auto kurtosis() const { return n * M4 / (M2 * M2) - 3; }

    /// get the sum of the values, or 0 if none were pushed
    [[nodiscard]] constexpr auto sum() const { return _sum; }

    /// get the minimum value, or NaN if none were pushed
    [[nodiscard]] constexpr auto min() const { return _min; }

    /// get the maximum value, or NaN if none were pushed
    [[nodiscard]] constexpr auto max() const { return _max; }

    /// get the sum of the absolute values, or 0 if none were pushed
    [[nodiscard]] constexpr auto sum_abs() const { return _sum_abs; }

    /// get the minimum absolute value, or NaN if none were pushed
    [[nodiscard]] constexpr auto min_abs() const { return _min_abs; }

    /// get the maximum absolute value, or NaN if none were pushed
    [[nodiscard]] constexpr auto max_abs() const { return _max_abs; }

    template <std::floating_point T2>
    friend running_stats<T2> operator+(const running_stats<T2>& a,
                                       const running_stats<T2>& b);

    /// Merge another \c running_stats into this one
    running_stats<T>& operator+=(const running_stats<T>& that)
    {
        const running_stats<T> combined = *this + that;
        *this = combined;
        return *this;
    }
};

/// Merge two \c running_stats objects into one
template <std::floating_point T>
[[nodiscard]] running_stats<T>
operator+(const running_stats<T>& a, const running_stats<T>& b)
{
    running_stats<T> combined;

    combined.n = a.n + b.n;

    const auto delta = b.M1 - a.M1;
    const auto delta2 = delta * delta;
    const auto delta3 = delta * delta2;
    const auto delta4 = delta2 * delta2;

    combined.M1 = (a.n * a.M1 + b.n * b.M1) / combined.n;

    combined.M2 = a.M2 + b.M2 + delta2 * a.n * b.n / combined.n;

    combined.M3 = a.M3 + b.M3 + delta3 * a.n * b.n * (a.n - b.n) / (combined.n * combined.n);
    combined.M3 += 3 * delta * (a.n * b.M2 - b.n * a.M2) / combined.n;

    combined.M4 = a.M4 + b.M4 +
                  delta4 * a.n * b.n * (a.n * a.n - a.n * b.n + b.n * b.n) /
                      (combined.n * combined.n * combined.n);
    combined.M4 +=
        6 * delta2 * (a.n * a.n * b.M2 + b.n * b.n * a.M2) / (combined.n * combined.n) +
        4 * delta * (a.n * b.M3 - b.n * a.M3) / combined.n;

    // The remaining members combine directly, the same way push() accumulates
    // them.  fmin and fmax carry the empty-object NaN sentinel through, so
    // merging an empty object leaves the other one's extrema.
    combined._sum = a._sum + b._sum;
    combined._min = std::fmin(a._min, b._min);
    combined._max = std::fmax(a._max, b._max);
    combined._sum_abs = a._sum_abs + b._sum_abs;
    combined._min_abs = std::fmin(a._min_abs, b._min_abs);
    combined._max_abs = std::fmax(a._max_abs, b._max_abs);

    return combined;
}
