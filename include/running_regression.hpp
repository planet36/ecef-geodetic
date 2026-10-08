// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

/// running regression class
/**
* \file
* \author John D. Cook
* \author Steven Ward
* \sa https://www.johndcook.com/blog/running_regression/
*/

#pragma once

#include "running_stats.hpp"

#include <concepts>

template <std::floating_point T = double>
class running_regression
{
private:
    running_stats<T> x_stats;
    running_stats<T> y_stats;
    T S_xy = 0;
    long long n = 0;

public:
    /*
    /// default ctor
    running_regression()
    {
        clear();
    }
    */

    constexpr void clear() noexcept
    {
        x_stats.clear();
        y_stats.clear();
        S_xy = 0;
        n = 0;
    }

    void push(const T x, const T y) noexcept
    {
        // The means of empty stats are NaN, and 0 * NaN is still NaN, so the
        // first value must skip this term rather than rely on n being zero.
        if (n > 0)
            S_xy += n * (x_stats.mean() - x) * (y_stats.mean() - y) / (n + 1);

        x_stats.push(x);
        y_stats.push(y);
        n++;
    }

    [[nodiscard]] constexpr auto num_data_values() const noexcept { return n; }

    [[nodiscard]] constexpr auto slope() const noexcept
    {
        const auto S_xx = x_stats.variance() * (n - 1);
        return S_xy / S_xx;
    }

    [[nodiscard]] constexpr auto intercept() const noexcept
    {
        return y_stats.mean() - slope() * x_stats.mean();
    }

    [[nodiscard]] auto correlation() const noexcept
    {
        const auto t = x_stats.standard_deviation() * y_stats.standard_deviation();
        return S_xy / ((n - 1) * t);
    }

    template <std::floating_point T2>
    friend running_regression<T2> operator+(const running_regression<T2>& a,
                                            const running_regression<T2>& b) noexcept;

    running_regression<T>& operator+=(const running_regression<T>& that) noexcept
    {
        const running_regression<T> combined = *this + that;
        *this = combined;
        return *this;
    }
};

template <std::floating_point T>
[[nodiscard]] running_regression<T>
operator+(const running_regression<T>& a, const running_regression<T>& b) noexcept
{
    // Merging an empty object returns the other one unchanged, because the
    // formula below divides by the combined count.
    if (a.n == 0)
        return b;

    if (b.n == 0)
        return a;

    running_regression<T> combined;

    combined.x_stats = a.x_stats + b.x_stats;
    combined.y_stats = a.y_stats + b.y_stats;
    combined.n = a.n + b.n;

    auto delta_x = b.x_stats.mean() - a.x_stats.mean();
    auto delta_y = b.y_stats.mean() - a.y_stats.mean();
    combined.S_xy = a.S_xy + b.S_xy + a.n * b.n * delta_x * delta_y / combined.n;

    return combined;
}
