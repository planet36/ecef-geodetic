// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// running regression class
/**
\file
\author John D. Cook
\author Steven Ward
\sa https://www.johndcook.com/blog/running_regression/
*/

#pragma once

#include "running_stats.hpp"

#include <concepts>

template <std::floating_point T>
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

    void clear()
    {
        x_stats.clear();
        y_stats.clear();
        S_xy = 0;
        n = 0;
    }

    void push(const T x, const T y)
    {
        S_xy += n * (x_stats.mean() - x) * (y_stats.mean() - y) / (n + 1);

        x_stats.push(x);
        y_stats.push(y);
        n++;
    }

    auto num_data_values() const { return n; }

    auto slope() const
    {
        const auto S_xx = x_stats.variance() * (n - 1);
        return S_xy / S_xx;
    }

    auto intercept() const { return y_stats.mean() - slope() * x_stats.mean(); }

    auto correlation() const
    {
        const auto t =
            x_stats.standard_deviation() * y_stats.standard_deviation();
        return S_xy / ((n - 1) * t);
    }

    template <std::floating_point T2>
    friend running_regression<T2> operator+(const running_regression<T2>& a,
                                            const running_regression<T2>& b);

    running_regression<T>& operator+=(const running_regression<T>& that)
    {
        const running_regression<T> combined = *this + that;
        *this = combined;
        return *this;
    }
};

template <std::floating_point T>
running_regression<T>
operator+(const running_regression<T>& a, const running_regression<T>& b)
{
    running_regression<T> combined;

    combined.x_stats = a.x_stats + b.x_stats;
    combined.y_stats = a.y_stats + b.y_stats;
    combined.n = a.n + b.n;

    auto delta_x = b.x_stats.mean() - a.x_stats.mean();
    auto delta_y = b.y_stats.mean() - a.y_stats.mean();
    combined.S_xy =
        a.S_xy + b.S_xy + a.n * b.n * delta_x * delta_y / combined.n;

    return combined;
}
