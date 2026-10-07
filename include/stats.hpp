// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

/// basic statistics functions for containers of floats
/**
* \file
* \author Steven Ward
*/

#pragma once

#include "container.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <iterator>
#include <limits>
#include <type_traits>
#include <utility>

// https://en.cppreference.com/w/cpp/named_req/Compare
constexpr auto compare_abs_less = [](const auto& a, const auto& b)
{
    return std::abs(a) < std::abs(b);
};

/// a container that keeps its values in ascending order, such as \c std::multiset
/**
* Its first and last values are its minimum and maximum, so they can be read
* in constant time instead of found by a linear search.
*/
template <typename C>
concept sorted_ascending =
    std::same_as<typename C::key_type, typename C::value_type> &&
    (std::same_as<typename C::key_compare, std::less<typename C::key_type>> ||
     std::same_as<typename C::key_compare, std::less<>>);

#define POW2(x) ((x) * (x))
#define POW3(x) ((x) * (x) * (x))
#define POW4(x) ((x) * (x) * (x) * (x))

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
min_val(const Container& c)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n == 0)
        return std::numeric_limits<T>::quiet_NaN();

    if constexpr (sorted_ascending<Container>)
        return *c.cbegin();
    else
        return *std::min_element(c.cbegin(), c.cend());
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
min_abs_val(const Container& c)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n == 0)
        return std::numeric_limits<T>::quiet_NaN();

    return *std::min_element(c.cbegin(), c.cend(), compare_abs_less);
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
max_val(const Container& c)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n == 0)
        return std::numeric_limits<T>::quiet_NaN();

    if constexpr (sorted_ascending<Container>)
        return *c.crbegin();
    else
        return *std::max_element(c.cbegin(), c.cend());
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
max_abs_val(const Container& c)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n == 0)
        return std::numeric_limits<T>::quiet_NaN();

    return *std::max_element(c.cbegin(), c.cend(), compare_abs_less);
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
minmax_vals(const Container& c)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n == 0)
        return std::make_pair(std::numeric_limits<T>::quiet_NaN(),
                              std::numeric_limits<T>::quiet_NaN());

    if constexpr (sorted_ascending<Container>)
    {
        return std::make_pair(*c.cbegin(), *c.crbegin());
    }
    else
    {
        const auto& [min_iter, max_iter] = std::minmax_element(c.cbegin(), c.cend());
        return std::make_pair(*min_iter, *max_iter);
    }
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
minmax_abs_vals(const Container& c)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n == 0)
        return std::make_pair(std::numeric_limits<T>::quiet_NaN(),
                              std::numeric_limits<T>::quiet_NaN());

    const auto& [min_iter, max_iter] =
        std::minmax_element(c.cbegin(), c.cend(), compare_abs_less);
    return std::make_pair(*min_iter, *max_iter);
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
sum_val(const Container& c)
{
    using T = Container::value_type;

    T sum{};
    for (const auto& x : c)
    {
        sum += x;
    }
    return sum;
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
sum_abs_val(const Container& c)
{
    using T = Container::value_type;

    T sum{};
    for (const auto& x : c)
    {
        sum += std::abs(x);
    }
    return sum;
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
arithmetic_mean_val(const Container& c)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n == 0)
        return std::numeric_limits<T>::quiet_NaN();

    return sum_val(c) / n;
}

// adapted from datamash
// https://git.savannah.gnu.org/cgit/datamash.git/tree/src/utils.c#n119

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
variance_val(const Container& c, const bool is_sample = false)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n <= 0)
        return std::numeric_limits<T>::quiet_NaN();

    if (is_sample && n <= 1)
        return std::numeric_limits<T>::quiet_NaN();

    const auto μ = arithmetic_mean_val(c);

    T m2{};
    for (const auto& x : c)
    {
        m2 += POW2(x - μ);
    }

    T variance = m2;

    if (is_sample)
        variance /= (n - 1);
    else
        variance /= n;

    return variance;
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
stdev_val(const Container& c, const bool is_sample = false)
{
    return std::sqrt(variance_val(c, is_sample));
}

// adapted from datamash
// https://git.savannah.gnu.org/cgit/datamash.git/tree/src/utils.c#n209

// https://brownmath.com/stat/shape.htm#Skewness
template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
skewness_val(const Container& c, const bool is_sample = false)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n <= 1)
        return std::numeric_limits<T>::quiet_NaN();

    if (is_sample && n <= 2)
        return std::numeric_limits<T>::quiet_NaN();

    const auto μ = arithmetic_mean_val(c);

    T m2{};
    T m3{};
    for (const auto& x : c)
    {
        m2 += POW2(x - μ);
        m3 += POW3(x - μ);
    }
    m2 /= n;
    m3 /= n;

    T skewness = m3 / std::sqrt(m2 * m2 * m2);

    if (is_sample)
        skewness *= std::sqrt(n * (n - 1)) / (n - 2);

    return skewness;
}

// adapted from datamash
// https://git.savannah.gnu.org/cgit/datamash.git/tree/src/utils.c#n269

// https://brownmath.com/stat/shape.htm#Kurtosis
template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
excess_kurtosis_val(const Container& c, const bool is_sample = false)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n <= 1)
        return std::numeric_limits<T>::quiet_NaN();

    if (is_sample && n <= 3)
        return std::numeric_limits<T>::quiet_NaN();

    const auto μ = arithmetic_mean_val(c);

    T m2{};
    T m4{};
    for (const auto& x : c)
    {
        m2 += POW2(x - μ);
        m4 += POW4(x - μ);
    }
    m2 /= n;
    m4 /= n;

    T excess_kurtosis = m4 / (m2 * m2) - 3;

    if (is_sample)
        excess_kurtosis =
            ((n + 1) * excess_kurtosis + 6) * static_cast<T>(n - 1) / ((n - 2) * (n - 3));

    return excess_kurtosis;
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
median_val(const Container& c)
{
    using T = Container::value_type;
    const auto n = c.size();

    if (n == 0)
        return std::numeric_limits<T>::quiet_NaN();

    if (n % 2 == 0) // even
        return (*std::next(c.cbegin(), n / 2 - 1) + *std::next(c.cbegin(), n / 2)) / 2;
    else // odd
        return *std::next(c.cbegin(), n / 2);
}

template <container Container>
requires std::is_floating_point_v<typename Container::value_type>
auto
range_val(const Container& c)
{
    const auto& [min_val, max_val] = minmax_vals(c);
    return max_val - min_val;
}

#undef POW2
#undef POW3
#undef POW4
