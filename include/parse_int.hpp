// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

/// Parse integers from strings, option arguments, and environment variables
/**
* \file
* \author Steven Ward
*/

#pragma once

#if defined(DEBUG)
#include <cassert>
#endif
#include <charconv>
#include <concepts>
#include <cstdlib>
#include <expected>
#include <format>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <system_error>

/// an integer type that \c std::from_chars can parse with a base
template <typename T>
concept from_chars_integral =
    std::integral<T> &&
    requires (const char* p, T& v) { std::from_chars(p, p, v, 10); };

/// Parse all of \a s as an integer of type \c T in <code>[min, max]</code>
/**
* \tparam T the integer type of the result
* \param s the string to parse
* \param min the minimum allowed value (inclusive)
* \param max the maximum allowed value (inclusive)
* \param base the numeric base (2 to 36), or \c 0 to detect a \c "0x"/"0X"
*        (hexadecimal) or \c "0" (octal) prefix like \c strtol
* \retval std::errc::invalid_argument \a s is not entirely an integer
* \retval std::errc::result_out_of_range the value is not representable in
*         \c T or is not in <code>[min, max]</code>
* \return the parsed value, or one of the above error values
* \pre \a min is at most \a max
* \pre \a base is \c 0, or in the interval <code>[2, 36]</code> that
*      \c std::from_chars accepts
* \note Unlike \c std::stoi, this rejects leading whitespace, a leading
*       \c '+', and trailing non-digit characters.
* \note With \a base \c 0, a \c '-' and a prefix interact unlike in
*       \c strtol.  A \c '-' before \c "0x" is rejected, and one before
*       \c "0" hides the prefix, so \c "-010" parses as decimal -10.  A
*       \c '-' after \c "0x" is accepted, so \c "0x-5" parses as -5.
*/
template <from_chars_integral T = int>
[[nodiscard]] constexpr std::expected<T, std::errc>
parse_int(std::string_view s,
          const T min = std::numeric_limits<T>::min(),
          const T max = std::numeric_limits<T>::max(),
          int base = 10)
{
    if (base == 0)
    {
        if (s.starts_with("0x") || s.starts_with("0X"))
        {
            s.remove_prefix(2);
            base = 16;
        }
        else if (s.size() > 1 && s.front() == '0')
        {
            base = 8;
        }
        else
        {
            base = 10;
        }
    }

#if defined(DEBUG)
    assert(min <= max);
    assert((base >= 2) && (base <= 36));
#endif

    T value{};
    const auto [ptr, ec] =
        std::from_chars(s.data(), s.data() + s.size(), value, base);

    if (ec != std::errc{})
        return std::unexpected{ec};

    if (ptr != s.data() + s.size())
        return std::unexpected{std::errc::invalid_argument};

    if (value < min || value > max)
        return std::unexpected{std::errc::result_out_of_range};

    return value;
}

/// Parse \a optarg as an int in <code>[min, max]</code>
/**
* \param optarg the option argument to parse
* \param min the minimum allowed value (inclusive)
* \param max the maximum allowed value (inclusive)
* \param option_name the option name, named by the exception message
* \return the parsed value
* \pre \a min is at most \a max
* \exception std::invalid_argument \a optarg is null, or is not entirely an
*            integer
* \exception std::out_of_range the value is not in <code>[min, max]</code>
*/
[[nodiscard]] inline int
parse_option_int(const char* optarg, const int min, const int max, const char* option_name)
{
    if (optarg == nullptr)
        throw std::invalid_argument{std::format("null argument: {}", option_name)};

    const auto value = parse_int<int>(optarg, min, max);

    if (!value.has_value())
    {
        if (value.error() == std::errc::result_out_of_range)
            throw std::out_of_range{
                std::format("out of range: {}: \"{}\"", option_name, optarg)};

        throw std::invalid_argument{
            std::format("invalid argument: {}: \"{}\"", option_name, optarg)};
    }

    return *value;
}

/// Parse \a optarg as an int
/**
* \param optarg the option argument to parse
* \param option_name the option name, named by the exception message
* \return the parsed value
* \note The parsed value is bounded only by the range of \c int.
* \exception std::invalid_argument \a optarg is null, or is not entirely an
*            integer
* \exception std::out_of_range the value is not representable in \c int
*/
[[nodiscard]] inline int
parse_option_int(const char* optarg, const char* option_name)
{
    return parse_option_int(optarg, std::numeric_limits<int>::min(),
                            std::numeric_limits<int>::max(), option_name);
}

/// Parse the environment variable \a name as an int in <code>[min, max]</code>
/**
* \param name the name of the environment variable
* \param min the minimum allowed value (inclusive)
* \param max the maximum allowed value (inclusive)
* \param default_value the value returned if the variable is not set, which is
*        not checked against <code>[min, max]</code>
* \return the parsed value, or \a default_value if the variable is not set
* \pre \a min is at most \a max
* \exception std::invalid_argument the value is empty, or is not entirely an
*            integer
* \exception std::out_of_range the value is not in <code>[min, max]</code>
*/
[[nodiscard]] inline int
parse_env_int(const char* name, const int min, const int max, const int default_value)
{
    const char* const value = std::getenv(name);
    if (value == nullptr)
        return default_value;

    return parse_option_int(value, min, max, name);
}
