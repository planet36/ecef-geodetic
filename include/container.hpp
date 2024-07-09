// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// Concept for a container type
/**
\file
\author Steven Ward
\sa https://en.cppreference.com/w/cpp/named_req/Container
\sa https://stackoverflow.com/a/70543294
*/

#pragma once

#include <concepts>
#include <iterator>

template <typename C>
concept container = requires(C a, const C b) {
	requires std::regular<C>;
	requires std::swappable<C>;
	requires std::destructible<typename C::value_type>;
	requires std::same_as<typename C::reference, typename C::value_type&>;
	requires std::same_as<typename C::const_reference,
	                      const typename C::value_type&>;
	requires std::forward_iterator<typename C::iterator>;
	requires std::forward_iterator<typename C::const_iterator>;
	requires std::signed_integral<typename C::difference_type>;
	requires std::same_as<
	    typename C::difference_type,
	    typename std::iterator_traits<typename C::iterator>::difference_type>;
	requires std::same_as<typename C::difference_type,
	                      typename std::iterator_traits<
	                          typename C::const_iterator>::difference_type>;
	{ a.begin() } -> std::same_as<typename C::iterator>;
	{ a.end() } -> std::same_as<typename C::iterator>;
	{ b.begin() } -> std::same_as<typename C::const_iterator>;
	{ b.end() } -> std::same_as<typename C::const_iterator>;
	{ a.cbegin() } -> std::same_as<typename C::const_iterator>;
	{ a.cend() } -> std::same_as<typename C::const_iterator>;
	{ a.size() } -> std::same_as<typename C::size_type>;
	{ a.max_size() } -> std::same_as<typename C::size_type>;
	{ a.empty() } -> std::convertible_to<bool>;
};
