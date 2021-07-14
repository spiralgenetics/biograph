#pragma once

// Replace this with std::integer_sequence once we switch to C++14.

template<unsigned... Ints>
struct int_seq
{};

template<unsigned N, unsigned... Ints>
struct generate_int_seq : generate_int_seq<N-1, N-1, Ints...>
{};

template<unsigned... Ints>
struct generate_int_seq<0, Ints...> : int_seq<Ints...>
{};

