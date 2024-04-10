/**
 * @file mtwister.hpp
 * @author Shai List
 * @brief Constexpr C++ implementation of the Mersenne Twister pseudo random number generation algorithm.
 * @version 1.0
 * @date 2024-04-10
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <type_traits>
#include <limits>

namespace mtwister
{
    namespace detail
    {
        template <typename UIntType, size_t w, bool = (w < static_cast<size_t>(std::numeric_limits<UIntType>::digits))>
        struct shift
        {
            static constexpr UIntType value = 0;
        };

        template <typename UIntType, size_t w>
        struct shift<UIntType, w, true>
        {
            static constexpr UIntType value = UIntType(1) << w;
        };

        template <typename UIntType, size_t w>
        constexpr UIntType shift_v = shift<UIntType, w>::value;

        template <typename T>
        using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

        template <typename SeedSeq>
        using seed_sequence_generate_t = decltype(std::declval<SeedSeq&>().generate(std::declval<uint_least32_t*>(), std::declval<uint_least32_t*>()));

        template<typename SeedSeq, typename Engine, typename Result, typename = seed_sequence_generate_t<SeedSeq>>
        using if_seed_sequence_for_t = std::enable_if_t<
            !std::is_same_v<remove_cvref_t<SeedSeq>, Engine> &&
            !std::is_unsigned_v<typename SeedSeq::result_type> &&
            !std::is_convertible_v<SeedSeq, Result>>;
    }

    /**
       * A generalized feedback shift register discrete random number generator.
       *
       * This algorithm avoids multiplication and division and is designed to be
       * friendly to a pipelined architecture.  If the parameters are chosen
       * correctly, this generator will produce numbers with a very long period and
       * fairly good apparent entropy, although still not cryptographically strong.
       *
       * The best way to use this generator is with the predefined mt19937 class.
       *
       * This algorithm was originally invented by Makoto Matsumoto and Takuji Nishimura.
       *
       * @tparam w  Word size, the number of bits in each element of the state vector.
       * @tparam n  The degree of recursion.
       * @tparam m  The period parameter.
       * @tparam r  The separation point bit index.
       * @tparam a  The last row of the twist matrix.
       * @tparam u  The first right-shift tempering matrix parameter.
       * @tparam d  The first right-shift tempering matrix mask.
       * @tparam s  The first left-shift tempering matrix parameter.
       * @tparam b  The first left-shift tempering matrix mask.
       * @tparam t  The second left-shift tempering matrix parameter.
       * @tparam c  The second left-shift tempering matrix mask.
       * @tparam l  The second right-shift tempering matrix parameter.
       * @tparam f  Initialization multiplier.
       */
    template <
        typename UIntType, size_t w, size_t n, size_t m, size_t r,
        UIntType a, size_t u, UIntType d, size_t s,
        UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    class mersenne_twister_engine
    {
        static_assert(std::is_unsigned_v<UIntType>, "result_type must be an unsigned integral type");
        static_assert(1u <= m && m <= n, "template argument substituting m out of bounds");
        static_assert(r <= w, "template argument substituting r out of bound");
        static_assert(u <= w, "template argument substituting u out of bound");
        static_assert(s <= w, "template argument substituting s out of bound");
        static_assert(t <= w, "template argument substituting t out of bound");
        static_assert(l <= w, "template argument substituting l out of bound");
        static_assert(w <= std::numeric_limits<UIntType>::digits, "template argument substituting w out of bound");
        static_assert(a <= (detail::template shift_v<UIntType, w> - 1), "template argument substituting a out of bound");
        static_assert(b <= (detail::template shift_v<UIntType, w> - 1), "template argument substituting b out of bound");
        static_assert(c <= (detail::template shift_v<UIntType, w> - 1), "template argument substituting c out of bound");
        static_assert(d <= (detail::template shift_v<UIntType, w> - 1), "template argument substituting d out of bound");
        static_assert(f <= (detail::template shift_v<UIntType, w> - 1), "template argument substituting f out of bound");

        template<typename SeedSeq>
        using if_seed_sequence_t = detail::if_seed_sequence_for_t<SeedSeq, mersenne_twister_engine, UIntType>;

    public:
        // The type of the generated random value.
        using result_type = UIntType;

        // parameter values
        static constexpr size_t      word_size = w;
        static constexpr size_t      state_size = n;
        static constexpr size_t      shift_size = m;
        static constexpr size_t      mask_bits = r;
        static constexpr result_type xor_mask = a;
        static constexpr size_t      tempering_u = u;
        static constexpr result_type tempering_d = d;
        static constexpr size_t      tempering_s = s;
        static constexpr result_type tempering_b = b;
        static constexpr size_t      tempering_t = t;
        static constexpr result_type tempering_c = c;
        static constexpr size_t      tempering_l = l;
        static constexpr result_type initialization_multiplier = f;
        static constexpr result_type default_seed = 5489u;

        // constructors and member functions

        constexpr mersenne_twister_engine();

        explicit constexpr mersenne_twister_engine(result_type value);

        /**
         * @brief Constructs a mersenne_twister_engine random number generator
         *        engine seeded from the seed sequence @p q.
         *
         * @param q the seed sequence.
         */
        template<typename SeedSeq, typename = if_seed_sequence_t<SeedSeq>>
        explicit constexpr mersenne_twister_engine(SeedSeq& seq);

        constexpr void seed(result_type value = default_seed);

        template<typename SeedSeq>
        constexpr if_seed_sequence_t<SeedSeq> seed(SeedSeq& seq);

        /**
         * @brief Gets the smallest possible value in the output range.
         */
        static constexpr result_type min();

        /**
         * @brief Gets the largest possible value in the output range.
         */
        static constexpr result_type max();

        /**
         * @brief Discard a sequence of random numbers.
         */
        constexpr void discard(unsigned long long z);

        constexpr result_type operator()();

    private:
        constexpr void _advance_state();

        UIntType m_state[state_size];
        size_t m_index;
    };
}

////////////////////
// Implementation //
////////////////////

namespace mtwister
{
    namespace detail
    {
        template<int s, int which = ((s <= 8 * sizeof(int)) + (s <= 8 * sizeof(long)) + (s <= 8 * sizeof(long long)))>
        struct select_uint_least_t
        {
            static_assert(which < 0, "sorry, would be too much trouble for a slow result");
        };

        template<int s>
        struct select_uint_least_t<s, 4>
        {
            using type = unsigned int;
        };

        template<int s>
        struct select_uint_least_t<s, 3>
        {
            using type = unsigned long;
        };

        template<int s>
        struct select_uint_least_t<s, 2>
        {
            using type = unsigned long long;
        };

        template<typename T>
        inline constexpr T lg(T n)
        {
            return (sizeof(+n) * 8 - 1) - (sizeof(+n) == sizeof(long long) ? __builtin_clzll(+n)
                : (sizeof(+n) == sizeof(long)
                    ? __builtin_clzl(+n)
                    : __builtin_clz(+n)));
        }

        template<typename T, T m, T a, T c, bool big_enough = (!(m& (m - 1)) || (T(-1) - c) / a >= m - 1), bool schrage_ok = m % a < m / a>
        struct Mod
        {
            static constexpr T calc(T x)
            {
                using T2 = typename select_uint_least_t<lg(a) + lg(m) + 2>::type;
                return static_cast<T>((static_cast<T2>(a) * x + c) % m);
            }
        };

        template <typename T, T m, T a, T c>
        struct Mod<T, m, a, c, false, true>
        {
            static constexpr T calc(T x)
            {
                if (a == 1)
                    x %= m;
                else
                {
                    static constexpr T q = m / a;
                    static constexpr T r = m % a;

                    T t1 = a * (x % q);
                    T t2 = r * (x / q);

                    if (t1 >= t2)
                    {
                        x = t1 - t2;
                    }
                    else
                    {
                        x = m - t2 + t1;
                    }
                }

                if (c != 0)
                {
                    const T d = m - x;

                    if (d > c)
                    {
                        x += c;
                    }
                    else
                    {
                        x = c - d;
                    }
                }
                return x;
            }
        };

        template<typename T, T m, T a, T c, bool s>
        struct Mod<T, m, a, c, true, s>
        {
            static constexpr T calc(T x)
            {
                T res = a * x + c;

                if constexpr (m != T{})
                {
                    res %= m;
                }

                return res;
            }
        };

        template<typename T, T m, T a = 1, T c = 0>
        inline constexpr T mod(T x)
        {
            if constexpr (a == 0)
            {
                return c;
            }
            else
            {
                return Mod<T, m, a, c>::calc(x);
            }
        }
    }

    template<typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    inline constexpr mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::mersenne_twister_engine() :
        mersenne_twister_engine(default_seed)
    {}

    template<typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    inline constexpr mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::mersenne_twister_engine(result_type value) :
        m_state(),
        m_index(0)
    {
        seed(value);
    }

    template<typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    template<typename SeedSeq, typename>
    inline constexpr mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::mersenne_twister_engine(SeedSeq& seq) :
        m_state(),
        m_index(0)
    {
        seed(seq);
    }

    template<typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    inline constexpr void mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::seed(result_type value)
    {
        m_state[0] = detail::mod<UIntType, detail::shift_v<UIntType, w>>(value);

        for (UIntType i = 1; i < state_size; ++i)
        {
            UIntType x = m_state[i - 1];
            x ^= x >> (w - 2);
            x *= f;
            x += detail::mod<UIntType, n>(i);
            m_state[i] = detail::mod<UIntType, detail::shift_v<UIntType, w>>(x);
        }

        m_index = state_size;
    }

    template<typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    template<typename SeedSeq>
    inline constexpr mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::if_seed_sequence_t<SeedSeq> mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::seed(SeedSeq& seq)
    {
        static constexpr UIntType upper_mask = ~((UIntType{ 1 } << r) - 1);
        static constexpr size_t k = (w + 31) / 32;

        uint_least32_t arr[n * k];
        seq.generate(arr + 0, arr + n * k);

        bool zero = true;
        for (size_t i = 0; i < state_size; ++i)
        {
            UIntType factor = 1u;
            UIntType sum = 0u;

            for (size_t j = 0; j < k; ++j)
            {
                sum += arr[k * i + j] * factor;
                factor *= detail::shift_v<UIntType, 32>;
            }

            m_state[i] = detail::mod<UIntType, detail::shift_v<UIntType, w>>(sum);

            if (zero)
            {
                if (i == 0)
                {
                    if ((m_state[0] & upper_mask) != 0u)
                    {
                        zero = false;
                    }
                }
                else if (m_state[i] != 0u)
                {
                    zero = false;
                }
            }
        }

        if (zero)
        {
            m_state[0] = detail::shift_v<UIntType, w - 1>;
        }

        m_index = state_size;

    }

    template<typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    inline constexpr typename mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::result_type mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::min()
    {
        return 0;
    }

    template<typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    inline constexpr typename mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::result_type mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::max()
    {
        return detail::shift_v<UIntType, w>::value - 1;
    }

    template<typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    inline constexpr void mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::discard(unsigned long long z)
    {
        while (z > state_size - m_index)
        {
            z -= state_size - m_index;
            _advance_state();
        }

        m_index += z;
    }

    template<typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    inline constexpr typename mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::result_type mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::operator()()
    {
        // Reload the vector - cost is O(n) amortized over n calls.
        if (m_index >= state_size)
            _advance_state();

        // Calculate o(x(i)).
        result_type z = m_state[m_index % state_size];
        ++m_index;

        z ^= (z >> u) & d;
        z ^= (z << s) & b;
        z ^= (z << t) & c;
        z ^= (z >> l);

        return z;
    }

    template<typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f>
    inline constexpr void mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f>::_advance_state()
    {
        constexpr UIntType upper_mask = ~((UIntType{ 1 } << r) - 1);
        constexpr UIntType lower_mask = ~upper_mask;

        for (size_t k = 0; k < (n - m); ++k)
        {
            UIntType y = ((m_state[k] & upper_mask) | (m_state[k + 1] & lower_mask));
            m_state[k] = (m_state[k + m] ^ (y >> 1) ^ ((y & 0x01) ? a : 0));
        }

        for (size_t k = (n - m); k < (n - 1); ++k)
        {
            UIntType y = ((m_state[k] & upper_mask) | (m_state[k + 1] & lower_mask));
            m_state[k] = (m_state[k + (m - n)] ^ (y >> 1) ^ ((y & 0x01) ? a : 0));
        }

        UIntType y = ((m_state[n - 1] & upper_mask) | (m_state[0] & lower_mask));
        m_state[n - 1] = (m_state[m - 1] ^ (y >> 1) ^ ((y & 0x01) ? a : 0));
        m_index = 0;
    }

    using mt19937 = mersenne_twister_engine<unsigned int, 32, 624, 397, 31, 0x9908b0df, 11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18, 1812433253>;
    using mt19937_64 = mersenne_twister_engine<unsigned long long, 64, 312, 156, 31, 0xb5026f5aa96619e9ULL, 29, 0x5555555555555555ULL, 17, 0x71d67fffeda60000ULL, 37, 0xfff7eee000000000ULL, 43, 6364136223846793005ULL>;
}

///////////
// Tests //
///////////

/**
 * @brief Set to 1 to include the tests in your build process.
 */
#ifndef CONSTEXPR_MTWISTER_TEST
#define CONSTEXPR_MTWISTER_TEST (0)
#endif

namespace mtwister
{
#if CONSTEXPR_MTWISTER_TEST
    namespace detail
    {
        namespace test
        {
            constexpr auto test32()
            {
                mt19937 gen32;
                gen32.discard(10000 - 1);
                return gen32();
            }

            constexpr auto test64()
            {
                mt19937_64 gen64;
                gen64.discard(10000 - 1);
                return gen64();
            }

            constexpr auto test32_result = test32();
            constexpr auto test64_result = test64();

            static_assert(test32_result == 4123659995);
            static_assert(test64_result == 9981545732273789042ull);
        }
    }
#endif
}
