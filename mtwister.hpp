/**
 * @file mtwister.hpp
 * @author Shai List
 * @brief Constexpr C++ implementation of the Mersenne Twister pseudo random number generation algorithm.
 * @version 1.0
 * @date 2024-04-07
 * 
 * @copyright 
 *     This is free and unencumbered software released into the public domain.
 *     Anyone is free to copy, modify, publish, use, compile, sell, or
 *     distribute this software, either in source code form or as a compiled
 *     binary, for any purpose, commercial or non-commercial, and by any
 *     means.
 *     
 *     In jurisdictions that recognize copyright laws, the author or authors
 *     of this software dedicate any and all copyright interest in the
 *     software to the public domain. We make this dedication for the benefit
 *     of the public at large and to the detriment of our heirs and
 *     successors. We intend this dedication to be an overt act of
 *     relinquishment in perpetuity of all present and future rights to this
 *     software under copyright law.
 *     
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *     IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 *     OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 *     ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *     
 *     For more information, please refer to <https://unlicense.org>
 */

#pragma once

#include <array>

namespace mtwister
{
    /**
     * @brief The result of `Rng::next`. Contains the generated value and the next RNG state. See below for the contents of the struct.
     */
    struct RngResult;

    /**
     * @brief Mersenne Twister RNG state and methods.
     */
    class Rng
    {
    public:
        explicit constexpr Rng(uint32_t seed);

        template <typename = void>
        constexpr RngResult next() const;

    private:
        static constexpr uint32_t STATE_VECTOR_LENGTH = 624;

    private:
        explicit constexpr Rng(std::array<uint32_t, STATE_VECTOR_LENGTH> mt, int index);

        static constexpr Rng _seed_rand(uint32_t seed, std::array<uint32_t, STATE_VECTOR_LENGTH> old_mt, int old_index);

        std::array<uint32_t, STATE_VECTOR_LENGTH> mt;
        int index;
    };

    struct RngResult
    {
        Rng rng;
        uint32_t value;
    };
}

////////////////////
// Implementation //
////////////////////

namespace mtwister
{
    inline constexpr Rng::Rng(uint32_t seed) :
        Rng(_seed_rand(seed, {}, 0))
    {
    }

    template <typename>
    inline constexpr RngResult Rng::next() const
    {
        constexpr uint32_t UPPER_MASK = 0x80000000;
        constexpr uint32_t LOWER_MASK = 0x7fffffff;
        constexpr uint32_t TEMPERING_MASK_B = 0x9d2c5680;
        constexpr uint32_t TEMPERING_MASK_C = 0xefc60000;
        constexpr uint32_t STATE_VECTOR_M = 397;

        auto rng = *this;

        uint32_t y = 0;
        uint32_t mag[2] = { 0x0, 0x9908b0df }; // mag[x] = x * 0x9908b0df for x = 0,1
        
        if (rng.index >= STATE_VECTOR_LENGTH || rng.index < 0) {
            // Generate STATE_VECTOR_LENGTH words at a time.
            int kk = 0;
            if (rng.index >= STATE_VECTOR_LENGTH + 1 || rng.index < 0) {
                rng = _seed_rand(4357, rng.mt, rng.index);
            }

            for (kk = 0; kk < STATE_VECTOR_LENGTH - STATE_VECTOR_M; kk++) {
                y = (rng.mt[kk] & UPPER_MASK) | (rng.mt[kk + 1] & LOWER_MASK);
                rng.mt[kk] = rng.mt[kk + STATE_VECTOR_M] ^ (y >> 1) ^ mag[y & 0x1];
            }

            for (; kk < STATE_VECTOR_LENGTH - 1; kk++) {
                y = (rng.mt[kk] & UPPER_MASK) | (rng.mt[kk + 1] & LOWER_MASK);
                rng.mt[kk] = rng.mt[kk + (STATE_VECTOR_M - STATE_VECTOR_LENGTH)] ^ (y >> 1) ^ mag[y & 0x1];
            }
            
            y = (rng.mt[STATE_VECTOR_LENGTH - 1] & UPPER_MASK) | (rng.mt[0] & LOWER_MASK);
            
            rng.mt[STATE_VECTOR_LENGTH - 1] = rng.mt[STATE_VECTOR_M - 1] ^ (y >> 1) ^ mag[y & 0x1];
            rng.index = 0;
        }

        y = rng.mt[rng.index++];
        y ^= (y >> 11);
        y ^= (y << 7) & TEMPERING_MASK_B;
        y ^= (y << 15) & TEMPERING_MASK_C;
        y ^= (y >> 18);

        return RngResult{ rng, y };
    }

    inline constexpr Rng::Rng(std::array<uint32_t, STATE_VECTOR_LENGTH> mt, int index) :
        mt(mt),
        index(index)
    {
    }

    inline constexpr Rng Rng::_seed_rand(uint32_t seed, std::array<uint32_t, STATE_VECTOR_LENGTH> old_mt, int old_index)
    {
        std::array<uint32_t, STATE_VECTOR_LENGTH> new_mt = old_mt;
        int new_index = old_index;

        new_mt[0] = seed & 0xffffffff;
        for (new_index = 1; new_index < STATE_VECTOR_LENGTH; new_index++) {
            new_mt[new_index] = (6069 * new_mt[new_index - 1]) & 0xffffffff;
        }

        return Rng(new_mt, new_index);
    }
}
