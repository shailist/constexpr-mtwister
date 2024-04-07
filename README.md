# constexpr-mtwister

Constexpr C++ implementation of the Mersenne Twister pseudo random number generation algorithm.

The implementation is pretty much a refactored version of https://github.com/ESultanik/mtwister, all credit for the implementation of the algorithm goes to them.

## Design Notes

Since `constexpr` (or rather `const`) variables cannot be modified, modifying the RNG's state is not possible. This means that instead of *advancing* the RNG's state, we need to *create* a new state instead.

Refer to the usage examples below to understand how to deal with this limitation.

## Usage

```cpp
#include "mtwister.hpp"

constexpr uint32_t seed = 0x1234;

// Create RNG state
constexpr auto rng = mtwister::Rng(seed);

// Generate a random number and the next rng state
constexpr auto random_result = rng.next();
constexpr uint32_t random_value = random_result.value;
constexpr auto& next_rng = random_result.rng;

// Usage within a constexpr function can make your life easier
constexpr uint32_t generate_random_number(uint32_t seed, size_t iterations)
{
    auto rng = mtwister::Rng(seed);
    uint32_t result = 0;

    for (size_t i = 0; i < iterations; ++i)
    {
        auto [next_rng, next_result] = rng.next();

        rng = next_rng;
        result = next_result;
    }

    return result;
}
```

## License

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org>
