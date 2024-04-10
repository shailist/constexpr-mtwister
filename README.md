# constexpr-mtwister

Constexpr C++ implementation of the Mersenne Twister pseudo random number generation algorithm.

The implementation is a highly modified and cleaned up version of the libstdc++ implementation of `std::mersenne_twister_engine`.

## Usage

The implementation is complient with the standard's requirements, which can demonstrated as follows:

```cpp
#include "mtwister.hpp"

constexpr auto test32()
{
    mtwister::mt19937 gen32;
    gen32.discard(10000 - 1);
    return gen32();
}

constexpr auto test64()
{
    mtwister::mt19937_64 gen64;
    gen64.discard(10000 - 1);
    return gen64();
}

constexpr auto test32_result = test32();
constexpr auto test64_result = test64();

static_assert(test32_result == 4123659995);
static_assert(test64_result == 9981545732273789042ull);
```

This test is also included in the source code at the bottom of `mtwister.hpp`.
You can define `CONSTEXPR_MTWISTER_TEST=1` to run the test at compile time.

## License

As this code is heavily based on a part of libstdc++, it is licensed under GPLv3.
Refer to [LICENSE](LICENSE) for full information.
