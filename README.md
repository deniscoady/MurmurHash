# Murmur Hash TMP Implementation

I created this to use for compile-time hashing of strings. The implementation follows the algorithm originally found at https://en.wikipedia.org/wiki/MurmurHash


## Usage

    #include <iostream>
    #include "murmur.hpp"

    int main(int argc, char* argv[]) {
        // calculate hash of string "Hi" at compile time.
        std::uint32_t hi = MurmurHash<'H','i'>::value;
        std::cout << hi << std::endl;
    
        return 0;
    }
    
## License

Released under MIT license.
