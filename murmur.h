// MurmurHash Template Metaprogramming Implementation 
// https://github.com/dcoded/MurmurHash

/*
The MIT License (MIT)

Copyright (c) 2015 Denis Coady

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <cstdint>

/*

This implementation follows the following algorithm originally found at:
https://en.wikipedia.org/wiki/MurmurHash

Murmur3_32(key, len, seed)
    // Note: In this version, all integer arithmetic is performed with
    // unsigned 32 bit integers.
    // In the case of overflow, the result is constrained by the
    // application of modulo 2^{32} arithmetic.
    
    c1 ← 0xcc9e2d51
    c2 ← 0x1b873593
    r1 ← 15
    r2 ← 13
    m ← 5
    n ← 0xe6546b64
 
    hash ← seed

    for each fourByteChunk of key
        k ← fourByteChunk

        k ← k × c1
        k ← (k ROL r1)
        k ← k × c2

        hash ← hash XOR k
        hash ← (hash ROL r2)
        hash ← hash × m + n

    with any remainingBytesInKey
        remainingBytes ← SwapEndianOrderOf(remainingBytesInKey)
        // Note: Endian swapping is only necessary on big-endian machines.
        //       The purpose is to place the meaningful digits towards the low
        //       end of the value, so that these digits have the greatest
        //       potential to affect the low range digits in the subsequent
        //       multiplication.  Consider that locating the meaningful digits
        //       in the high range would produce a greater effect upon the high
        //       digits of the multiplication, and notably, that such high
        //       digits are likely to be discarded by the modulo arithmetic
        //       under overflow.  We don't want that.
        
        remainingBytes ← remainingBytes × c1
        remainingBytes ← (remainingBytes ROL r1)
        remainingBytes ← remainingBytes × c2

        hash ← hash XOR remainingBytes
 
    hash ← hash XOR len

    hash ← hash XOR (hash >> 16)
    hash ← hash × 0x85ebca6b
    hash ← hash XOR (hash >> 13)
    hash ← hash × 0xc2b2ae35
    hash ← hash XOR (hash >> 16)

*/

/* Define Magic Constants */
#define C1 0xcc9e2d51
#define C2 0x1b873593
#define R1 15
#define R2 13
#define M 5
#define N 0xe6546b64

// Forward declaration of implementation
template<std::uint32_t Hash, std::uint32_t Length, char... Args>
struct MurmurHash32Impl;

// Public interface function.
//
// Usage:
//      std::uint32_t x = MurmurHash<'H','i'>::value;
//
template<std::uint32_t Seed, char...Chars>
struct MurmurHash32 {
    enum { value = MurmurHash32Impl<Seed,0, Chars...>::value };
};

// Covers final hash calculations, no more characters
template <std::uint32_t Hash, std::uint32_t Length>
struct MurmurHash32Impl<Hash, Length> {
    // hash ^= len;
    enum { Hash1 = Hash ^ Length };
    // hash ^= (hash >> 16);
    enum { Hash2 = Hash1 ^ (Hash1 >> 16) };
    // hash *= 0x85ebca6b;
    enum { Hash3 = Hash2 * 0x85ebca6b };
    // hash ^= (hash >> 13);
    enum { Hash4 = Hash3 ^ (Hash3 >> 13) };
    // hash *= 0xc2b2ae35;
    enum { Hash5 = Hash4 * 0xc2b2ae35 };
    // hash ^= (hash >> 16);
    enum { Hash6 = Hash5 ^ (Hash5 >> 16) };

    enum { value = Hash6 };
};

// Calculates 1 remaining byte
template <std::uint32_t Hash, std::uint32_t Length, char A>
struct MurmurHash32Impl<Hash,Length,A> {
    // remainingBytes ← remainingBytes × c1
    enum { K0 = 0x00000000 ^ A };
    // remainingBytes ← remainingBytes × c1
    enum { K1 = K0 * C1 };
    // remainingBytes ← (remainingBytes ROL r1)
    enum { K2 = (K1 << R1) | (K1 >> (32-R1)) };
    // remainingBytes ← remainingBytes × c2
    enum { K3 = K2 * C2 };
    // hash ← hash XOR remainingBytes
    enum { Hash1 = Hash ^ K3 };
    // No more string, just final hash manipulations
    enum { value = MurmurHash32Impl<Hash1, Length+1>::value };
};

// Calculates 2 remaining bytes
template <std::uint32_t Hash, std::uint32_t Length, char A, char B>
struct MurmurHash32Impl<Hash,Length,A,B> {
    enum { K0 = 0x00000000  ^ (B << 8) };
    enum { K1 = K0 ^ A };
    // remainingBytes ← remainingBytes × c1
    enum { K2 = K1 * C1 };
    // remainingBytes ← (remainingBytes ROL r1)
    enum { K3 = (K2 << R1) | (K2 >> (32-R1)) };
    // remainingBytes ← remainingBytes × c2
    enum { K4 = K3 * C2 };
    // hash ← hash XOR remainingBytes
    enum { Hash1 = Hash ^ K4 };
    // No more string, just final hash manipulations
    enum { value = MurmurHash32Impl<Hash1, Length+2>::value };
};

// Calculates 3 remaining bytes
template <std::uint32_t Hash, std::uint32_t Length, char A, char B, char C>
struct MurmurHash32Impl<Hash,Length,A,B,C> {
    enum { K0 = 0  ^ (C << 16) };
    enum { K1 = K0 ^ (B << 8) };
    enum { K2 = K1 ^ A };
    // remainingBytes ← remainingBytes × c1
    enum { K3 = K2 * C1 };
    // remainingBytes ← (remainingBytes ROL r1)
    enum { K4 = (K3 << R1) | (K3 >> (32-R1)) };
    // remainingBytes ← remainingBytes × c2
    enum { K5 = K4 * C2 };
    // hash ← hash XOR remainingBytes
    enum { Hash1 = Hash ^ K5 };
    // No more string, just final hash manipulations
    enum { value = MurmurHash32Impl<Hash1, Length+3>::value };
};

// Calculates a full block (4 bytes)
template <std::uint32_t Hash, std::uint32_t Length, char A, char B, char C, char D, char...Tail>
struct MurmurHash32Impl<Hash,Length,A,B,C,D,Tail...> {
    // k ← fourByteChunk
    enum { K0 = (D << 24) | (C << 16) | (B << 8) | A };
    // k ← k × c1
    enum { K1 = K0 * C1 };
    // k ← (k ROL r1)
    enum { K2 = (K1 << R1) | (K1 >> (32-R1)) };
    // k ← k × c2
    enum { K3 = K2 * C2 };
    // hash ← hash XOR k
    enum { Hash1 = Hash ^ K3 };
    // hash ← (hash ROL r2)
    enum { Hash2 = (Hash1 << R2) | (Hash1 >> (32-R2)) };
    // hash ← hash × m + n
    enum { Hash3 = Hash2 * M + N };
    // Calculate rest of string and final hash manipulations
    enum { value = MurmurHash32Impl<Hash3, Length+4, Tail...>::value };
};