// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_UTILS_HPP
#define BOOST_GENETICS_UTILS_HPP

#include <iostream>
#include <cstdlib>

#include <string.h>
#include <string>
#include <iosfwd>
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <chrono>

#if !defined(_CRAYC) && !defined(__CUDACC__) && (!defined(__GNUC__) || (__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ > 3)))
    #if (defined(_M_IX86_FP) && (_M_IX86_FP >= 2)) || defined(__SSE2__)
        #include <mmintrin.h>
    #endif
#endif

#if defined(_MSC_VER) && defined(_M_AMD64)
    #define BOOST_GENETICS_IS_WIN64 1
#else
    #define BOOST_GENETICS_IS_WIN64 0
#endif

namespace boost { namespace genetics {
    typedef unsigned char uint8_t;
    typedef unsigned short uint16_t;
    typedef unsigned int uint32_t;
    typedef unsigned long long uint64_t;
    typedef long long int64_t;

    static inline bool is_base(unsigned int chr) {
        static const uint64_t vals = (
            (1 << ('A'-'A')) |
            (1 << ('C'-'A')) |
            (1 << ('G'-'A')) |
            (1 << ('T'-'A'))
        );
        unsigned int idx = (chr & ~0x20) - 'A';
        return idx <= 'T'-'A' && ((vals >> idx) & 1) != 0;
    }

    static inline int base_to_code(unsigned int chr) {
        static const uint64_t vals = (
            (1ull << ('C'-'C')*2) |
            (2ull << ('G'-'C')*2) |
            (3ull << ('T'-'C')*2)
        );
        unsigned int idx = (chr & ~0x20) - 'C';
        return idx <= 'T'-'C' ? (int)(vals >> (idx*2)) & 0x03 : 0x00;
    }
    
    static inline char code_to_base(int code) {
        return (char)((('A' | 'C' << 8 | 'G' << 16 | 'T' << 24) >> (code * 8)));
    }

    static inline bool is_whitespace(int chr) {
        return (unsigned)(chr & 0xff) <= ' ';
    }

    template <class Ptr>
    void touch_nta(Ptr ptr) {
        #if BOOST_GENETICS_IS_WIN64
            _mm_prefetch((const char *)ptr, _MM_HINT_NTA);
        #endif
    }

    template <class Ptr>
    void touch_stream(Ptr ptr) {
        #if BOOST_GENETICS_IS_WIN64
            _mm_prefetch((const char *)ptr, _MM_HINT_T1);
        #endif
    }

    template <class StringType>
    static inline uint64_t get_index(const StringType &str, size_t pos, size_t num_index_chars) {
        uint64_t result = 0;
        for (size_t i = 0; i != num_index_chars; ++i) {
            result = result * 4 + base_to_code(str[pos + i]);
        }
        return result;
    }

    // Some older hardware treats lzcnt like bsr.
    static inline bool has_lzcnt() {
        #if BOOST_GENETICS_IS_WIN64
            int result[4];
            __cpuidex(result, 0x80000001, 0x000000000);
            return (result[2] & (1 << 5)) != 0;
        #else
            return false;
        #endif
    }

    // Some older X86 hardware does not have popcnt.
    static inline bool has_popcnt() {
        #if BOOST_GENETICS_IS_WIN64
            int result[4];
            __cpuidex(result, 0x00000001, 0x000000000);
            return (result[2] & (1 << 23)) != 0;
        #elif defined(__GNUC__)
            /*unsigned result[4];
            __get_cpuid (0, result+0, result+1, result+2, result+3);
            return (result[2] & (1 << 23)) != 0;*/
            return true; // rash assumption!
        #else
            return false;
        #endif
    }

    static inline int soft_lzcnt(uint64_t value) {
        int result = 0;
        result = (value >> 32) ? result : result + 32;
        value = (value >> 32) ? (value >> 32) : value;
        result = (value >> 16) ? result : result + 16;
        value = (value >> 16) ? (value >> 16) : value;
        result = (value >> 8) ? result : result + 8;
        value = (value >> 8) ? (value >> 8) : value;
        static const std::uint8_t lzcnt_table[] = {
            8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
            2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        };
        return result + lzcnt_table[value];
    }

    // Leading zero count: either use machine instruction or C version.
    static inline int lzcnt(uint64_t value, bool has_lzcnt) {
        #if defined(_MSC_VER) && defined(_M_X64)
            return (int)__lzcnt64(value) ^ (has_lzcnt ? 0x00 : 0x1f);
        #elif defined(__GNUC__)
            // sadly the builtin is a pile of bovine excrement on both
            // clang and GCC
            //return (int)__builtin_clzll(value);
            if (has_lzcnt) {
                int64_t res;
                __asm__ volatile ("lzcnt %1, %0" : "=r"(res) : "r"(value));
                return (int)res;
            }
        #endif
        return soft_lzcnt(value);
    }

    static inline int soft_popcnt2(uint64_t value) {
        value = (value & 0x3333333333333333ull) + ((value >> 2) & 0x3333333333333333ull);
        value = (value & 0x0f0f0f0f0f0f0f0full) + ((value >> 4) & 0x0f0f0f0f0f0f0f0full);
        value = (value & 0x00ff00ff00ff00ffull) + ((value >> 8) & 0x00ff00ff00ff00ffull);
        value = (value & 0x0000ffff0000ffffull) + ((value >> 16) & 0x0000ffff0000ffffull);
        return (int)(value + (value>>32));
    }

    static inline int soft_popcnt(uint64_t value) {
        value = (value & 0x5555555555555555ull) + ((value >> 1) & 0x5555555555555555ull);
        return soft_popcnt2(value);
    }

    // Non-zero bit population count: either use machine instruction or C version.
    static inline int popcnt(uint64_t value, bool has_popcnt) {
        #if defined(_MSC_VER) && defined(_M_X64)
            if (has_popcnt) {
                return (int)__popcnt64(value);
            }
        #elif defined(__GNUC__)
            if (has_popcnt) {
                // sadly the builtin is a pile of bovine excrement on both
                // clang and GCC
                //return (int)__builtin_popcountll(value);
                int64_t res;
                __asm__ volatile ("popcnt %1, %0" : "=r"(res) : "r"(value));
                return (int)res;
            }
        #endif
        return soft_popcnt(value);
    }

    static inline uint64_t rev_comp_word(uint64_t x) {
        x = ((x & 0x00000000ffffffff) << 32) | ((x >> 32) & 0x00000000ffffffff);
        x = ((x & 0x0000ffff0000ffff) << 16) | ((x >> 16) & 0x0000ffff0000ffff);
        x = ((x & 0x00ff00ff00ff00ff) << 8)  | ((x >> 8) & 0x00ff00ff00ff00ff);
        x = ((x & 0x0f0f0f0f0f0f0f0f) << 4)  | ((x >> 4) & 0x0f0f0f0f0f0f0f0f);
        x = ((x & 0x3333333333333333) << 2)  | ((x >> 2) & 0x3333333333333333);
        return ~x;
    }

    static inline size_t count_word(uint64_t x, bool has_popcnt) {
        x |= x >> 1;  // count if either bit is 1.
        x &= 0x5555555555555555;
        return popcnt(x, has_popcnt);
    }

    template<class OutIter>
    OutIter make_int(OutIter &dest, uint64_t val) {
        static const uint64_t p10[] = {
            10000000000000000000ull, 1000000000000000000ull, 100000000000000000ull, 10000000000000000ull, 1000000000000000ull,
            100000000000000ull, 10000000000000ull, 1000000000000ull, 100000000000ull, 10000000000ull,
            1000000000ull, 100000000ull, 10000000ull, 1000000ull, 100000ull,
            10000ull, 1000ull, 100ull, 10ull, 1ull,
        };
        int start = val < p10[4] ? 4 : 0;
        start += val < p10[start + 8] ? 8 : 0;
        start += val < p10[start + 4] ? 4 : 0;
        start += val < p10[start + 2] ? 2 : 0;
        int lz = 19;
        for (int i = start; i != 20; ++i) {
            uint64_t one = p10[i];
            if (i >= lz || val >= one) {
                char c = '0';
                while (val >= one) {
                    val -= one;
                    c++;
                }
                *dest++ = c;
                lz = i;
            }
        }
        return dest;
    }

    // format SAM MD field.
    template<class OutIter, class String>
    OutIter make_MD_field(OutIter dest, String &s1, String &s2) {
        size_t len = s1.size();
        if (s2.size() != len || len == 0) { throw std::runtime_error("make_MD_field size mismatch"); }

        std::uint32_t matches = 0;
        for (size_t i = 0; i != len; ++i) {
            if (s1[i] == s2[i]) {
                ++matches;
            } else {
                dest = make_int(dest, matches);
                *dest++ = s1[i];
                matches = 0;
            }
        }

        if (matches) {
            dest = make_int(dest, matches);
        }

        return dest;
    }

    template<class OutIter, class String>
    OutIter make_rev_comp(OutIter dest, String &string) {
        auto i = string.end(), b = string.begin();
        for (; i != b; --i) {
            int chr = i[-1];
            if (is_base(chr)) {
                *dest++ = code_to_base(3-base_to_code(chr));
            } else {
                *dest++ = chr;
            }
        }
        return dest;
    }

    template<class OutIter>
    OutIter make_str(OutIter dest, const char *str) {
        while (*str) *dest++ = *str++;
        return dest;
    }

    template<class OutIter>
    OutIter make_rev_str(OutIter dest, const char *str) {
        const char *b = str;
        while (*str) str++;
        while (str != b) {
            *dest++ = *--str;
        }
        return dest;
    }

    // Consistent reverse complement interface:
    // A <-> T  C <-> G and reverse string
    // This is because DNA has two strands in opposite directions.

    // This one is for std::string and other strings
    template <class Type>
    Type rev_comp(const Type &x, typename Type::iterator *b = 0) {
        Type result = x;
        std::reverse(result.begin(), result.end());
        for (typename Type::iterator i = result.begin(), e = result.end(); i != e; ++i) {
          int chr = *i;
          if (is_base(chr)) { *i = code_to_base(3-base_to_code(chr)); }
        }
        return result;
    }

    /// Consistent interface to get a code (0..3) from a string.
    /// Note there is a specialisation for dna_string
    template <class StringType>
    static inline int get_code(const StringType &str, size_t index) {
        return base_to_code(str[index]);
    }

    const char *to_dna(uint64_t value, size_t len) {
        static char buf[sizeof(uint64_t)*4+1];
        size_t i = 0;
        for (; i != len; ++i) {
            if (i+1 == sizeof(buf)) break;
            buf[i] = code_to_base((value >> (len - i - 1)*2)&3);
        }
        buf[i] = 0;
        return buf;
    }

    class writer {
    public:
        writer(char *begin=nullptr, char *end=nullptr) :
            begin(begin), ptr(begin), end(end)
        {
        }
        
        template <class Type>
        void write(const Type *src, size_t size, size_t align) {
            char *aptr = begin + ((ptr - begin + align - 1) & (0-align));
            if (end != nullptr && aptr + size <= end) {
                if (aptr != ptr) memset(ptr, 0, aptr - ptr);
                memcpy(aptr, src, sizeof(Type) * size);
            }
            ptr = aptr + sizeof(Type) * size;
        }

        template <class A, class B> struct exists { typedef B type; };

        // Todo: in C++11 use alignof.
        template <class VecType>
        void write(const VecType &vec, size_t align = sizeof(typename VecType::value_type)) {
            write64(sizeof(typename VecType::value_type));
            write64(vec.size());
            write(vec.data(), vec.size(), align);
        }
        
        void write64(uint64_t value) {
            write(&value, 1, sizeof(value));
        }
        
        void write(const std::string &str) {
            write(str.c_str(), str.size(), 1);
        }
        
        bool is_end() const {
            return ptr == end;
        }
        
        char *get_ptr() const {
            return ptr;
        }

        size_t get_size() const {
            return (size_t)(ptr - begin);
        }
    private:
        char *begin;
        char *ptr;
        char *end;
    };
    
    class mapper {
    public:
        typedef void is_mapper;

        mapper(const char *begin, const char *end) :
            begin(begin), ptr(begin), end(end)
        {
        }
        
        uint64_t read64() {
            return *map<uint64_t>(1, sizeof(uint64_t));
        }
        
        template <class Type>
        Type *map(size_t size, size_t align) {
            const char *aptr = begin + (((ptr - begin) + align - 1) & (0-align));
            Type *res = (Type*)aptr;
            ptr = aptr + sizeof(Type) * size;
            if (ptr > end) throw(std::runtime_error("mapped file: overflow reading vector (check file format)"));
            return res;
        }
        
        bool is_end() const {
            return ptr == end;
        }

        const char *get_ptr() const {
            return ptr;
        }
    public:
        const char *begin;
        const char *ptr;
        const char *end;
    };
    
    
    /// Read only mapped vector.
    /// Todo: investigate using boost::interprocess containers.
    template<class Type>
    class mapped_vector {
    public:
        typedef Type value_type;
        
        mapped_vector() {
            sz = 0;
            dat = nullptr;
        }
        
        mapped_vector(mapper &map) {
            size_t value_size = (size_t)map.read64();
            if (value_size != sizeof(value_type)) {
                throw(std::runtime_error("mapped file: item size mismatch (check file format)"));
            }
            sz = (size_t)map.read64();
            dat = map.map<value_type>(sz, sizeof(value_type));
        }
        
        size_t size() const {
            return sz;
        }
        
        value_type &operator[](size_t idx) {
            return dat[idx];
        }
        
        const value_type &operator[](size_t idx) const {
            return dat[idx];
        }
        
        void resize(size_t new_size) {
            // todo: throw something
        }
        
        void push_back(const value_type &val) {
            // todo: throw something
        }
        
        const value_type *data() const {
            return dat;
        }

        value_type &back() {
            return dat[sz-1];
        }

        bool empty() const {
            return sz == 0;
        }
    private:
        size_t sz;
        value_type *dat;
    };

    struct chromosome {
        char name[80]; /// Note: these need to be fixed length strings for binary mapping.
        char info[80];
        size_t start;
        size_t end;
        size_t num_leading_N;
        size_t num_trailing_N;

        bool operator <(size_t pos) const {
            return end < pos;
        }

        chromosome() {
            memset(this, 0, sizeof(*this));
        }
    };

    //! Feedback from the algorithms.
    struct search_stats {
        size_t merges_done = 0;
        size_t compares_done = 0;
    };

    //! Parameters for inexact searches.
    struct search_params {
        //! max allowable errors
        size_t max_distance = 0;

        //! maximum gaps allowed (introns)
        size_t max_gap = 0;

        //! always do a linear scan of the indexed string
        bool always_brute_force = false;

        //! never do a linear scan of the indexed string
        bool never_brute_force = true;

        //! Limit of size of results returned.
        size_t max_results = 100;

        //! Search reverse complement strand also.
        bool search_rev_comp = true;
    };

    
    struct common_traits {
        typedef uint64_t DnaWordType;
    };

    //! \brief traits for unmapped classes (std::vector)
    struct unmapped_traits : common_traits {
        typedef std::vector<DnaWordType> DnaArrayType; 
        typedef std::vector<uint32_t> IndexArrayType;
        typedef std::vector<uint32_t> RleArrayType;
        typedef std::vector<uint32_t> TsiIndexArrayType;
        typedef std::vector<uint32_t> TsiAddrArrayType;
        typedef std::vector<uint32_t> SuffixArrayType;
        typedef std::vector<chromosome> FastaChromosomeType;
        typedef bool unmapped;
    };

    //! \brief traits for file mapped classes (mapped_vector)
    struct mapped_traits : common_traits {
        typedef mapped_vector<DnaWordType> DnaArrayType; 
        typedef mapped_vector<uint32_t> IndexArrayType;
        typedef mapped_vector<uint32_t> RleArrayType;
        typedef mapped_vector<uint32_t> TsiIndexArrayType;
        typedef mapped_vector<uint32_t> TsiAddrArrayType;
        typedef mapped_vector<uint32_t> SuffixArrayType;
        typedef mapped_vector<chromosome> FastaChromosomeType;
        typedef bool mapped;
    };
} }


#endif
