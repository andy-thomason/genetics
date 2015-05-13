#ifndef BOOST_GENETICS_UTILS_HPP
#define BOOST_GENETICS_UTILS_HPP

#include <iostream>
#include <cstdlib>

#include <string>
#include <iosfwd>
#include <vector>
#include <algorithm>

#if !defined(_CRAYC) && !defined(__CUDACC__) && (!defined(__GNUC__) || (__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ > 3)))
    #if (defined(_M_IX86_FP) && (_M_IX86_FP >= 2)) || defined(__SSE2__)
        #include <mmintrin.h>
    #endif
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
        #if defined(_MSC_VER) && defined(_M_AMD64)
            _mm_prefetch((const char *)ptr, _MM_HINT_NTA);
        #endif
    }

    template <class Ptr>
    void touch_stream(Ptr ptr) {
        #if defined(_MSC_VER) && defined(_M_AMD64)
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

    // some older hardware treats lzcnt like bsr
    static inline bool has_lzcnt() {
        #if defined(_MSC_VER) && defined(_M_AMD64)
            int result[4];
            __cpuidex(result, 0x80000001, 0x000000000);
            return (result[2] & (1 << 5)) != 0;
        #else
            return false;
        #endif
    }

    // some older X86 hardware does not have popcnt
    static inline bool has_popcnt() {
        #if defined(_MSC_VER) && defined(_M_AMD64)
            int result[4];
            __cpuidex(result, 0x00000001, 0x000000000);
            return (result[2] & (1 << 23)) != 0;
        #else
            return false;
        #endif
    }

    // leading zero count: either use machine instruction or C version.
    static inline int lzcnt(uint64_t value, bool has_lzcnt) {
        #if defined(MSC_VER) && defined(_M_X64)
            return __lzcnt64(value) ^ (has_lzcnt ? 0x00 : 0x1f);
        #else
            int result = value ? 0 : 1;
            result = (value >> 32) ? result : result + 32;
            value = (value >> 32) ? (value >> 32) : value;
            result = (value >> 16) ? result : result + 16;
            value = (value >> 16) ? (value >> 16) : value;
            result = (value >> 8) ? result : result + 8;
            value = (value >> 8) ? (value >> 8) : value;
            result = (value >> 4) ? result : result + 4;
            value = (value >> 4) ? (value >> 4) : value;
            result = (value >> 2) ? result : result + 2;
            value = (value >> 2) ? (value >> 2) : value;
            result = (value >> 1) ? result : result + 1;
            value = (value >> 1) ? (value >> 1) : value;
            return result;
        #endif
    }

    // non-zero bit population count: either use machine instruction or C version.
    static inline int popcnt(uint64_t value, bool has_popcnt) {
        #if defined(MSC_VER) && defined(_M_X64)
            if (has_popcnt) {
                return (int)__popcnt64(value);
            }
        #endif
        value = (value & 0x5555555555555555ull) + ((value >> 1) & 0x5555555555555555ull);
        value = (value & 0x3333333333333333ull) + ((value >> 2) & 0x3333333333333333ull);
        value = (value & 0x0f0f0f0f0f0f0f0full) + ((value >> 4) & 0x0f0f0f0f0f0f0f0full);
        value = (value & 0x00ff00ff00ff00ffull) + ((value >> 8) & 0x00ff00ff00ff00ffull);
        value = (value & 0x0000ffff0000ffffull) + ((value >> 16) & 0x0000ffff0000ffffull);
        return (int)(value + (value>>32));
    }

    // consistent reverse complement interface
    // A <-> T  C <-> G and reverse string
    // this is because DNA has two strands in opposite directions.
    // this one is for dna_string and augmented_string
    template <class Type>
    Type rev_comp(const Type &x, typename Type::word_type *y = 0) {
      return x.substr(0, x.size(), true);
    }

    // this one is for std::string and other strings
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

    // return true if x[pos..] matches y to max_distance mismatches)
    template <class Lhs, class Rhs>
    bool match(const Lhs &x, size_t pos, const Rhs &y, size_t max_distance) {
        size_t n = y.size();
        size_t result = 0;
        //TCGAGACCATCCTGGCTAACACGGGGAAACCCCGTCTCCACTAAAAATACAAAAAGTTAG
        //TCGAGACCATCCTGGCTAACACGGGGAAACCCCGTCTCCACTAAAAATACAAAAAGTTAG
        //std::cout << x.substr(pos, y.size()) << "\n";
        //std::cout << y << "\n";
        for (size_t i = 0; i != n && result <= max_distance; ++i) {
            int xchr = x[pos + i];
            int ychr = y[i];
            result += ychr != xchr && ychr != 'N';
        }
        return result <= max_distance;
    }

    // consistent interface to get a code (0..3) from a string
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
        writer(char *begin, char *end) :
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

        template <class Type>
        void write(const std::vector<Type> &vec) {
            write64(vec.size());
            // todo: in C++11 use alignof
            write(vec.data(), vec.size(), sizeof(Type));
        }
        
        void write64(uint64_t value) {
            write(&value, 1, sizeof(value));
        }
        
        bool is_end() const {
            return ptr == end;
        }
        
        char *get_ptr() const {
            return ptr;
        }
    private:
        char *begin;
        char *ptr;
        char *end;
    };
    
    class mapper {
    public:
        mapper(const char *begin, const char *end) :
            begin(begin), ptr(begin), end(end)
        {
        }
        
        uint64_t read64() {
            return *map<uint64_t>(1, sizeof(uint64_t));
        }
        
        template <class Type>
        const Type *map(size_t size, size_t align) {
            const char *aptr = begin + (((ptr - begin) + align - 1) & (0-align));
            Type *res = (Type*)aptr;
            ptr = aptr + sizeof(Type) * size;
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
    
    
    // read only mapped vector
    template<class Type>
    class mapped_vector {
    public:
        typedef Type value_type;
        
        mapped_vector(mapper &map) {
            sz = (size_t)map.read64();
            dat = map.map<value_type>(sz, sizeof(value_type));
        }
        
        size_t size() const {
            return sz;
        }
        
        value_type operator[](size_t idx) const {
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
    private:
        size_t sz;
        const value_type *dat;
    };
} }


#endif