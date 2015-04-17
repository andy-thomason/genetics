#include <iostream>
#include <cstdlib>
#include "stdint.h"

#include <string>
#include <iosfwd>
#include <type_traits>

namespace boost { namespace genetics {
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

    static inline unsigned int code_to_base(int code) {
        return (('A' | 'C' << 8 | 'G' << 16 | 'T' << 24) >> (code * 8)) & 0xff;
    }

    /// fixed length vector of packed bases
    template <size_t N> class bases {
    public:
        typedef size_t is_bases;
        static const size_t num_bases = N;
        typedef unsigned long long word_type;
        static const size_t bases_per_value = sizeof(word_type) * 4;
        static const size_t num_values = (num_bases + bases_per_value-1) / bases_per_value;

    public:
        bases(unsigned long long val = 0) {
            values[0] = val;
            for (size_t i = 1; i < num_values; ++i) {
                values[i] = 0;
            }
        }

        template<class InIter>
        bases(InIter b, InIter e) {
            init(b, e);
        }

        template<class String>
        bases(
            const String &str,
            size_t pos = 0, size_t n = ~(size_t)0
        ) {
            init(str.begin() + pos, str.begin() + std::min(n, str.size()));
        }

        template <class charT>
        bases(
            const charT* str,
            size_t pos = 0, size_t n = ~(size_t)0
        ) {
            const charT *b = str + pos;
            const charT *e = str + pos;
            while (e - b != n && *e) ++e;
            init(b, e);
        }

        unsigned long long to_ullong() const {
            return values[0];
        }

        unsigned char operator[](size_t index) const {
            size_t sh = ((bases_per_value - 1 - index) % bases_per_value) * 2;
            size_t off = index/bases_per_value;
            return (unsigned char)((values[off] >> sh) & 0x03);
        }

        template<class String=std::string>
        String to_string() const {
            String res;
            res.resize(num_bases);
            for (size_t i = 0; i != num_bases; ++i) {
                res[i] = code_to_base((*this)[i]);
            }
            return res;
        }

        template <class Type>
        void rev_comp(Type &result, size_t offset=0, size_t length=~(size_t)0) const {
            length = std::min(length, result.size());
            length = std::min(length, size() - offset);
            size_t nv = (length + bases_per_value - 1) / bases_per_value;
            for (size_t i = 0; i < num_values; ++i) {
                word_type word = window((ptrdiff_t)(num_bases - bases_per_value - i * bases_per_value));
                result.values[i] = rev_comp_word(word);
            }
        }

        template <class Type>
        void substring(Type &result, size_t offset=0, size_t length=~(size_t)0) const {
            size_t nv = (length + bases_per_value - 1) / bases_per_value;
            for (size_t i = 0; i < nv; ++i) {
                result.values[i] = window((ptrdiff_t)(offset + i * bases_per_value));
            }
        }

        size_t distance(const bases &rhs) const {
            size_t diff = 0;
            for (size_t i = 0; i < num_values; ++i) {
                diff += count_word(values[i] ^ rhs.values[i]);
            }
            return diff;
        }

        size_t size() const {
          return num_bases;
        }
    private:
        template<class InIter>
        void init(InIter b, InIter e) {
            size_t i = 0;
            size_t d = 0;
            word_type acc = 0;
            for (; b != e && i != num_bases; ++i) {
                acc = acc * 4 + base_to_code((int)*b++);
                if ((i+1) % bases_per_value == 0) {
                    values[d++] = acc;
                    acc = 0;
                }
            }
            if (i % bases_per_value && d != num_values) {
                acc <<= (bases_per_value - i % bases_per_value) * 2;
                values[d++] = acc;
            }
            while (d < num_values) {
                values[d++] = 0;
            }
        }
    private:
        word_type window(ptrdiff_t base) const {
          word_type v0 = base < 0 ? 0 : values[base/bases_per_value];
          if (base % bases_per_value) {
            return v0;
          } else {
            word_type v1 = base/bases_per_value+1 < num_values ? values[base/bases_per_value+1] : 0;
            return v0 << (base % bases_per_value) | v1 >> ((0-base) % bases_per_value);
          }
        }

        static word_type rev_comp_word(word_type x) {
          x = _byteswap_uint64(x);
          x = ((x & 0x3333333333333333) << 2) | ((x >> 2) & 0x3333333333333333);
          return ~x;
        }

        static size_t count_word(word_type x) {
          x |= x >> 1;
          x &= 0x5555555555555555;
          return (size_t)__popcnt64(x);
        }

        unsigned long long values[num_values];
    };

    template <class charT, class traits, size_t N>
    std::basic_istream<charT, traits>&
    operator>>(std::basic_istream<charT, traits>& is, bases<N>& x) {
        std::string str;
        is >> str;
        x = bases<N>(str);
        return is;
    }

    template <class charT, class traits, size_t N>
    std::basic_ostream<charT, traits>&
    operator<<(std::basic_ostream<charT, traits>& os, const bases<N>& x) {
        os << x.to_string();
        return os;
    }

    template <class Type>
    Type rev_comp(const Type &x, typename Type::word_type *y = 0) {
      Type result;
      x.rev_comp(result);
      return result;
    }

    template <class Type>
    Type rev_comp(const Type &x, typename Type::iterator *b = 0) {
        Type result = x;
        std::reverse(result.begin(), result.end());
        for (Type::iterator i = result.begin(), e = result.end(); i != e; ++i) {
          int chr = *i;
          if (is_base(chr)) { *i = code_to_base(3-base_to_code(chr)); }
        }
        return result;
    }

    //template <class Type, Type::num_bases_type NB = Type::num_bases>
} }


