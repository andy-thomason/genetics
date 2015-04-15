#include <iostream>
#include <cstdlib>
#include "stdint.h"

#include <string>
#include <iosfwd>

namespace boost { namespace genetics {
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
        return "ACGT"[code & 3];
    }

    template <size_t N> class bases {
    public:
        static const size_t num_bases = N;
        typedef unsigned long long value_type;
        static const size_t bases_per_value = sizeof(value_type) * 4;
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

        /*template <class charT = char, class traits = std::char_traits<charT>,
            class Allocator = std::allocator<charT> >
        std::basic_string<charT, traits, Allocator>*/

        std::string
        to_string() const {
            //std::basic_string<charT, traits, Allocator> res(num_bases);
            std::string res;
            res.resize(num_bases);
            for (size_t i = 0; i != num_bases; ++i) {
                res[i] = "ACGT"[(*this)[i]];
            }
            return res;
        }
    private:
        template<class InIter>
        void init(InIter b, InIter e) {
            size_t i = 0;
            size_t d = 0;
            value_type acc = 0;
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

/*
    template <class Stream, int N> Stream &operator>>(Stream &s, bases<N> &b) {
        
    }

    template <class Stream, int N> Stream &operator<<(Stream &s, bases<N> &b) {
        
    }
*/
} }


