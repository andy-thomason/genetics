#include <iostream>
#include <cstdlib>

#include <string>
#include <iosfwd>
#include <type_traits>
#include <vector>
#include <algorithm>

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

    /// variable length vector of packed bases
    class bases {
    public:
        typedef uint64_t word_type;
        static const size_t bases_per_value = sizeof(word_type) * 4;

    public:
        bases(size_t size=0) {
            num_bases = size;
            values.resize(num_bases * bases_per_value);
        }

        template<class InIter>
        bases(InIter b, InIter e) {
            num_bases = 0;
            append(b, e);
        }

        template<class String>
        bases(
            const String &str,
            size_t pos = 0, size_t n = ~(size_t)0
        ) {
            num_bases = 0;
            append(str.begin() + pos, str.begin() + std::min(n, str.size()));
        }

        template <class charT>
        bases(
            const charT* str,
            size_t pos = 0, size_t n = ~(size_t)0
        ) {
            num_bases = 0;
            const charT *b = str + pos;
            const charT *e = str + pos;
            while (e - b != n && *e) ++e;
            append(b, e);
        }

        char operator[](size_t index) const {
            size_t sh = ((bases_per_value - 1 - index) % bases_per_value) * 2;
            size_t off = index/bases_per_value;
            return index >= num_bases ? 'N' : code_to_base((values[off] >> sh) & 0x03);
        }

        operator std::string() const {
            std::string res;
            res.resize(num_bases);
            for (size_t i = 0; i != num_bases; ++i) {
                res[i] = (*this)[i];
            }
            return res;
        }

        bases substr(size_t offset=0, size_t length=~(size_t)0, bool rev_comp=false) const {
            bases result;
            length = std::min(length, size() - offset);
            size_t nv = (length + bases_per_value - 1) / bases_per_value;
            result.num_bases = length;
            result.values.resize(nv);
            if (!rev_comp) {
              for (size_t i = 0; i < nv; ++i) {
                  result.values[i] = window((int64_t)(offset + i * bases_per_value));
              }
            } else {
              for (size_t i = 0; i < nv; ++i) {
                  result.values[i] = rev_comp_word(
                    window((int64_t)(offset + length - (i+1) * bases_per_value))
                  );
              }
            }
            if (result.num_bases % bases_per_value) {
                result.values[nv-1] &= ~(word_type)0 << (((0-length) % bases_per_value) * 2);
            }
            return result;
        }

        size_t distance(const bases &rhs) const {
            size_t diff = 0;
            for (size_t i = 0; i < values.size(); ++i) {
                diff += count_word(values[i] ^ rhs.values[i]);
            }
            return diff;
        }

        size_t size() const {
            return num_bases;
        }

        void reserve(size_t size) {
            values.reserve((size + bases_per_value - 1) / bases_per_value);
        }

        void resize(size_t size) {
            values.resize((size + bases_per_value - 1) / bases_per_value);
            if (size < num_bases && size % bases_per_value) {
                size_t last = size / bases_per_value;
                values[last] &= ~(word_type)0 << (((0-size) % bases_per_value) * 2);
            }
            num_bases = size;
        }

        void append(const char *str) {
            const char *e = str;
            while (*e) ++e;
            append(str, e);
        }

        template<class InIter>
        void append(InIter b, InIter e) {
            size_t max_bases = values.size() * bases_per_value;
            word_type acc = num_bases < max_bases ? values.back() >> (max_bases - num_bases) * 2 : 0;
            while (b != e) {
                acc = acc * 4 + base_to_code((int)*b++);
                num_bases++;
                if (num_bases % bases_per_value == 0) {
                    size_t index = (num_bases-1) / bases_per_value;
                    if (index >= values.size()) {
                        values.push_back(acc);
                    } else {
                        values[index] = acc;
                    }
                    acc = 0;
                }
            }

            if (num_bases % bases_per_value != 0) {
              acc <<= (max_bases - num_bases) * 2;
              size_t index = (num_bases-1) / bases_per_value;
              if (index >= values.size()) {
                  values.push_back(acc);
              } else {
                  values[index] = acc;
              }
            }
        }

        bool operator==(const bases &rhs) const {
            return values == rhs.values && size() == rhs.size();
        }

        bool operator!=(const bases &rhs) const {
            return values != rhs.values || size() != rhs.size();
        }

        bool operator>(const bases &rhs) const {
            return values > rhs.values || values == rhs.values && size() > rhs.size();
        }

        bool operator<(const bases &rhs) const {
            return values < rhs.values || values == rhs.values && size() < rhs.size();
        }

        bool operator>=(const bases &rhs) const {
            return values > rhs.values || values == rhs.values && size() >= rhs.size();
        }

        bool operator<=(const bases &rhs) const {
            return values < rhs.values || values == rhs.values && size() <= rhs.size();
        }

    private:
        word_type window(int64_t base) const {
            word_type v0 = base < 0 ? 0 : values[base/bases_per_value];
            if (base % bases_per_value == 0) {
                return v0;
            } else {
                word_type v1 = base/bases_per_value+1 < values.size() ? values[base/bases_per_value+1] : 0;
                return v0 << ((base % bases_per_value)*2) | v1 >> (((0-base) % bases_per_value)*2);
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

        size_t num_bases;
        std::vector<word_type> values;
    };

    /// containter for bases ACGT and occasional runs of 'N' and other letters.
    class bases_and_letters : public bases {
        static const size_t lg_bases_per_index = 16;
        static const size_t bases_per_index = (size_t)1 << lg_bases_per_index;
        typedef uint32_t index_type;
        typedef uint32_t rle_type;
    public:
        bases_and_letters() {
        }

        bases_and_letters(const char *str) {
            const char *e = str;
            while (*e) ++e;
            append(str, e);
        }

        template<class String>
        bases_and_letters(
            const String &str,
            size_t pos = 0, size_t n = ~(size_t)0
        ) : bases() {
            append(str.begin() + pos, str.begin() + std::min(n, str.size()));
        }

        char operator[](int64_t base) const {
            if ((uint64_t)base/bases_per_index > index.size()) return 'N';

            index_type i0 = index[(uint64_t)base/bases_per_index];

            index_type i1 =
              (uint64_t)base/bases_per_index+1 >= index.size() ?
              (index_type)rle.size() :
              index[(uint64_t)base/bases_per_index+1]
            ;

            rle_type search = (rle_type)((uint64_t)((base) % bases_per_index) << 8) | 0xff;
            const rle_type *b = rle.data() + i0;
            const rle_type *e = rle.data() + i1;
            const rle_type *p = std::lower_bound(b, e, search);
            return p > b && (p[-1] & 0xff) ? (p[-1] & 0xff) : (*(const bases*)this)[base];
        }

        void append(const char *str) {
            const char *e = str;
            while (*e) ++e;
            append(str, e);
        }

        void resize(size_t new_size, char chr='A') {
            size_t num_bases = size();
            bases::resize(new_size);
            if (new_size > num_bases) {
                while (num_bases < new_size) {
                    internal_append(num_bases, &chr, &chr+1);
                    num_bases++;
                }
            } else {
                index.resize(new_size / bases_per_index);
                rle.resize(index.empty() ? 0 : index.back());
            }
        }

        operator std::string() const {
            std::string res;
            res.resize(size());
            for (size_t i = 0; i != size(); ++i) {
                res[i] = (*this)[i];
            }
            return res;
        }

    private:

        template<class InIter>
        void append(InIter b, InIter e) {
            size_t num_bases = size();
            bases::append(b, e);
            internal_append(num_bases, b, e);
        }

        template<class InIter>
        void internal_append(size_t num_bases, InIter b, InIter e) {
            int prev_val = rle.empty() ? 0 : rle.back();
            while (b != e) {
                int chr = (int)*b++;
                int val = (is_base(chr)) ? 0 : chr;
                if (num_bases % bases_per_index == 0 || val != prev_val) {
                    if (num_bases % bases_per_index == 0) {
                        index.push_back((index_type)rle.size());
                    }
                    rle.push_back((rle_type)((num_bases % bases_per_index) * 256 + val));
                    prev_val = val;
                }
                num_bases++;
            }
        }

        std::vector<index_type> index;
        std::vector<rle_type> rle;
    };

    class annotation {
    public:
    private:
        typedef uint32_t index_type;
        std::vector<index_type> index;
        std::vector<uint8_t> data;
    };

    template <class charT, class traits>
    std::basic_istream<charT, traits>&
    operator>>(std::basic_istream<charT, traits>& is, bases& x) {
        std::string str;
        is >> str;
        x = bases(str);
        return is;
    }

    template <class charT, class traits>
    std::basic_ostream<charT, traits>&
    operator<<(std::basic_ostream<charT, traits>& os, const bases& x) {
        os << (std::string)x;
        return os;
    }

    template <class Type>
    Type rev_comp(const Type &x, typename Type::word_type *y = 0) {
      Type result;
      result.resize(x.size());
      x.substring(result, 0, x.size(), true);
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


