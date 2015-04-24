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

    template <class String>
    static inline int get_code(const String &str, size_t index) {
        return base_to_code(str[index]);
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

    template <class String>
    static inline uint64_t get_index(const String &str, size_t pos, size_t num_index_chars) {
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

    // some older hardware does not have popcnt
    static inline bool has_popcnt() {
        #if defined(_MSC_VER) && defined(_M_AMD64)
            int result[4];
            __cpuidex(result, 0x00000001, 0x000000000);
            return (result[2] & (1 << 23)) != 0;
        #else
            return false;
        #endif
    }

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

    static inline int popcnt(uint64_t value, bool has_popcnt) {
        value = (value & 0x5555555555555555ull) + ((value >> 1) & 0x5555555555555555ull);
        value = (value & 0x3333333333333333ull) + ((value >> 2) & 0x3333333333333333ull);
        value = (value & 0x0f0f0f0f0f0f0f0full) + ((value >> 4) & 0x0f0f0f0f0f0f0f0full);
        value = (value & 0x00ff00ff00ff00ffull) + ((value >> 8) & 0x00ff00ff00ff00ffull);
        value = (value & 0x0000ffff0000ffffull) + ((value >> 16) & 0x0000ffff0000ffffull);
        return (int)(value + (value>>32));
    }

    /// variable length vector of packed bases
    class dna_string {
    public:
        typedef uint64_t word_type;
        static const size_t bases_per_value = sizeof(word_type) * 4;

    public:
        dna_string(size_t size=0) {
            num_bases = size;
            values.resize(num_bases * bases_per_value);
        }

        template<class InIter>
        dna_string(InIter b, InIter e) {
            num_bases = 0;
            append(b, e);
        }

        template<class String>
        dna_string(
            const String &str,
            size_t pos = 0, size_t n = ~(size_t)0
        ) {
            num_bases = 0;
            append(str.begin() + pos, str.begin() + std::min(n, str.size()));
        }

        template <class charT>
        dna_string(
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

        dna_string substr(size_t offset=0, size_t length=~(size_t)0, bool rev_comp=false) const {
            dna_string result;
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
                  word_type w = window((int64_t)(offset + length - (i+1) * bases_per_value));
                  result.values[i] = rev_comp_word(w);
              }
            }
            if (result.num_bases % bases_per_value) {
                result.values[nv-1] &= ~(word_type)0 << (((0-length) % bases_per_value) * 2);
            }
            return result;
        }

        size_t distance(const dna_string &rhs) const {
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

        bool operator==(const dna_string &rhs) const {
            return values == rhs.values && size() == rhs.size();
        }

        bool operator!=(const dna_string &rhs) const {
            return values != rhs.values || size() != rhs.size();
        }

        bool operator>(const dna_string &rhs) const {
            return values > rhs.values || values == rhs.values && size() > rhs.size();
        }

        bool operator<(const dna_string &rhs) const {
            return values < rhs.values || values == rhs.values && size() < rhs.size();
        }

        bool operator>=(const dna_string &rhs) const {
            return values > rhs.values || values == rhs.values && size() >= rhs.size();
        }

        bool operator<=(const dna_string &rhs) const {
            return values < rhs.values || values == rhs.values && size() <= rhs.size();
        }

        size_t find(const dna_string& str, size_t pos = 0, int max_distance = 0) const {
            size_t ssz = str.size();
            if (ssz == 0) {
                return pos;
            }

            size_t sz = size();
            if (pos + ssz > sz || values.size() == 0) {
                return std::string::npos;
            }

            const size_t bpv = bases_per_value;
            size_t nv = values.size() - 1;
            bool cpu_has_lzcnt = has_lzcnt();
            if (ssz >= 3 && pos < nv * bpv) {
                size_t max_base = ssz > size() - ssz;
                word_type r1c = 0x5555555555555555ull;
                word_type s0 = str.values[0];
                word_type rep0 = ~((s0 >> (bpv*2-2)) * r1c);
                word_type rep1 = ~(((s0 >> (bpv*2-4)) & 3) * r1c);
                word_type rep2 = ~(((s0 >> (bpv*2-6)) & 3) * r1c);
                for (size_t i = pos/bpv; i < nv; ++i) {
                    word_type v0 = values[i];
                    word_type v1 = values[i+1];
                    word_type mask = v0 ^ rep0;
                    mask &= ((v0 << 2 ) | (v1 >> (bpv*2-2))) ^ rep1;
                    mask &= ((v0 << 4 ) | (v1 >> (bpv*2-4))) ^ rep2;
                    mask &= mask << 1;
                    mask &= ~r1c;
                    size_t bit_pos = 0;
                    while (mask << bit_pos) {
                        int lz = (int)lzcnt(mask << bit_pos, cpu_has_lzcnt);
                        bit_pos += lz;
                        size_t search_pos = i * bpv + bit_pos/2;
                        ++bit_pos;
                        if (search_pos >= pos && compare(search_pos, ssz, str, max_distance) == 0) {
                            return search_pos;
                        }
                    }
                }

                pos = nv * bpv;
            }

            while (pos < sz) {
                if (compare(pos, ssz, str, max_distance) == 0) {
                    return pos;
                }
                ++pos;
            }
            return std::string::npos;
        }

        int compare(size_t pos, size_t ssz, const dna_string &str, int max_distance=0) const {
            pos = std::min(pos, num_bases);
            ssz = std::min(ssz, str.num_bases);
            ssz = std::min(ssz, num_bases - pos);

            const size_t bpv = bases_per_value;
            size_t nv = std::min(str.values.size(), (ssz+bpv-1)/bpv);
            for (size_t i = 0; i != nv; ++i) {
                word_type w = window((int64_t)pos);
                word_type s = str.values[i];
                if (s != w) {
                    if (i == ssz/bpv) {
                        s &= ~(word_type)0 << (((0-ssz) % bases_per_value) * 2);
                        w &= ~(word_type)0 << (((0-ssz) % bases_per_value) * 2);
                    }
                    return s == w ? 0 : s < w ? -1 : 1;
                }
                pos += bpv;
            }
            return 0;
        }

        int get_code(size_t index) const {
            size_t sh = ((bases_per_value - 1 - index) % bases_per_value) * 2;
            size_t off = index/bases_per_value;
            return index >= num_bases ? 0 : ((values[off] >> sh) & 0x03);
        }
    private:
        word_type window(int64_t base) const {
            int64_t max_size = (int64_t)values.size()*(int64_t)bases_per_value;
            word_type v0 = base >= 0 && base < max_size ? values[base/bases_per_value] : 0;
            if (base % bases_per_value == 0) {
                return v0;
            } else {
                base += bases_per_value;
                word_type v1 = base >= 0 && base < max_size ? values[base/bases_per_value] : 0;
                return v0 << ((base % bases_per_value)*2) | v1 >> (((0-base) % bases_per_value)*2);
            }
        }

        static word_type rev_comp_word(word_type x) {
            x = _byteswap_uint64(x);
            x = ((x & 0x0f0f0f0f0f0f0f0f) << 4) | ((x >> 4) & 0x0f0f0f0f0f0f0f0f);
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

    template <>
    static inline int get_code<class dna_string>(const dna_string &str, size_t index) {
        return str.get_code(index);
    }

    /// containter for bases ACGT and occasional runs of 'N' and other letters.
    class augmented_string : public dna_string {
        static const size_t lg_bases_per_index = 16;
        static const size_t bases_per_index = (size_t)1 << lg_bases_per_index;
        typedef uint32_t index_type;
        typedef uint32_t rle_type;
    public:
        augmented_string() {
        }

        augmented_string(const char *str) {
            const char *e = str;
            while (*e) ++e;
            append(str, e);
        }

        template<class String>
        augmented_string(
            const String &str,
            size_t pos = 0, size_t n = ~(size_t)0
        ) : dna_string() {
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
            return p > b && (p[-1] & 0xff) ? (p[-1] & 0xff) : (*(const dna_string*)this)[base];
        }

        void append(const char *str) {
            const char *e = str;
            while (*e) ++e;
            append(str, e);
        }

        void resize(size_t new_size, char chr='A') {
            size_t num_bases = size();
            dna_string::resize(new_size);
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

        std::string substr(size_t offset=0, size_t length=~(size_t)0, bool rev_comp=false) const {
            length = std::min(length, size() - offset);
            std::string result(length, ' ');
            if (!rev_comp) {
                for (size_t i = 0; i != length; ++i) {
                    result[i] = (*this)[offset+i];
                }
            } else {
                for (size_t i = 0; i != length; ++i) {
                    char chr = (*this)[offset + length - 1 - i];
                    result[i] = is_base(chr) ? code_to_base(3-base_to_code(chr)) : chr;
                }
            }
            return result;
        }
    private:

        template<class InIter>
        void append(InIter b, InIter e) {
            size_t num_bases = size();
            dna_string::append(b, e);
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
    operator>>(std::basic_istream<charT, traits>& is, dna_string& x) {
        std::string str;
        is >> str;
        x = dna_string(str);
        return is;
    }

    template <class charT, class traits>
    std::basic_ostream<charT, traits>&
    operator<<(std::basic_ostream<charT, traits>& os, const dna_string& x) {
        os << (std::string)x;
        return os;
    }

    template <class Type>
    Type rev_comp(const Type &x, typename Type::word_type *y = 0) {
      return x.substr(0, x.size(), true);
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

    /// two stage index, first index ordered by value, second by address.
    template <class String>
    class two_stage_index {
    public:
        typedef uint32_t index_type;
        typedef uint32_t addr_type;

        two_stage_index(
            String &string,
            size_t num_indexed_chars
        ) : string(string), num_indexed_chars(num_indexed_chars) {
            reindex();
        }

        void reindex() {
            if (
                num_indexed_chars <= 1 ||
                num_indexed_chars > 32 ||
                string.size() < num_indexed_chars ||
                (addr_type)string.size() != string.size()
            ) {
                throw std::out_of_range("reindex() failed");
            }

            size_t str_size = string.size();
            size_t index_size = (size_t)1 << (num_indexed_chars*2);

            index.resize(0);
            index.resize(index_size+1);
            addr.resize(str_size);

            uint64_t acc0 = 0;
            for (size_t i = 0; i != num_indexed_chars-1; ++i) {
                acc0 = acc0 * 4 + get_code(string, i);
            }

            uint64_t acc = acc0;
            for (size_t i = num_indexed_chars; i != str_size - num_indexed_chars; ++i) {
                acc = (acc * 4 + get_code(string, i)) & (index_size-1);
                index[acc]++;
            }

            addr_type cur = 0;
            for (size_t i = 0; i != index_size; ++i) {
                addr_type val = index[i];
                index[i] = cur;
                cur += val;
            }
            index[index_size] = cur;

            acc = acc0;
            for (size_t i = num_indexed_chars; i != str_size - num_indexed_chars; ++i) {
                acc = (acc * 4 + get_code(string, i)) & (index_size-1);
                addr[index[acc]++] = (addr_type)(i - num_indexed_chars + 1);
            }

            addr_type prev = 0;
            for (size_t i = 0; i != index_size; ++i) {
                std::swap(prev, index[i]);
            }
        }

        size_t find(const String& str, size_t pos = 0, int max_distance = 0) const {
            size_t num_seeds = str.size() / num_indexed_chars;
            active.resize(num_seeds);
            if (num_seeds == 0) {
                return std::string::npos;
            }

            // get seed values from string
            // touch the index in num_seeds places (with NTA hint)
            for (size_t i = 0; i != num_seeds; ++i) {
                active[i].idx = (index_type)get_index(str, i * num_indexed_chars, num_indexed_chars);
                touch_nta(addr.data() + index[i]);
            }

            // get seed values from string
            // touch the address vector in num_seeds places (with stream hint)
            for (size_t i = 0; i != num_seeds; ++i) {
                const addr_type *ptr = addr.data() + index[active[i].idx];
                const addr_type *end = addr.data() + index[active[i].idx+1];
                active[i].ptr = ptr;
                active[i].end = end;
                if (ptr != end) touch_stream(ptr);
            }

            for (size_t i = 0; i != num_seeds; ++i) {
                const addr_type *ptr = active[i].ptr;
                const addr_type *end = active[i].end;
                while (ptr != end && *ptr < pos) {
                   ++ptr;
                }
                active[i].ptr = ptr;
            }

            for (;;) {
              addr_type next = (addr_type)-1;
              for (size_t i = 0; i != num_seeds; ++i) {
                  std::cout << std::hex << active[i].idx << ":" << str.substr(i * num_indexed_chars, num_indexed_chars) << " ";
              }
              std::cout << std::endl;

              for (size_t i = 0; i != num_seeds; ++i) {
                  const addr_type *ptr = active[i].ptr;
                  active[i].start = (addr_type)-1;
                  if (ptr != active[i].end) {
                      addr_type pos = *ptr;
                      if (pos > i * num_indexed_chars) {
                          addr_type start = pos - (addr_type)(i * num_indexed_chars);
                          active[i].start = start;
                          next = std::min(next, start);
                      }
                  }
              }

              std::cout.width(10);
              for (size_t i = 1; i != num_seeds; ++i) {
                  addr_type start = active[i].start;
                  std::cout << (int)start << " ";
              }
              std::cout << "\n";

              if (next == (addr_type)-1) {
                  break;
              }

              size_t count = 0;
              for (size_t i = 1; i != num_seeds; ++i) {
                  addr_type start = active[i].start;
                  count += start == next ? 1 : 0;
              }

              if (count + max_distance >= num_seeds) {
                  // test some more
              }

              for (size_t i = 1; i != num_seeds; ++i) {
                  addr_type start = active[i].start;
                  if (start == next) {
                      active[i].ptr++;
                  }
              }
            }
            /*std::make_heap(active.begin(), active.end());
            for (;;) {
                active_state f = active.front();
                if (f.value == (addr_type)-1) {
                    break;
                }
                BOOST_MESSAGE(f.value);
                const addr_type *ptr = f.ptr + 1;
                f.value = ptr == f.end ? (addr_type)-1 : *ptr;
                std::pop_heap(active.begin(), active.end());
                active.back() = f;
                std::push_heap(active.begin(), active.end());
            }*/
            return std::string::npos;
        }

        template <class charT, class traits, class String>
        void dump(std::basic_ostream<charT, traits>& os) const {
            std::basic_ostream<charT, traits>::fmtflags save = os.flags();
            size_t str_size = string.size();
            size_t index_size = (size_t)1 << (num_indexed_chars*2);
            for (size_t i = 0; i != index_size; ++i) {
                for (index_type j = index[i]; j != index[i+1]; ++j) {
                    addr_type a = addr[j];
                    os << "[" << a << " " << string.substr(a, num_indexed_chars) << "]";
                }
                os << "\n";
            }
            os.setstate( save );
        }
    private:
        String &string;
        size_t num_indexed_chars;
        std::vector<index_type> index;
        std::vector<addr_type> addr;

        // search state
        struct active_state {
            const addr_type *ptr;
            const addr_type *end;
            index_type idx;
            addr_type start;

            bool operator<(const active_state &rhs) {
                return value > rhs.value;
            }
        };
        mutable std::vector<active_state> active;
    };

    template <class charT, class traits, class String>
    std::basic_ostream<charT, traits>&
    operator<<(std::basic_ostream<charT, traits>& os, const two_stage_index<String>& x) {
        x.dump<charT, traits, String>(os);
        return os;
    }

} }


