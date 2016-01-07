// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//! \file
//! \brief DNA string Container class definitions.


#ifndef BOOST_GENETICS_DNA_STRING_HPP
#define BOOST_GENETICS_DNA_STRING_HPP

#include <boost/genetics/utils.hpp>

#include <string>

namespace boost { namespace genetics {
    //! \brief This class stores DNA strings compactly allowing 32 or more bases to
    //! be accessed in a single instruction.

    //! Like many of the container classes in this library it can be specialised
    //! into a standard (`std::vector`) version for construction and a read-only
    //! mapped (`mapped_vector`) version for high performance use.

    //! \tparam WordType Integer word type, typically 64-bit unsigned integer `uint64_t`.
    //! \tparam ArrayType Container array type, typically `std::vector<uint64_t`.
    template<class Traits>
    class basic_dna_string {
    public:
        typedef typename Traits::DnaWordType word_type;
        typedef typename Traits::DnaArrayType array_type;
        static const size_t bases_per_value = sizeof(word_type) * 4;
        static const size_t npos = (size_t)-1;
        typedef basic_dna_string<Traits> this_type;

        typedef typename Traits::SuffixArrayType addr_array_type;
        typedef typename Traits::SuffixArrayType::value_type addr_type;

        struct bwt_sort_type { addr_type group, next_group, addr; };
        struct ibwt_sort_type { addr_type addr, chr; };

        static_assert(
            sizeof(word_type) >= sizeof(addr_type),
            "suffix array address type must be smaller than the string word type"
        );

    public:
        //! \brief Default constructor.
        basic_dna_string() {
            num_bases = 0;
        }

        //! \brief Construct and empty dna_string with size elements (all 'A')
        basic_dna_string(size_t size) {
            num_bases = size;
            values.resize((num_bases+bases_per_value-1)/bases_per_value);
        }

        //! \brief Construct a dna_string from a range of memory.
        template<class InIter>
        basic_dna_string(InIter b, InIter e, bool randomise=false) {
            num_bases = 0;
            append(b, e, randomise);
        }

        //! \brief Construct a dna_string from a C++ string
        template<class StrChar, class StrTraits, class StrAllocator>
        basic_dna_string(
            const std::basic_string<StrChar, StrTraits, StrAllocator> &str,
            size_t pos = 0, size_t n = ~(size_t)0, bool randomise = false
        ) {
            num_bases = 0;
            append(str.data() + pos, str.data() + std::min(n, str.size()), randomise);
        }

        //! \brief Construct a dna_string from a substring
        template <class charT>
        basic_dna_string(
            const charT* str,
            size_t pos = 0, size_t n = ~(size_t)0, bool randomise=false
        ) {
            num_bases = 0;
            const charT *b = str + pos;
            const charT *e = str + pos;
            while ((size_t)(e - b) != n && *e) ++e;
            append(b, e, randomise);
        }

        //! \brief Construct a dna_string from a mapper object (mapped_dna_string only).
        template <class Mapper>
        basic_dna_string(Mapper &map, typename Mapper::is_mapper *p=0) :
            num_bases((size_t)map.read64()),
            values(map)
        {
        }

        //! \brief Read a single base (in ASCII) from a dna_string
        char operator[](size_t index) const {
            size_t sh = ((bases_per_value - 1 - index) % bases_per_value) * 2;
            size_t off = index / bases_per_value;
            return index >= num_bases ? 'N' : code_to_base((values[off] >> sh) & 0x03);
        }

        //! \brief Convert to a C++ string.
        operator std::string() const {
            std::string res;
            res.resize(num_bases);
            for (size_t i = 0; i != num_bases; ++i) {
                res[i] = (*this)[i];
            }
            return res;
        }

        //! \brief Get a substring from a dna_string. The rev_comp parameter allows
        //! this to be the reverse complement of the substring.
        basic_dna_string substr(size_t offset=0, size_t length=~(size_t)0, bool rev_comp=false) const {
            basic_dna_string result;
            length = std::min(length, size() - offset);
            size_t nv = (length + bases_per_value - 1) / bases_per_value;
            result.num_bases = length;
            result.values.resize(nv);
            if (!rev_comp) {
              for (size_t i = 0; i < nv; ++i) {
                  size_t addr = offset + i * bases_per_value;
                  result.values[i] = window(addr);
              }
            } else {
              for (size_t i = 0; i < nv; ++i) {
                  size_t addr = offset + length - (i+1) * bases_per_value;
                  word_type w = window(addr);
                  result.values[i] = rev_comp_word(w);
              }
            }
            if (result.num_bases % bases_per_value) {
                result.values[nv-1] &= ~(word_type)0 << (((0-length) % bases_per_value) * 2);
            }
            return result;
        }

        //! \brief Return the number of bases in this string.
        size_t size() const {
            return num_bases;
        }

        //! \brief Reserve extra bases in the array for appending.
        //! This makes appending a little faster.
        void reserve(size_t size) {
            values.reserve((size + bases_per_value - 1) / bases_per_value);
        }

        //! \brief Resize the string, appending 'A' to the end if expanding.
        void resize(size_t size) {
            values.resize((size + bases_per_value - 1) / bases_per_value);
            if (size < num_bases && size % bases_per_value) {
                size_t last = size / bases_per_value;
                values[last] &= ~(word_type)0 << (((0-size) % bases_per_value) * 2);
            }
            num_bases = size;
        }

        //! \brief Append a C string.
        void append(const char *str, bool randomise=false) {
            const char *e = str;
            while (*e) ++e;
            append(str, e, randomise);
        }

        //! \brief Append ascii characters (A, C, G, T) to the string.
        //! Whitespace will be ignored and other characters (eg. N) may be mapped to a random value.
        template<class InIter>
        void append(InIter b, InIter e, bool randomise) {
            size_t max_bases = values.size() * bases_per_value;
            word_type acc = num_bases < max_bases ? values.back() >> (max_bases - num_bases) * 2 : 0;
            // todo: check the repeat of these numbers (and spectral purity if pedantry demands!)
            std::int32_t random_acc = 0xd3ac6435;
            const std::int32_t random_xor = 0xa9831bc5;
            while (b != e) {
                int chr = (int)*b++;
                if (!is_whitespace(chr)) {
                    acc = acc * 4 + base_to_code(chr);
                    if (randomise && !is_base(chr)) {
                        // acme branch-free rotate and xor random number generator
                        random_acc = ((random_acc >> 31) & random_xor ) ^ (random_acc << 1);
                        acc |= random_acc & 3;
                    }
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
            }

            if (num_bases % bases_per_value != 0) {
              //printf("%016llx sh=%d\n", (long long)acc, (int)((0-num_bases)%bases_per_value * 2));
              acc <<= (0-num_bases) % bases_per_value * 2;
              //printf("%016llx..\n", (long long)acc);
              size_t index = (num_bases-1) / bases_per_value;
              if (index >= values.size()) {
                  values.push_back(acc);
              } else {
                  values[index] = acc;
              }
            }
        }

        //! \brief Comparison operator
        bool operator==(const basic_dna_string &rhs) const {
            return values == rhs.values && size() == rhs.size();
        }

        //! \brief Comparison operator
        bool operator!=(const basic_dna_string &rhs) const {
            return values != rhs.values || size() != rhs.size();
        }

        //! \brief Comparison operator
        bool operator>(const basic_dna_string &rhs) const {
            return values > rhs.values || (values == rhs.values && size() > rhs.size());
        }

        //! \brief Comparison operator
        bool operator<(const basic_dna_string &rhs) const {
            return values < rhs.values || (values == rhs.values && size() < rhs.size());
        }

        //! \brief Comparison operator
        bool operator>=(const basic_dna_string &rhs) const {
            return values > rhs.values || (values == rhs.values && size() >= rhs.size());
        }

        //! \brief Comparison operator
        bool operator<=(const basic_dna_string &rhs) const {
            return values < rhs.values || (values == rhs.values && size() <= rhs.size());
        }

        //! \brief Brute force string search. For a more refined aproach, use two_stage_index.
        //! This method performs a linear search on DNA data taking O(max_bases) time.
        //! \tparam String DNA sequence string type, typically default dna_string.
        //! \param search_str DNA string to search.
        //! \param start_pos Zero-based offset to start the search.
        //! \param max_bases maxiumum number of bases to search.
        size_t find(const this_type& str, size_t start_pos = 0, size_t max_bases = ~(size_t)0) const {
            return find_inexact(str, start_pos, max_bases, 0);
        }


        //! \brief Brute force string search. For a more refined aproach, use two_stage_index.
        //! This method performs a linear search on DNA data taking O(max_bases) time.
        //! \tparam String DNA sequence string type, typically default dna_string.
        //! \param search_str DNA string to search.
        //! \param start_pos Zero-based offset to start the search.
        //! \param max_bases maxiumum number of bases to search.
        //! \param max_distance number of allowable errors in the search.
        size_t find_inexact(
            const std::string &search_str,
            size_t start_pos = 0,
            size_t max_bases = ~(size_t)0,
            size_t max_distance = 0
        ) const {
            basic_dna_string<unmapped_traits> dna_str(search_str);
            size_t pos = start_pos;
            size_t ssz = search_str.size();
            if (ssz == 0) {
                return pos;
            }

            size_t sz = size();
            if (pos >= sz) {
                return basic_dna_string::npos;
            }

            max_bases = std::min((size_t)(sz - pos), max_bases);
            size_t last = pos + max_bases;
            if (pos + ssz > last || values.size() == 0) {
                return basic_dna_string::npos;
            }

            bool cpu_has_popcnt = has_popcnt();
            bool cpu_has_lzcnt = has_lzcnt();
            //printf("%d %d\n", cpu_has_popcnt, cpu_has_lzcnt);

            // Come gentle pedantry, shine upon this code.
            const size_t bpv = bases_per_value;
            size_t nv = std::min(values.size() - 1, last/bpv);
            const word_type *search_values = dna_str.get_values().data();
            word_type s0 = search_values[0];
            word_type s0mask = ssz >= bpv ? ~(word_type)0 : ~(word_type)0 << (bpv*2-ssz*2);
            if (ssz >= 4 && max_distance == 0 && pos < nv * bpv) {
                // Search for four characters, bpv at a time.
                // This should be about 8-32 times faster than strstr.
                // Note: these strings can be several GB and have billions of
                // operations which makes this an essential optimisation.

                // The variables rep0-3 contain inverse duplicates of the first four characters.
                word_type r1c = 0x5555555555555555ull;
                word_type rep0 = ~((s0 >> (bpv*2-2)) * r1c);
                word_type rep1 = ~(((s0 >> (bpv*2-4)) & 3) * r1c);
                word_type rep2 = ~(((s0 >> (bpv*2-6)) & 3) * r1c);
                word_type rep3 = ~(((s0 >> (bpv*2-8)) & 3) * r1c);
                for (size_t i = pos/bpv; i < nv; ++i) {
                    // do this every bpv characters (usually 32)
                    word_type v0 = values[i];
                    word_type v1 = values[i+1];

                    // Each of the reps is xored with the next 32 characters.
                    // example:
                    //  ~AAAAAAAAAAAAAAAAAAAA
                    //   ACGAAGGGTATCGAGGTTAC
                    //   ====================
                    //   *  **    *   *    *  (*=11)
                    // ie. we get two '1' bits where characters match.
                    //
                    // So if we were searching for 'ATCG' the string above would give us:
                    //   ACGAAGGGT*TCGAGGTTAC
                    //   CGAAGGGTA*CGAGGTTAC
                    //   GAAGGGTAT*GAGGTTAC
                    //   AAGGGTATC*AGGTTAC
                    //
                    // ie. All four bits are '11'
                    // After grouping the 1's mask is left with a '1' only in bit positions
                    // where all four values match.
                    word_type mask = v0 ^ rep0;
                    mask &= ((v0 << 2 ) | (v1 >> (bpv*2-2))) ^ rep1;
                    mask &= ((v0 << 4 ) | (v1 >> (bpv*2-4))) ^ rep2;
                    mask &= ((v0 << 6 ) | (v1 >> (bpv*2-6))) ^ rep3;
                    mask &= mask << 1;
                    mask &= 0xaaaaaaaaaaaaaaaaull;

                    // usually this while loop is not executed.
                    size_t bit_pos = 0;
                    while (mask << bit_pos) {
                        int lz = (int)lzcnt(mask << bit_pos, cpu_has_lzcnt);
                        bit_pos += lz;
                        size_t search_pos = i * bpv + bit_pos/2;
                        word_type v = bit_pos == 0 ? v0 : (v0 << bit_pos) | (v1 >> (64-bit_pos));
                        ++bit_pos;
                        if (
                            ((v ^ s0) & s0mask) == 0 &&
                            search_pos >= pos && search_pos <= last - ssz &&
                            compare_inexact(search_pos, ssz, dna_str, max_distance) == 0
                        ) {
                            return search_pos;
                        }
                    }
                }

                pos = nv * bpv;
            } else {
                if (cpu_has_popcnt) {
                    pos = inexact_search<unmapped_traits, true>(dna_str, pos, nv, s0, s0mask, max_distance, ssz, last);
                } else {
                    pos = inexact_search<unmapped_traits, false>(dna_str, pos, nv, s0, s0mask, max_distance, ssz, last);
                }
                if (pos != basic_dna_string::npos) {
                    return pos;
                }
                pos = nv * bpv;
            }

            while (pos <= last - ssz) {
                word_type w0 = window(pos);
                if (count_word((w0 ^ s0) & s0mask, cpu_has_popcnt) <= max_distance) {
                    if (compare_inexact(pos, ssz, dna_str, max_distance) == 0) {
                        return pos;
                    }
                }
                ++pos;
            }
            return basic_dna_string::npos;
        }


        //! \brief Compare two substrings exactly.
        //! \tparam StringTraits Traits of other string to compare with.
        //! \param start_pos Zero-based offset to start the search.
        //! \param max_bases maxiumum number of bases to search.
        //! \param str dna_string to compare with.
        template <class StringTraits>
        int compare(size_t start_pos, size_t max_bases, const basic_dna_string<StringTraits> &str) const {
            return compare_inexact(start_pos, max_bases, str);
        }

        //! \brief return the distance (number of errors) between a string and a substring.
        //! \tparam StringTraits Traits of other string to compare with.
        //! \param start_pos Zero-based offset to start the search.
        //! \param max_bases maxiumum number of bases to search.
        //! \param str dna_string to compare with.
        template <class StringTraits>
        size_t distance(size_t start_pos, size_t max_bases, const basic_dna_string<StringTraits> &str) const {
            size_t pos = std::min(start_pos, num_bases);
            max_bases = std::min(max_bases, str.size());
            max_bases = std::min(max_bases, num_bases - pos);

            bool cpu_has_popcnt = has_popcnt();
            const auto str_values = str.get_values();
            const size_t bpv = bases_per_value;
            size_t nv = std::min(str_values.size(), (max_bases+bpv-1)/bpv);
            size_t error = 0;
            for (size_t i = 0; i != nv; ++i) {
                word_type w = window(pos);
                word_type s = str_values[i];
                if (i == max_bases/bpv) {
                    s &= ~(word_type)0 << (((0-max_bases) % bases_per_value) * 2);
                    w &= ~(word_type)0 << (((0-max_bases) % bases_per_value) * 2);
                }
                if (s != w) {
                    error += count_word(s^w, cpu_has_popcnt);
                }
                pos += bpv;
            }
            return error;
        }

        //! \brief Compare two substrings with errors.
        //! \tparam StringTraits Traits of other string to compare with.
        //! \param start_pos Zero-based offset to start the search.
        //! \param max_bases maxiumum number of bases to search.
        //! \param str dna_string to compare with.
        //! \param max_distance number of allowable errors in the search.
        template <class StringTraits>
        int compare_inexact(size_t start_pos, size_t max_bases, const basic_dna_string<StringTraits> &str, size_t max_distance=0) const {
            size_t pos = std::min(start_pos, num_bases);
            max_bases = std::min(max_bases, str.size());
            max_bases = std::min(max_bases, num_bases - pos);

            bool cpu_has_popcnt = has_popcnt();
            const auto str_values = str.get_values();
            const size_t bpv = bases_per_value;
            size_t nv = std::min(str_values.size(), (max_bases+bpv-1)/bpv);
            size_t error = 0;
            for (size_t i = 0; i != nv; ++i) {
                word_type w = window(pos);
                word_type s = str_values[i];
                if (i == max_bases/bpv) {
                    s &= ~(word_type)0 << (((0-max_bases) % bases_per_value) * 2);
                    w &= ~(word_type)0 << (((0-max_bases) % bases_per_value) * 2);
                }
                if (s != w) {
                    if (max_distance == 0) {
                        return s == w ? 0 : s < w ? -1 : 1;
                    } else {
                        error += count_word(s^w, cpu_has_popcnt);
                        if (error > max_distance) {
                            return s < w ? -1 : 1;
                        }
                    }
                }
                pos += bpv;
            }
            return 0;
        }

        //! \brief Get a value from 0..3 for a single base at offset "index".
        int get_code(size_t index) const {
            size_t sh = ((bases_per_value - 1 - index) % bases_per_value) * 2;
            size_t off = index/bases_per_value;
            return index >= num_bases ? 0 : ((values[off] >> sh) & 0x03);
        }

        //! \brief Set a value from 0..3 for a single base at offset "index".
        void set_code(size_t index, int code) {
            size_t sh = ((bases_per_value - 1 - index) % bases_per_value) * 2;
            size_t off = index/bases_per_value;
            word_type mask = (word_type)0x03 << sh;
            if (index < num_bases) {
                values[off] = (values[off] & ~mask) | (code << sh);
            }
        }

        //! \brief Get a right-justified word of values limited by num_index_chars.
        //! eg. get_index(pos, 3) gives AAA...AAAXXX
        //! Used by two_stage_index to index values
        word_type get_index(size_t pos, size_t num_index_chars) const {
            return (word_type)(window(pos) >> (bases_per_value - num_index_chars)*2);
        }

        //! \brief Get an unaligned word from the centre of the string (typically 32 values).
        word_type window(size_t base) const {
            size_t offset = base / bases_per_value;
            size_t sh = (base % bases_per_value) * 2;
            word_type v0 = offset < values.size() ? values[offset] : 0;
            if (sh == 0) {
                return v0;
            } else {
                size_t offset1 = (base + bases_per_value) / bases_per_value;
                word_type v1 = offset1 < values.size() ? values[offset1] : 0;
                return v0 << sh | v1 >> (64-sh);
            }
        }

        //! \brief Swap two dna_strings.
        void swap(basic_dna_string &rhs) {
            std::swap(num_bases, rhs.num_bases);
            values.swap(rhs.values);
        }

        //! \brief Write the structure in binary for later mapping.
        void write_binary(writer &wr) const {
            wr.write64(num_bases);
            wr.write(values);
        }

        //! \brief Back-door access to the values.
        const array_type &get_values() const {
            return values;
        }

        //! \brief Calculate the Burrows Wheeler Transform
        void bwt(basic_dna_string<unmapped_traits> &result, size_t &inverse_sa0) const {
            const bool debug = false;

            result.resize(size()+1);

            // The suffix array is a sort of all the substrings of a string.
            // Example "hello" -> "", "ello", "hello", llo", "lo", "o"
            // with the first character in alphabetical order
            // giving the positions 5, 1, 0, 2, 3, 4
            
            // We use the suffix array to generate a Burrows Wheeler transform
            // of the string, which is the character preceding each substring.
            // In this case this is 'o', 'h', '$', 'e', 'l', 'l'
            // '$' is a special character which is an 'A' in the result but has
            // location inverse_sa0 in the result.
            
            // see: https://en.wikipedia.org/wiki/FM-index
        
            size_t str_size = size();
            if (str_size == 0) {
                return;
            }
            if (
                size() >= std::numeric_limits<addr_type>::max()
            ) {
                throw std::invalid_argument(
                    "bwt(): string too large for address type"
                );
            }
            
            // Big note: it may seem strange to computer scientists to use sorts in this way
            // when you could invert an array using indices, but using x[y[i]]
            // constructs in real code is highly inefficient (~3000 cycles per
            // element for large arrays)

            // construction sorter.
            // note: this is impractical on real (non-academic) machines of the 2015
            // generation as it uses at least 12.25 bytes per base (~32GB for ENSEMBL human).
            // 64GB machines will be the norm fairly soon and in the meantime, SSE backed swap
            // space works reasonably well with std::sort.
            //
            // Richard Durbin's aproach of using BWT merging seems better and we should
            // investigate adopting it, especially as cloud computing is the target
            // with limited memory available at reasonable price.
            std::vector<bwt_sort_type> sorter;

            // allocate memory.
            sorter.reserve(str_size+1);
            result.resize(str_size+1);

            auto sort_by_address = [](const bwt_sort_type &a, const bwt_sort_type &b) { return a.addr < b.addr; };

            auto sort_by_group = [](const bwt_sort_type &a, const bwt_sort_type &b) {
                return a.group < b.group;
            };

            auto sort_by_group_and_next = [](const bwt_sort_type &a, const bwt_sort_type &b) {
                return a.group < b.group + (a.next_group < b.next_group);
            };
            
            auto debug_dump = [this,&sorter](const char *msg, size_t h) {
                printf("\n%s h=%d\n", msg, (int)h);
                for (size_t i = 0; i != sorter.size(); ++i) {
                    auto &s = sorter[i];
                    const char *str = std::string(substr(s.addr, 40)).c_str();
                    printf("%4d %12d %4d %4d %s\n", (int)i, (int)s.group, (int)s.next_group, (int)s.addr, str);
                }
            };

            // use the upper few bits of the dna string as a key for the initial sort
            const size_t key_chars = sizeof(addr_type) * 4;
            const int sh = sizeof(word_type) * 8 - key_chars * 2;
            for (size_t i = 0; i != str_size; ++i) {
                addr_type group = (addr_type)(window(i) >> sh);
                addr_type next_group = i + key_chars < str_size ? key_chars : str_size - i;
                sorter.push_back(bwt_sort_type{ group, next_group, (addr_type)i });
            }
            
            {
                // the $ value
                sorter.push_back(bwt_sort_type{ 0, 0, (addr_type)str_size });
            }

            for (size_t h = key_chars; ; h *= 2) {
                // sort by group    
                std::sort(sorter.begin(), sorter.end(), sort_by_group_and_next);
    
                if (debug) debug_dump("sorted group and next_group", h);

                // build the group numbers
                auto t = sorter[0];
                addr_type group = 0;
                sorter[0].group = group++;
                int finished = 1;
                for (size_t i = 1; i != str_size+1; ++i) {
                    auto &s = sorter[i];
                    int different = t.group != s.group | t.next_group != s.next_group;
                    finished &= different;
                    group = different ? (addr_type)i : group;
                    t = s;
                    s.group = group;
                }
                
                if (debug) debug_dump("built group numbers", h);

                if (finished) break;
    
                // invert the sort, restoring the original order of the sequence
                std::sort(sorter.begin(), sorter.end(), sort_by_address);

                // avoid random access to the string by using a sort.                
                for (size_t i = 0; i != str_size+1; ++i) {
                    bwt_sort_type &s = sorter[i];
                    s.next_group = i < str_size - h ? sorter[i+h].group : 0;
                }
    
                if (debug) debug_dump("set next_group", h);
            }

            // linearise the addresses for fast access.
            // this may be faster than using a peek/poke aproach
            std::sort(sorter.begin(), sorter.end(), sort_by_address);

            // this bwt entry (for address 0) would usually have '$',
            // but we have only four characters in the alphabet.
            inverse_sa0 = sorter[0].group;
            sorter[0].next_group = (addr_type)0;

            // put the previous char in next_group for the bwt.
            for (size_t i = 1; i != str_size+1; ++i) {
                bwt_sort_type &s = sorter[i];
                s.next_group = get_code(i-1);
            }

            if (debug) debug_dump("set next_group for bwt", 0);

            // restore the group order
            std::sort(sorter.begin(), sorter.end(), sort_by_group);

            // copy the chars into the BWT
            for (size_t i = 0; i != str_size+1; ++i) {
                bwt_sort_type &s = sorter[i];
                result.set_code(i, s.next_group);
            }

            if (debug) debug_dump("set bwt data", 0);
            
            if (debug) printf("inverse_sa0=%d\n", (int)inverse_sa0);
            
            // make sure we fail the tests if debug is on.
            if (debug) result.resize(0);
        }

        //! \brief Calculate the inverse Burrows Wheeler Transform
        void ibwt(basic_dna_string<unmapped_traits> &result, size_t inverse_sa0) const {
            size_t str_size = size();
            result.resize(str_size-1);
            
            std::vector<ibwt_sort_type> sorter;
            sorter.resize(str_size);
            for (size_t i = 0; i != str_size; ++i) {
                sorter[i].addr = (addr_type)i;
                sorter[i].chr = get_code(i) + 1;
            }
            sorter[inverse_sa0].chr = 0;

            // Ok, we should be able to bucket sort this one!
            auto sort_by_chr = [](const ibwt_sort_type &a, const ibwt_sort_type &b) { return a.chr < b.chr + (a.addr < b.addr); };
            std::sort(sorter.begin(), sorter.end(), sort_by_chr);

            // Sadly, very little choice but to pointer chase this loop.
            for (size_t addr = sorter[0].addr, i = 0; i != str_size; ++i, addr = sorter[addr].addr) {
                result.set_code(i, sorter[i].chr - 1);
            }
        }
    private:
        //! Inexact search using popcnt (if we have one!)
        template <class StringTraits, bool cpu_has_popcnt>
        size_t inexact_search(basic_dna_string<StringTraits> &search_str, size_t pos, size_t nv, word_type s0, word_type s0mask, size_t max_distance, size_t max_bases, size_t last) const {
            const size_t bpv = bases_per_value;
            for (size_t i = pos/bpv; i < nv; ++i) {
                word_type v0 = values[i];
                word_type v1 = values[i+1];
                word_type s0x = s0;
                int have_hits = 0;
                #define BOOST_GENETICS_UNROLL \
                    have_hits |= popcnt((v0 ^ s0x) & s0mask, cpu_has_popcnt) - (int)(max_distance*2+1); \
                    v0 = (v0 << 2) | v1 >> (bpv*2-2); \
                    v1 <<= 2;
                for (size_t j = 0; j != bpv/4; ++j) {
                    BOOST_GENETICS_UNROLL BOOST_GENETICS_UNROLL BOOST_GENETICS_UNROLL BOOST_GENETICS_UNROLL
                }
                #undef BOOST_GENETICS_UNROLL
                if (have_hits < 0) {
                    for (size_t j = 0; j != bpv; ++j) {
                        size_t search_pos = i * bpv + j;
                        if (
                            search_pos >= pos && search_pos <= last - max_bases &&
                            compare_inexact(search_pos, max_bases, search_str, max_distance) == 0
                        ) {
                            return search_pos;
                        }
                    }
                }
            }

            return basic_dna_string::npos;
        }

        // Note: order matters as constructor builds first value first.
        size_t num_bases;
        array_type values;
    };

    //! \brief Write the dna_string as a stream of ASCII characters.
    template <class charT, class traits, class DnaTraits>
    std::basic_istream<charT, traits>&
    operator>>(std::basic_istream<charT, traits>& is, basic_dna_string<DnaTraits>& x) {
        std::string str;
        is >> str;
        x = basic_dna_string<DnaTraits>(str);
        return is;
    }

    //! \brief Read the dna_string from a stream of ASCII characters.
    template <class charT, class traits, class DnaTraits>
    std::basic_ostream<charT, traits>&
    operator<<(std::basic_ostream<charT, traits>& os, const basic_dna_string<DnaTraits>&x) {
        os << (std::string)x;
        return os;
    }

    //! \brief Reverse complement the string. Converts T<->A C<->G and reverses the string.
    template <class DnaTraits>
    basic_dna_string<DnaTraits> rev_comp(const basic_dna_string<DnaTraits> &x) {
        return x.substr(0, x.size(), true);
    }

    //! \brief Conventionally allocated dna string used for construction.
    typedef basic_dna_string<unmapped_traits> dna_string;

    //! \brief File mapped dna string used for searches.
    typedef basic_dna_string<mapped_traits> mapped_dna_string;

    template <>
    inline int get_code<dna_string>(const dna_string &str, size_t index) {
        return str.get_code(index);
    }
} }


#endif
