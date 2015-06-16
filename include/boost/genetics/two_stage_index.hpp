// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_TWO_STAGE_INDEX_HPP
#define BOOST_GENETICS_TWO_STAGE_INDEX_HPP

#include <boost/genetics/augmented_string.hpp>

#include <stdexcept>

namespace boost { namespace genetics {
    /// two stage index, first index ordered by value, second by address.
    template <class StringType, class IndexArrayType, class AddrArrayType, bool Writable>
    class basic_two_stage_index {
    public:
        typedef typename IndexArrayType::value_type index_type;
        typedef typename AddrArrayType::value_type addr_type;

        basic_two_stage_index(
        ) : string(nullptr), num_indexed_chars(0) {
        }
        
        void write_binary(writer &wr) const {
            wr.write64(num_indexed_chars);
            wr.write(index);
            wr.write(addr);
        }

        basic_two_stage_index &operator =(basic_two_stage_index &&rhs) {
            string = rhs.string;
            index = std::move(rhs.index);
            addr = std::move(rhs.addr);
            num_indexed_chars = rhs.num_indexed_chars;
            return *this;
        }

        template <class Mapper>
        basic_two_stage_index(
            StringType &string,
            Mapper &map,
            typename Mapper::is_mapper *p=0
        ) :
            string(&string),
            num_indexed_chars((size_t)map.read64()),
            index(map),
            addr(map)
        {
        }
        
        basic_two_stage_index(
            StringType &string,
            size_t num_indexed_chars
        ) : string(&string), num_indexed_chars(num_indexed_chars) {
            reindex();
        }

        void reindex() {
            reindex_impl<Writable>();
        }
        
        size_t end() const {
            return (size_t)-1;
        }

        class iterator {
        public:
            iterator() {
            }

            iterator(const basic_two_stage_index *tsi, const dna_string& search_str, size_t min_pos, size_t max_distance, size_t max_gap, bool is_brute_force) :
                tsi(tsi), search_str(search_str), max_distance(max_distance), max_gap(max_gap), is_brute_force(is_brute_force)
            {
                num_indexed_chars = tsi->num_indexed_chars;
                num_seeds = search_str.size() / num_indexed_chars;
                this->is_brute_force = is_brute_force || num_seeds <= max_distance;
                merges_done = 0;
                compares_done = 0;

                if (this->is_brute_force) {
                    pos = min_pos;
                    find_next(true);
                    return;
                } else {
                    active.resize(num_seeds);
                    if (num_seeds == 0) {
                        pos = dna_string::npos;
                        return;
                    }

                    //long long t0 = __rdtsc();
                    // get seed values from string
                    // touch the index in num_seeds places (with NTA hint)
                    for (size_t i = 0; i != num_seeds; ++i) {
                        active_state &s = active[i];
                        s.idx = (index_type)get_index(search_str, i * num_indexed_chars, num_indexed_chars);
                        touch_nta(tsi->addr.data() + tsi->index[i]);
                    }

                    // get seed values from string
                    // touch the address vector in num_seeds places (with stream hint)
                    for (size_t i = 0; i != num_seeds; ++i) {
                        active_state &s = active[i];
                        const addr_type *ptr = tsi->addr.data() + tsi->index[s.idx];
                        const addr_type *end = tsi->addr.data() + tsi->index[s.idx+1];
                        s.ptr = ptr;
                        s.end = end;
                        if (ptr != end) touch_stream(ptr);
                    }

                    for (size_t i = 0; i != num_seeds; ++i) {
                        active_state &s = active[i];
                        const addr_type *ptr = s.ptr;
                        const addr_type *end = s.end;
                        s.elem = (addr_type)i;
                        addr_type skip = (addr_type)(i * num_indexed_chars + min_pos);
                        while (ptr != end && *ptr < skip) {
                            ++ptr;
                        }
                        s.start = ptr == end ? (addr_type)-1 : (addr_type)(*ptr - (s.elem * num_indexed_chars));
                        s.prev = (addr_type)-1;
                        s.ptr = ptr;
                    }

                    //long long t1 = __rdtsc();
                    //printf("touch %lld\n", t1 - t0);

                    std::make_heap(active.begin(), active.end());
                    find_next(true);
                }
            }

            operator size_t() const {
                return pos;
            }

            iterator &operator++() {
                find_next(false);
                return *this;
            }

            iterator &operator++(int) {
                find_next(false);
                return *this;
            }
        private:
            void find_next(bool is_start) {
                if (pos == dna_string::npos) {
                    // if we have already reached the end, stop.
                    return;
                } else if (is_brute_force) {
                    size_t start = is_start ? pos : pos + 1;
                    if (start + search_str.size() > tsi->string->size()) {
                        pos = dna_string::npos;
                    } else {
                        pos = tsi->string->find_inexact(search_str, start, ~(size_t)0, max_distance);
                    }
                    return;
                } else {
                    // for a small number of unknowns, use a merge to find potential starts.
                    addr_type prev_start = (addr_type)-1;
                    size_t repeat_count = 0;
                    pos = dna_string::npos;

                    for (;;) {
                        merges_done++;
                        active_state s = active.front();

                        if (s.start != prev_start) {
                            size_t erc = 0;
                            if (repeat_count >= num_seeds - max_distance) {
                                compares_done++;
                                if (tsi->string->compare_inexact(prev_start, search_str.size(), search_str, max_distance) == 0) {
                                  pos = prev_start;
                                  return;
                                }
                            }
                            repeat_count = 0;
                            prev_start = s.start;
                        }

                        if (s.start == (addr_type)-1) {
                            //printf("merges_done=%lld\n", (long long)merges_done);
                            //printf("compares_done=%lld\n", (long long)compares_done);
                            return;
                        }

                        // printf("xxx %5d %2d/%2d [%2d]\n", (int)s.start, (int)repeat_count, (int)(num_seeds - max_distance), (int)s.elem);

                        const addr_type *ptr = s.ptr + 1;
                        const addr_type *end = s.end;
                        s.prev = s.start;
                        s.start = ptr == end ? (addr_type)-1 : (addr_type)(*ptr - s.elem * num_indexed_chars);
                        s.ptr = ptr;
                        std::pop_heap(active.begin(), active.end());
                        active.back() = s;
                        std::push_heap(active.begin(), active.end());
                        
                        // printf("xxx next: %d\n", (int)active.front().start);
                        repeat_count++;
                    }
                }
            }

            // search state
            struct active_state {
                const addr_type *ptr;
                const addr_type *end;
                index_type idx;
                addr_type start;
                addr_type prev;
                addr_type elem;

                bool operator<(const active_state &rhs) {
                    return start > rhs.start;
                }
            };

            // index used
            const basic_two_stage_index *tsi;

            // array of active pointers for each seed
            std::vector<active_state> active;

            // dna string
            const dna_string &search_str;

            // max allowable errors
            size_t max_distance;

            // maximum gaps allowed (introns)
            size_t max_gap;

            // current search position.
            size_t pos;

            // do a linear scan of the indexed string (for testing)
            bool is_brute_force;

            // number of seeds in search string.
            size_t num_seeds;

            // number of chars per index location
            size_t num_indexed_chars;

            // measure work done by the algorithm.
            size_t merges_done;
            size_t compares_done;
        };

        /// find the next dna string which is close to the search string allowing max_distance errors and max_gap gaps between exons.
        iterator find_inexact(const dna_string& search_str, size_t pos, size_t max_distance, size_t max_gap, bool is_brute_force) const {
            return iterator(this, search_str, pos, max_distance, max_gap, is_brute_force);
        }

        template <class charT, class traits>
        void write_ascii(std::basic_ostream<charT, traits>& os) const {
            typename std::basic_ostream<charT, traits>::fmtflags save = os.flags();
            size_t str_size = string->size();
            size_t index_size = (size_t)1 << (num_indexed_chars*2);
            for (size_t i = 0; i != index_size; ++i) {
                for (index_type j = index[i]; j != index[i+1]; ++j) {
                    addr_type a = addr[j];
                    os << "[" << a << " " << string->substr(a, num_indexed_chars) << "]";
                }
                os << "\n";
            }
            os.setstate( save );
        }
        
        void swap(basic_two_stage_index &rhs) {
            std::swap(string, rhs.string);
            std::swap(num_indexed_chars, rhs.num_indexed_chars);
            index.swap(rhs.index);
            addr.swap(rhs.addr);
        }
    private:
        // mapped version has no reindex() as we can't write to the vectors
        template<bool B>
        void reindex_impl() {
            throw (std::exception("two_stage_index::reindex(): attempt to reindex a read-only index"));
        }

        // std::vector version has reindex()
        template<>
        void reindex_impl<true>() {
            // todo: use threads to speed this up.
            //       consider using std::sort for better cache performance.
            if (
                num_indexed_chars <= 1 ||
                num_indexed_chars > 32 ||
                string->size() < num_indexed_chars ||
                (addr_type)string->size() != string->size()
            ) {
                throw std::invalid_argument("two_stage_index::reindex()");
            }

            size_t str_size = string->size();
            size_t index_size = (size_t)1 << (num_indexed_chars*2);

            index.resize(0);
            index.resize(index_size+1);
            addr.resize(str_size);

            // lead-in: fill acc0 with first DNA codes
            size_t acc0 = 0;
            for (size_t i = 0; i != num_indexed_chars-1; ++i) {
                acc0 = acc0 * 4 + get_code(*string, i);
            }

            // count phase: count bucket sizes
            size_t acc = acc0;
            for (size_t i = num_indexed_chars; i != str_size; ++i) {
                acc = (acc * 4 + get_code(*string, i)) & (index_size-1);
                index[acc]++;
                if (i % 0x100000 == 0) printf("c %5.2f\n", (double)i * (100.0/str_size));
            }
            
            // compute running sum of buckets to convert counts to offsets.
            addr_type cur = 0;
            for (size_t i = 0; i != index_size; ++i) {
                addr_type val = index[i];
                index[i] = cur;
                cur += val;
            }
            index[index_size] = cur;

            // store phase: fill "addr" with addresses of values.
            acc = acc0;
            for (size_t i = num_indexed_chars; i != str_size; ++i) {
                acc = (acc * 4 + get_code(*string, i)) & (index_size-1);
                addr[index[acc]++] = (addr_type)(i - num_indexed_chars + 1);
                if (i % 0x100000 == 0) printf("s %5.2f\n", (double)i * (100.0/str_size));
            }

            // shift the index up one so ends become starts.
            addr_type prev = 0;
            for (size_t i = 0; i != index_size; ++i) {
                std::swap(prev, index[i]);
            }
        }

        // note: order matters
        StringType *string;
        size_t num_indexed_chars;
        IndexArrayType index;
        AddrArrayType addr;
    };

    typedef basic_two_stage_index<augmented_string, std::vector<uint32_t>, std::vector<uint32_t>, true > two_stage_index;
    typedef basic_two_stage_index<mapped_augmented_string, mapped_vector<uint32_t>, mapped_vector<uint32_t>, false > mapped_two_stage_index;

    template <class charT, class traits, class StringType, class IndexArrayType, class AddrArrayType, bool Writable>
    std::basic_ostream<charT, traits>&
    operator<<(std::basic_ostream<charT, traits>& os, const basic_two_stage_index<StringType, IndexArrayType, AddrArrayType, Writable>& x) {
        x.write_ascii(os);
        return os;
    }
} }


#endif
