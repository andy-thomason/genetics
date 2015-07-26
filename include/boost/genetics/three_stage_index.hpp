// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_THREE_STAGE_INDEX_HPP
#define BOOST_GENETICS_THREE_STAGE_INDEX_HPP

#include <stdexcept>
#include <type_traits>
#include <boost/genetics/augmented_string.hpp>


namespace boost { namespace genetics {
    /// Three stage index, first index ordered by value, second variable length, third by address.
    /// Note: work in progress!
    template<class Traits>
    class basic_three_stage_index {
    public:
        typedef basic_augmented_string<Traits> string_type;
        typedef typename Traits::TsiIndexArrayType index_array_type;
        typedef typename Traits::TsiAddrArrayType addr_array_type;
        typedef typename Traits::TsiIndexArrayType::value_type index_type;
        typedef typename Traits::TsiAddrArrayType::value_type addr_type;

        basic_three_stage_index(
        ) : string(nullptr), num_indexed_chars(0) {
        }
        
        void write_binary(writer &wr) const {
            wr.write64(num_indexed_chars);
            wr.write(index);
            wr.write(addr);
        }

        basic_three_stage_index &operator =(basic_three_stage_index &&rhs) {
            string = rhs.string;
            index = std::move(rhs.index);
            addr = std::move(rhs.addr);
            num_indexed_chars = rhs.num_indexed_chars;
            return *this;
        }

        template <class Mapper>
        basic_three_stage_index(
            string_type &string,
            Mapper &map,
            typename Mapper::is_mapper *p=0
        ) :
            string(&string),
            num_indexed_chars((size_t)map.read64()),
            index(map),
            addr(map)
        {
        }
        
        basic_three_stage_index(
            string_type &str,
            size_t num_indexed_chars
        ) : string(&str), num_indexed_chars(num_indexed_chars) {
            // Todo: use threads to speed this up.
            //       consider using std::sort for better cache performance.
            if (
                num_indexed_chars <= 1 ||
                num_indexed_chars > 32 ||
                string->size() < num_indexed_chars ||
                (addr_type)string->size() != string->size()
            ) {
                throw std::invalid_argument("three_stage_index::reindex()");
            }

            size_t str_size = string->size();
            size_t index_size = (size_t)1 << (num_indexed_chars*2);

            index.resize(0);
            index.resize(index_size+1);

            std::uint32_t acc0 = 0;
            for (size_t i = 0; i != 32/2-1; ++i) {
                acc0 = acc0 * 4 + get_code(*string, i);
            }

            std::vector<std::uint64_t> sorter;
            for (size_t prefix = 1; prefix != 32; ++prefix) {
                std::cerr << "p" << prefix << "\n";
                std::uint32_t acc = acc0;
                const int sh = 32-5;
                sorter.resize(0);
                for (size_t i = 32/2; i != str_size; ++i) {
                    acc = (acc * 4 + get_code(*string, i));
                    if ((acc >> sh) == prefix) {
                        sorter.push_back(((std::uint64_t)acc << (64-sh)) | (i - 32/2));
                    }
                }
                std::cerr << "s" << prefix << "\n";
                std::sort(sorter.data(), sorter.data() + sorter.size());

                for (size_t i = 0, imax = (sorter.size(),0x100); i != imax; ++i) {
                    printf("%llx\n", sorter[i]);
                }
            exit(1);
            }
            exit(1);
        }

        size_t end() const {
            return (size_t)-1;
        }

        class iterator {
        public:
            iterator() {
            }

            iterator(const basic_three_stage_index *tsi, const std::string& search_str, size_t min_pos, search_params &params) :
                tsi(tsi), search_str(search_str), params(params)
            {
                dna_search_str = search_str;
                num_indexed_chars = tsi->num_indexed_chars;
                size_t max_seeds = search_str.size() / num_indexed_chars;
                if (max_seeds <= params.max_distance) {
                    is_brute_force = true;
                    if (params.never_brute_force) {
                        pos = dna_string::npos;
                        return;
                    }
                } else {
                    is_brute_force = params.always_brute_force;
                }

                if (this->is_brute_force) {
                    pos = min_pos;
                    find_next(true);
                    return;
                } else {
                    active.resize(0);
                    if (max_seeds == 0) {
                        pos = dna_string::npos;
                        return;
                    }

                    //long long t0 = __rdtsc();
                    // Get seed values from string.
                    // Touch the index in num_seeds places (with NTA hint).
                    // If there are any 'N's in a seed, 
                    const char *str = search_str.data();
                    size_t total_N = 0;
                    size_t index_size = (size_t)1 << (num_indexed_chars*2);
                    size_t poly_A = 0;
                    for (size_t i = 0; i != max_seeds; ++i) {
                        const char *b = str + i * num_indexed_chars;
                        const char *e = b + num_indexed_chars;
                        size_t num_N = std::count(b, e, 'N');
                        total_N += num_N;
                        if (num_N == 0) {
                            active_state s;
                            s.elem = (addr_type)i;
                            s.idx = (index_type)get_index(dna_search_str, i * num_indexed_chars, num_indexed_chars);
                            if (s.idx != poly_A) {
                                //touch_nta(tsi->addr.data() + tsi->index[i]);
                                active.push_back(s);
                            }
                        }
                    }

                    // Get seed values from string.
                    // Touch the address vector in num_seeds places (with stream hint).
                    for (size_t i = 0; i != active.size(); ++i) {
                        active_state &s = active[i];
                        const addr_type *ptr = tsi->addr.data() + tsi->index[s.idx];
                        const addr_type *end = tsi->addr.data() + tsi->index[s.idx+1];
                        s.ptr = ptr;
                        s.end = end;
                        //if (ptr != end) touch_stream(ptr);
                    }

                    std::sort(
                        active.begin(), active.end(),
                        [](active_state &a, active_state &b) {
                            return a.end - a.ptr < b.end - b.ptr; 
                        }
                    );

                    for (size_t i = 0; i != active.size(); ++i) {
                        active_state &s = active[i];
                        if (s.end - s.ptr > 100) {
                            active.resize(i);
                            break;
                        }
                    }

                    /*if (active.size() > params.max_distance + 1) {
                        active.resize(params.max_distance + 1);
                    }*/

                    /*char tmp[1024], *p = tmp;
                    for (size_t i = 0; i != active.size(); ++i) {
                        active_state &s = active[i];
                        p += sprintf(p, "%d,", (int)(s.end - s.ptr));
                    }
                    puts(tmp);*/

                    /*for (size_t i = 0; i != active.size(); ++i) {
                        active_state &s = active[i];
                        if (s.end - s.ptr > 100) {
                            pos = dna_string::npos;
                            return;
                        }
                    }*/

                    for (size_t i = 0; i != active.size(); ++i) {
                        active_state &s = active[i];
                        const addr_type *ptr = s.ptr;
                        const addr_type *end = s.end;
                        addr_type skip = (addr_type)(s.elem * num_indexed_chars + min_pos);
                        while (ptr != end && *ptr < skip) {
                            ++ptr;
                        }
                        s.start = ptr == end ? (addr_type)-1 : (addr_type)(*ptr - (s.elem * num_indexed_chars));
                        s.prev = (addr_type)-1;
                        s.ptr = ptr;
                    }

                    if (active.size() <= params.max_distance) {
                        pos = dna_string::npos;
                        return;
                    }

                    required_seed_matches = active.size() - params.max_distance;
                    max_error = search_str.size() - params.max_distance - total_N;

                    std::make_heap(active.begin(), active.end());
                    pos = min_pos;
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

            size_t distance() const {
                return distance_;
            }
        private:
            void find_next(bool is_start) {
                if (pos == dna_string::npos) {
                    // If we have already reached the end, stop.
                    return;
                } else if (is_brute_force) {
                    size_t start = is_start ? pos : pos + 1;
                    if (start + search_str.size() > tsi->string->size()) {
                        pos = dna_string::npos;
                    } else {
                        pos = tsi->string->find_inexact(search_str, start, ~(size_t)0, params.max_distance);
                    }
                    return;
                } else {
                    // For a small number of unknowns, use a merge to find potential starts.
                    addr_type prev_start = (addr_type)-1;
                    size_t repeat_count = 0;
                    pos = dna_string::npos;

                    for (;;) {
                        params.stats.merges_done++;
                        active_state s = active.front();

                        if (s.start != prev_start) {
                            if (repeat_count >= required_seed_matches) {
                                params.stats.compares_done++;
                                distance_ = tsi->string->distance(prev_start, dna_search_str.size(), dna_search_str);
                                if (distance_ <= max_error) {
                                    // todo: check search_str also and don't count 'N's as error.
                                    pos = prev_start;
                                    return;
                                }
                            }
                            repeat_count = 0;
                            prev_start = s.start;
                        }

                        if (s.start == (addr_type)-1) {
                            return;
                        }

                        const addr_type *ptr = s.ptr + 1;
                        const addr_type *end = s.end;
                        s.prev = s.start;
                        s.start = ptr == end ? (addr_type)-1 : (addr_type)(*ptr - s.elem * num_indexed_chars);
                        s.ptr = ptr;
                        std::pop_heap(active.begin(), active.end());
                        active.back() = s;
                        std::push_heap(active.begin(), active.end());
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
            const basic_three_stage_index *tsi;

            // array of active pointers for each seed
            std::vector<active_state> active;

            // dna string
            const std::string &search_str;
            dna_string dna_search_str;

            // how to do this search.
            search_params &params;

            // using brute force search.
            bool is_brute_force;

            // distance (ie. error) between search string and current position.
            size_t distance_;

            // current search position.
            size_t pos;

            // Total maximum error, including 'N's.
            size_t max_error;

            // Required number of seed matches.
            size_t required_seed_matches;

            // number of chars per index location
            size_t num_indexed_chars;
        };

        /// find the next dna string which is close to the search string allowing max_distance errors and max_gap gaps between exons.
        iterator find_inexact(const std::string& search_str, size_t pos, search_params &params) const {
            return iterator(this, search_str, pos, params);
        }

        template <class charT, class traits>
        void write_ascii(std::basic_ostream<charT, traits>& os) const {
            typename std::basic_ostream<charT, traits>::fmtflags save = os.flags();
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
        
        void swap(basic_three_stage_index &rhs) {
            std::swap(string, rhs.string);
            std::swap(num_indexed_chars, rhs.num_indexed_chars);
            index.swap(rhs.index);
            addr.swap(rhs.addr);
        }

    private:
        // given a number of entries, how many extra words do we need?
        static int num_extra_words(int lg_count) {
            int lg_bytes_per_entry = std::min(lg_count >> 3, 2);
            int extra_chars = (lg_count-2) >> 1;
            int lg_entries = extra_chars * 2;
            return lg_count < 3 ? 0 : (1 << (lg_bytes_per_entry+lg_entries-2));
        }


        // Note: order matters
        string_type *string;
        size_t num_indexed_chars;
        index_array_type index;
        addr_array_type addr;
    };

    typedef basic_three_stage_index<unmapped_traits> three_stage_index;
    typedef basic_three_stage_index<mapped_traits> mapped_three_stage_index;

    template <class charT, class traits, class Traits>
    std::basic_ostream<charT, traits>&
    operator<<(std::basic_ostream<charT, traits>& os, const basic_three_stage_index<Traits>& x) {
        x.write_ascii(os);
        return os;
    }
} }


#endif
