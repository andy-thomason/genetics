// Copyright Andy Thomason 2016
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_SUFFIX_ARRAY
#define BOOST_GENETICS_SUFFIX_ARRAY

#include <stdexcept>
#include <type_traits>
#include <limits>
#include <boost/genetics/dna_string.hpp>

namespace boost { namespace genetics {
    //! \brief This class is a suffix array that may be used to index a DNA string.

    //! Like many of the container classes in this library it can be specialised
    //! into a standard (`std::vector`) version for construction and a read-only
    //! mapped (`mapped_vector`) version for high performance use.

    //! \tparam Traits Typically one of unmapped_traits or mapped_traits.
    //! These may be customised for larger address lengths, for example.

    template<class Traits>
    class basic_suffix_array {
    public:
        typedef basic_dna_string<Traits> string_type;
        typedef typename Traits::SuffixArrayType addr_array_type;
        typedef typename Traits::SuffixArrayType::value_type addr_type;
        typedef typename string_type::word_type word_type;

        struct sorter_t { addr_type group, next_group, addr; };

        static_assert(
            sizeof(word_type) >= sizeof(addr_type),
            "address type must be smaller than the string word type"
        );

        //! \brief Construct an empty suffix array.
        basic_suffix_array(
        ) : string_(nullptr) {
        }
        
        //! \brief Write to a binary stream for subsequent mapping.
        void write_binary(writer &wr) const {
            wr.write(addr_);
        }

        //! \brief rvalue move operator
        basic_suffix_array &operator =(basic_suffix_array &&rhs) {
            string_ = rhs.string_;
            addr_ = std::move(rhs.addr_);
            return *this;
        }

        //! \brief Construct a suffix array from a persistant block of memory.
        //! This will liely be a memory mapped file.
        template <class Mapper>
        basic_suffix_array(
            string_type &string,
            Mapper &map,
            typename Mapper::is_mapper *p=0
        ) :
            string_(&string),
            addr_(map)
        {
        }
        
        //! \brief Construct a suffix array from a dna_string or mapped_dna_string.
        basic_suffix_array(
            string_type &str
        ) : string_(&str) {
            size_t str_size = string_->size();
            if (str_size == 0) {
                return;
            }
            if (
                string_->size() >= std::numeric_limits<addr_type>::max()
            ) {
                throw std::invalid_argument(
                    "basic_suffix_array() string too large for address type"
                );
            }

            // construction sorter.
            std::vector<sorter_t> sorter;
            sorter.reserve(str_size);

            // use the upper few bits of the dna string as a key for the initial sort
            const size_t key_chars = sizeof(addr_type) * 4;
            const int sh = sizeof(word_type) * 8 - key_chars * 2;
            for (size_t i = 0; i != str_size; ++i) {
                addr_type group = (addr_type)(str.window(i) >> sh);
                addr_type next_group = i + key_chars < str_size ? key_chars : str_size - i;
                sorter_t rec = { group, next_group, (addr_type)i };
                sorter.push_back(rec);
            }

            for (size_t h = key_chars; ; h *= 2) {
                // sort by group    
                auto fwd_two_value_sort = [](const sorter_t &a, const sorter_t &b) {
                    return a.group < b.group + (a.next_group < b.next_group);
                };
                std::sort(sorter.begin(), sorter.end(), fwd_two_value_sort);
    
                if (debug) debug_dump(sorter, "sorted group and next_group", h);

                // build the group numbers
                auto t = sorter[0];
                addr_type group = 0;
                sorter[0].group = group;
                int finished = 1;
                for (size_t i = 1; i != str_size; ++i) {
                    auto &s = sorter[i];
                    int different = t.group != s.group | t.next_group != s.next_group;
                    finished &= different;
                    group = different ? (addr_type)i : group;
                    t = s;
                    s.group = group;
                }
                
                if (debug) debug_dump(sorter, "built group numbers", h);

                if (finished) break;
    
                // invert the sort, restoring the original order of the sequence
                auto sort_by_address = [](const sorter_t &a, const sorter_t &b) { return a.addr < b.addr; };
                std::sort(sorter.begin(), sorter.end(), sort_by_address);
                
                for (size_t i = 0; i != str_size; ++i) {
                    sorter_t &s = sorter[i];
                    s.next_group = i < str_size - h ? sorter[i+h].group + 1 : 0;
                }
    
                if (debug) debug_dump(sorter, "set next_group", h);
            }

            addr_.reserve(str_size);
            for (size_t i = 0; i != str_size; ++i) {
                sorter_t &s = sorter[i];
                addr_.push_back(s.addr);
            }
        }

        void swap(basic_suffix_array &rhs) {
            std::swap(string_, rhs.string_);
            addr_.swap(rhs.addr_);
        }

        //! Write the suffix array to a stream for debugging.
        template <class charT, class traits>
        void write_ascii(std::basic_ostream<charT, traits>& os) const {
            for (size_t i = 0; i != addr_.size(); ++i) {
                addr_type a = addr_[i];
                os << string_->substr(a, 40) << "\n";
            }
        }
        
        //! Verify that the suffix array is more or less correct.
        bool verify() {
            size_t asize = addr_.size();
            if (string_->size() != asize) {
                std::cout << "fail: size mismatch\n";
                return false;
            }
            
            if (addr_.size() == 0) {
                return true;
            }
            
            for (size_t i = 1; i != asize; ++i) {
                addr_type a0 = addr_[i-1];
                addr_type a1 = addr_[i];
                word_type w0 = string_->window(a0);
                word_type w1 = string_->window(a1);
                if (w0 > w1) {
                    size_t l0 = std::min((size_t)(asize - a0), (size_t)32);
                    size_t l1 = std::min((size_t)(asize - a1), (size_t)32);
                    std::cout << "fail: order mismatch" << to_dna(w0, l0) << " " << to_dna(w1, l1) << "\n";
                    return false;
                }
                // todo: test more values!
            }

            if (debug) {
                std::cout << "fail: debug left on\n";
                return false;
            }
            
            return true;
        }
        
        //! /brief Indexing operator, gives the offset of the i'th substring
        //! in lexical order. The zero'th element is not the end of the string
        //! but the lowest non-terminal substring.
        addr_type operator[](size_t index) const {
            return addr_[index];
        }
        
        //! /brief Size of the suffix array not counting the implicit zero's element.
        //! Thus this is also the same length as the string.
        size_t size() const {
            return addr_.size();
        }
    private:
        // change this to "true" to diagnose errors.
        static const bool debug = false;
        void debug_dump(const std::vector<sorter_t> &sorter, const char *msg, size_t h) {
            printf("\n%s h=%d\n", msg, (int)h);
            for (size_t i = 0; i != sorter.size(); ++i) {
                auto &s = sorter[i];
                const char *str = std::string(string_->substr(s.addr, 40)).c_str();
                printf("%4d %4d %4d %4d %s\n", (int)i, (int)s.group, (int)s.next_group, (int)s.addr, str);
            }
        }

        // Note: order matters
        string_type *string_;
        addr_array_type addr_;
    };

    //! \brief Unmapped suffix array, typically used for construction.
    typedef basic_suffix_array<unmapped_traits> suffix_array;

    //! \brief Mapped suffix array, typically used for high performance loading.
    typedef basic_suffix_array<mapped_traits> mapped_suffix_array;

    //! \brief Stream write operator
    template <class charT, class traits, class Traits>
    std::basic_ostream<charT, traits>&
    operator<<(std::basic_ostream<charT, traits>& os, const basic_suffix_array<Traits>& x) {
        x.write_ascii(os);
        return os;
    }
} }


#endif
