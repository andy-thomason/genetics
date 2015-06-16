// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_AUGMENTED_STRING_HPP
#define BOOST_GENETICS_AUGMENTED_STRING_HPP

#include <boost/genetics/dna_string.hpp>

namespace boost { namespace genetics {
    /// containter for bases ACGT and occasional runs of 'N' and other letters.
    template<
        class WordType, class ParentType,
        class IndexArrayType, class RleArrayType
    >
    class basic_augmented_string : public ParentType {
        static const size_t lg_bases_per_index = 16;
        static const size_t bases_per_index =
            (size_t)1 << lg_bases_per_index;
    public:
        typedef typename IndexArrayType::value_type index_type;
        typedef typename RleArrayType::value_type rle_type;
        typedef ParentType parent;

        basic_augmented_string() {
        }

        basic_augmented_string(const char *str) {
            const char *e = str;
            while (*e) ++e;
            append(str, e);
        }

        template<class StringType>
        basic_augmented_string(
            const StringType &str,
            size_t pos = 0, size_t n = ~(size_t)0
        ) : parent() {
            append(str.begin() + pos, str.begin() + std::min(n, str.size()));
        }

        template <class Mapper>
        basic_augmented_string(Mapper &map, typename Mapper::is_mapper *p=0) :
            parent(map),
            index(map),
            rle(map)
        {
        }

        char operator[](size_t base) const {
            size_t offset = base / bases_per_index;
            if (offset > index.size()) return 'N';

            index_type i0 = index[offset];

            index_type i1 =
              offset+1 >= index.size() ?
              (index_type)rle.size() :
              index[offset+1]
            ;

            rle_type search = (rle_type)((base % bases_per_index) << 8) | 0xff;
            const rle_type *b = rle.data() + i0;
            const rle_type *e = rle.data() + i1;
            const rle_type *p = std::lower_bound(b, e, search);
            return p > b && (p[-1] & 0xff) ? (p[-1] & 0xff) : ((const parent&)*this)[base];
        }

        void append(const char *str) {
            const char *e = str;
            while (*e) ++e;
            append(str, e);
        }

        void resize(size_t new_size, char chr='A') {
            size_t old_size = parent::size();
            parent::resize(new_size);
            if (new_size > old_size) {
                while (old_size < new_size) {
                    internal_append(old_size, &chr, &chr+1);
                    old_size++;
                }
            } else {
                index.resize(new_size / bases_per_index);
                rle.resize(index.empty() ? 0 : index.back());
            }
        }

        operator std::string() const {
            std::string res;
            res.resize(parent::size());
            for (size_t i = 0; i != parent::size(); ++i) {
                res[i] = (*this)[i];
            }
            return res;
        }

        std::string substr(
            size_t offset=0, size_t length=~(size_t)0, bool rev_comp=false
        ) const {
            length = std::min(length, parent::size() - offset);
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

        template<class InIter>
        void append(InIter b, InIter e) {
            size_t num_bases = parent::size();
            parent::append(b, e);
            internal_append(num_bases, b, e);
        }

        void swap(basic_dna_string &rhs) {
            parent::swap(rhs);
            index.swap(rhs.index);
            rle.swap(rhs.rle);
        }

        void write_binary(writer &wr) const {
            parent::write_binary(wr);
            wr.write(index);
            wr.write(rle);
        }
    private:
        template<class InIter>
        void internal_append(size_t num_bases, InIter b, InIter e) {
            int prev_val = rle.empty() ? 0 : rle.back();
            while (b != e) {
                int chr = (int)*b++;
                if (!is_whitespace(chr)) {
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
        }

        // note: order matters
        IndexArrayType index;
        RleArrayType rle;
    };

    typedef basic_augmented_string<
        uint64_t,
        dna_string,
        std::vector<uint32_t>,
        std::vector<uint32_t>
    > augmented_string;

    typedef basic_augmented_string<
        uint64_t,
        mapped_dna_string,
        mapped_vector<uint32_t>,
        mapped_vector<uint32_t>
    > mapped_augmented_string;

    template <>
    static inline int get_code<augmented_string>(const augmented_string &str, size_t index) {
        return str.get_code(index);
    }
} }


#endif
