// Copyright Andy Thomason 2016
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_FM_INDEX
#define BOOST_GENETICS_FM_INDEX

#include <stdexcept>
#include <type_traits>
#include <limits>
#include <boost/genetics/dna_string.hpp>

namespace boost { namespace genetics {
    //! \brief This class is an FM Index which may be used to index a dna string.
    //! https://en.wikipedia.org/wiki/FM-index
    //! FM-indices are use in many common aligners, such as BWA and Bowtie.

    //! Like many of the container classes in this library it can be specialised
    //! into a standard (`std::vector`) version for construction and a read-only
    //! mapped (`mapped_vector`) version for high performance use.

    //! \tparam Traits Typically one of unmapped_traits or mapped_traits.
    //! These may be customised for larger address lengths, for example.

    template<class Traits>
    class basic_fm_index {
    public:
        typedef basic_dna_string<Traits> string_type;
        typedef basic_dna_string<Traits> dna_string_type;

        typedef typename Traits::SuffixArrayType addr_array_type;
        typedef typename Traits::SuffixArrayType::value_type addr_type;

        //! \brief Construct an empty suffix array.
        basic_fm_index(
        ) : string_(nullptr) {
        }
        
        //! \brief Write to a binary stream for subsequent mapping.
        void write_binary(writer &wr) const {
            wr.write(addr_);
        }

        //! \brief rvalue move operator
        basic_fm_index &operator =(basic_fm_index &&rhs) {
            string_ = rhs.string_;
            addr_ = std::move(rhs.addr_);
            return *this;
        }

        //! \brief Construct a suffix array from a persistant block of memory.
        //! This will liely be a memory mapped file.
        template <class Mapper>
        basic_fm_index(
            string_type &string,
            Mapper &map,
            typename Mapper::is_mapper *p=0
        ) :
            string_(&string),
            addr_(map),
            bwt_(map),
            inverse_sa0_((size_t)map.read64())
        {
        }
        
        //! \brief Construct an FM index for a dna string https://en.wikipedia.org/wiki/FM-index
        basic_fm_index(
            string_type &str
        ) : string_(&str) {
            str.bwt(bwt_, inverse_sa0_);
        }

        //! Swap two suffix arrays
        void swap(basic_fm_index &rhs) {
            std::swap(string_, rhs.string_);
            addr_.swap(rhs.addr_);
        }

        //! Write the suffix array to a stream for debugging.
        template <class charT, class traits>
        void write_ascii(std::basic_ostream<charT, traits>& os) const {
            os << "bwt = " << bwt_ << "\n";
            os << "inverse_sa0 = " + inverse_sa0_ << "\n";
        }
        
        //! Get the Burrows Wheeler transform
        const dna_string_type &bwt() const {
            return bwt_;
        }
        
        //! Get the Burrows Wheeler transform '$' index
        const size_t &inverse_sa0() const {
            return inverse_sa0_;
        }
        
        bool verify() const {
            if (bwt_.size() != string_->size() + 1) {
                std::cerr << "fm_index verify fail: wrong bwt size\n";
                return false;
            }
            dna_string_type ibwt;
            bwt_.ibwt(ibwt, inverse_sa0_);
            
            std::cout << *string_ << "\n";
            std::cout << bwt_ << "\n";
            std::cout << ibwt << "\n";
            
            if (ibwt.compare(0, ibwt.size(), *string_) != 0) {
                std::cerr << "fm_index verify fail: inverse bwt not same as string\n";
                return false;
            }
            return true;
        }
    private:
        // Note: order matters
        string_type *string_ = nullptr;
        addr_array_type addr_;
        dna_string_type bwt_;
        size_t inverse_sa0_ = 0;
        std::array<addr_type, 6> cumulative_index_;
    };

    //! \brief Unmapped suffix array, typically used for construction.
    typedef basic_fm_index<unmapped_traits> fm_index;

    //! \brief Mapped suffix array, typically used for high performance loading.
    typedef basic_fm_index<mapped_traits> mapped_fm_index;

    //! \brief Stream write operator
    template <class charT, class traits, class Traits>
    std::basic_ostream<charT, traits>&
    operator<<(std::basic_ostream<charT, traits>& os, const basic_fm_index<Traits>& x) {
        x.write_ascii(os);
        return os;
    }
} }


#endif
