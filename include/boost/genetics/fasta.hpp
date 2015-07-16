// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_FASTA_HPP
#define BOOST_GENETICS_FASTA_HPP

#include <type_traits>
#include <cstdint>

#include <boost/genetics/dna_string.hpp>
#include <boost/genetics/augmented_string.hpp>
#include <boost/genetics/two_stage_index.hpp>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

namespace boost { namespace genetics {
    //! Result of queries of the FASTA file
    struct fasta_result {
        size_t location;
        size_t distance;
        bool reverse_complement;
    };
   
    //! Interface to the various incarantions of the reference
    struct fasta_file_interface {
        //! find a list of results that match the string dstr with up to
        //! max_distance errors and up to max_gap bases of introns.
        //! If the flag is_brute_force is set, do the search on every base
        //! in the file using popcnt if possible.
        virtual void find_inexact(std::vector<fasta_result> &result, const std::string &dstr, search_params &params) = 0;

        //! Get chromosome data for one chromosome.
        virtual const chromosome &get_chromosome(size_t index) const = 0;

        //! Return the number of chromosomes.
        virtual size_t get_num_chromosomes() const = 0;

        //! Called after append to index the file.
        virtual void make_index(size_t num_indexed_chars) = 0;

        //! Append a FASTA file to this reference.
        virtual void append(const std::string &filename) = 0;

        //! Give chromosome data or a null entry for this linear location in the reference.
        virtual const chromosome &find_chromosome(size_t location) const = 0;

        //! Write the reference data and the index to a binary writer wr.
        virtual void write_binary(writer &wr) const = 0;

        //! Write the FASTA data in ASCII to an ostream.
        virtual void write_ascii(std::ostream &str) const = 0;

        //! Destroy the reference.
        virtual ~fasta_file_interface() {}
    };
    
    template <class Traits>
    class basic_fasta_file : public fasta_file_interface { 
    public:
        typedef typename Traits::FastaChromosomeType chromosome_type;
        typedef basic_augmented_string<Traits> string_type;
        typedef basic_two_stage_index<Traits> index_type;

        //! Create an empty FASTA reference file. Use append() to add files.
        basic_fasta_file() {
        }

        //! Create a FASTA reference from a single text file.
        template <class T = Traits::unmapped>
        basic_fasta_file(const std::string &filename) {
            append(filename);
        }

        //! Use a mapper to instantly load from a mapped file.
        template <class T = Traits::mapped>
        basic_fasta_file(mapper &map) :
            chromosomes(map),
            str(map),
            idx(str, map)
        {
        }

        //! Move from another reference of the same type.
        basic_fasta_file &operator=(basic_fasta_file &&rhs) {
            chromosomes = std::move(rhs.chromosomes);
            str = std::move(rhs.str);
            idx = std::move(rhs.idx);
            return *this;
        }

        virtual ~basic_fasta_file() {
        }

        //! Append a single FASTA file to this file.
        void append(const std::string &filename) {
            using namespace boost::interprocess;
            file_mapping fm(filename.c_str(), read_only);
            mapped_region region(fm, read_only);
            const char *p = (const char*)region.get_address();
            const char *end = p + region.get_size();
            append(p, end);
        }

        //! Append a randomly generated chromosome "chr" of "size" bases.
        void append_random(const std::string &chr, size_t size, std::uint32_t seed) {
            std::string random;
            random.resize(size + chr.size() + 2);
            auto p = random.begin(), e = random.end();
            *p++ = '@';
            for (auto c : chr) if (p != e) *p++ = c;
            if (p != e) *p++ = '\n';

            auto b = p;
            std::int32_t r = seed;
            while (p != b + size/16) {
                *p++ = 'N';
            }
            while (p != b + size*15/16) {
                r = ( ( r >> 31 ) & 0xa67b9c35 ) ^ ( r << 1 );
                r = ( ( r >> 31 ) & 0xcb731943 ) ^ ( r << 1 );
                if (p != e) *p++ = (r & 0xfff0000) == 0 ? 'N' : "ACGT"[r & 3];
            }
            while (p != e) {
                *p++ = 'N';
            }
            append(random.data(), random.data() + size);
        }

        //! Append the image of a FASTA file to this file.
        void append(const char *p, const char *end) {
            chromosome g;
            while (p != end) {
                g.start = str.size();

                //if (*p != '>') throw(std::runtime_error("bad fasta"));
                const char *b = p;
                while (p != end && *p != '\n' && *p != '\r') ++p;
                size_t size = std::min(sizeof(g.info)-1, (size_t)(p-(b+1)));
                memcpy(g.info, b+1, size);

                g.info[size] = 0;
                size_t i = 0;
                for (; i != size && i != sizeof(g.name)-1 && g.info[i] != ' '; ++i) {
                    g.name[i] = g.info[i];
                }
                g.name[i] = 0;

                b = p;
                while (p != end && *p != '>') ++p;
                str.append(b, p);
                g.end = str.size();

                chromosomes.push_back(g);
            }
        }
        
        //! Swap references.
        void swap(basic_fasta_file &rhs) {
            std::swap(str, rhs.str);
            std::swap(chromosomes, rhs.chromosomes);
            std::swap(idx, rhs.idx);
        }
        
        //! copy the bytes in this file to an image.
        void write_binary(writer &wr) const {
            wr.write(chromosomes);
            str.write_binary(wr);
            idx.write_binary(wr);
        }
        
        //! Write as an ASCII FASTA file.
        void write_ascii(std::ostream &os) const {
            std::string line;
            line.reserve(256);
            for (size_t i = 0; i != chromosomes.size(); ++i) {
                const chromosome &c = chromosomes[i];
                size_t info_size = strlen(c.info);
                line.resize(info_size + 2);
                size_t j = 0;
                line[j++] = '>';
                for (const char *p = c.info; *p; ++p) line[j++] = *p;
                line[j++] = '\n';
                os << line;
                for (size_t b = c.start; b < c.end; ) {
                    size_t e = std::min(b + 60, c.end);
                    line.resize(e - b + 1);
                    for (size_t i = 0; i != (e-b); ++i) {
                        line[i] = str[b + i];
                    }
                    line[e-b] = '\n';
                    os << line;
                    b = e;
                }
            }

        }
        
        //! Search the FASTA file for strings with some allowable errors.
        void find_inexact(std::vector<fasta_result> &result, const std::string &dstr, search_params &params) {
            result.resize(0);
            int max_pass = params.search_rev_comp ? 2 : 1;
            for (int pass = 0; pass != max_pass; ++pass) {
                bool reverse_complement = pass == 1;
                std::string search_str = reverse_complement ? rev_comp(dstr) : dstr;
                for (
                    auto i = idx.find_inexact(search_str, 0, params);
                    i != idx.end();
                    ++i
                ) {
                    fasta_result r;
                    r.location = (size_t)i;
                    r.reverse_complement = reverse_complement;
                    r.distance = i.distance();
                    if (r.distance <= params.max_distance) {
                        result.push_back(r);
                        if (result.size() == params.max_results) {
                            return;
                        }
                    }
                }
            }
        }
        
        //! Get the chromosomes in this file.
        const chromosome &get_chromosome(size_t index) const {
            return chromosomes[index];
        }

        size_t get_num_chromosomes() const {
            return chromosomes.size();
        }

        const string_type &get_string() const {
            return str;
        }

        //! Must be called after appending FASTA data.
        void make_index(size_t num_indexed_chars) {
            idx = index_type(str, num_indexed_chars);
        }

        //! Number of base pairs in this reference.
        size_t size() const {
            return str.size();
        }

        //! Find a chomosome for a location.
        const chromosome &find_chromosome(size_t location) const {
            const chromosome *end = chromosomes.data() + chromosomes.size();
            const chromosome *i = std::lower_bound(chromosomes.data(), end, location);
            if (i != end && location >= i[0].start && location < i[0].end) {
                return i[0];
            } else {
                return null_chr;
            }
        }
    private:
        // Note: for these data members, order matters because the map constructor requires this.

        //! Empty chromosome for out-of-range queries.
        chromosome null_chr;

        //! Vector of chromosomes.
        chromosome_type chromosomes;

        //! Usually an augmented_string type to store bases.
        string_type str;

        //! Index on str.
        index_type idx;
    };

    //! This container is a writable type for conversion from ASCII files.
    typedef basic_fasta_file<unmapped_traits> fasta_file;

    //! This container is read only for mapped files.
    typedef basic_fasta_file<mapped_traits> mapped_fasta_file;

} }

#endif
