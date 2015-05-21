// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_FASTA_HPP
#define BOOST_GENETICS_FASTA_HPP

#include <boost/genetics/dna_string.hpp>
#include <boost/genetics/augmented_string.hpp>
#include <boost/genetics/two_stage_index.hpp>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

namespace boost { namespace genetics {
    struct chromosome {
        std::string name;
        std::string info;
        size_t start;
        size_t end;

        bool operator <(size_t pos) const {
            return end < pos;
        }

        chromosome() : start(0), end(0) {
        }
    };
    
    struct fasta_result {
        size_t location;
        size_t distance;
    };
    
    struct fasta_file_interface {
        virtual void find_inexact(std::vector<fasta_result> &result, const std::string &str, size_t max_distance, size_t max_results) = 0;
        virtual void get_chromosomes(std::vector<chromosome> &result) = 0;
        virtual void make_index(size_t num_indexed_chars) = 0;
        virtual void append(const std::string &filename) = 0;
        virtual const chromosome &find_chromosome(size_t location) const = 0;
        virtual ~fasta_file_interface() {}
    };
    
    template <class ChromosomeType, class StringType, class IndexType>
    class basic_fasta_file : public fasta_file_interface { 
    public:
        typedef ChromosomeType chromosome_type;
        typedef StringType string_type;
        typedef IndexType index_type;

        /// create an empty FASTA reference file. Use append() to add files.
        basic_fasta_file() {
        }

        /// create a FASTA reference from a single text file
        basic_fasta_file(const std::string &filename) {
            append(filename);
        }
        
        basic_fasta_file &operator=(basic_fasta_file &&rhs) {
            dat = std::move(rhs.dat);
            str = std::move(rhs.str);
            idx = std::move(rhs.idx);
            return *this;
        }

        basic_fasta_file(mapper &map) :
            dat(map),
            str(map),
            idx(map)
        {
        }
        
        virtual ~basic_fasta_file() {
        }

        /// append a single FASTA file to this file.
        void append(const std::string &filename) {
            using namespace boost::interprocess;
            file_mapping fm(filename.c_str(), read_only);
            mapped_region region(fm, read_only);
            const char *p = (const char*)region.get_address();
            const char *end = p + region.get_size();
            append(p, end);
        }

        /// append the image of a FASTA file to this file.
        void append(const char *p, const char *end) {
            chromosome g;
            g.start = 0;
            size_t start = 0;
            while (p != end) {
                //if (*p != '>') throw(std::exception("bad fasta"));
                const char *b = p;
                while (p != end && *p != '\n' && *p != '\r') ++p;
                g.info.assign(b+1, p);
                size_t sp = g.info.find(" ");
                if (sp == std::string::npos) {
                    g.name = g.info;
                } else {
                    g.name.assign(g.info.begin(), g.info.begin() + sp);
                }

                b = p;
                while (p != end && *p != '>') ++p;
                str.append(b, p);
                g.end = str.size();
                dat.push_back(g);
                g.start = g.end;
            }
        }
        
        /// swap FASTA files.
        void swap(basic_fasta_file &rhs) {
            std::swap(str, rhs.str);
            std::swap(dat, rhs.dat);
            std::swap(idx, rhs.idx);
        }
        
        /// copy the bytes in this file to an image.
        void write_binary(writer &wr) const {
            wr.write(dat);
            wr.write(str);
            wr.write(idx);
        }
        
        /// search the FASTA file for strings with some allowable errors.
        void find_inexact(
            std::vector<fasta_result> &result,
            const std::string &dstr,
            size_t max_distance,
            size_t max_results
        ) /* overload */ {
            result.resize(0);
            dna_string dna_str = dstr;
            for (
                index_type::iterator i = idx.find_inexact(dna_str, 0, max_distance);
                i != idx.end();
                ++i
            ) {
                fasta_result r;
                r.location = (size_t)i;
                std::string substr = str.substr(r.location, dstr.size());
                r.distance = 0;
                for (size_t i = 0; i != dstr.size(); ++i) {
                    r.distance += (dstr[i] != substr[i]);
                }
                if (r.distance <= max_distance) {
                    result.push_back(r);
                    if (result.size() == max_results) {
                        return;
                    }
                }
            }
        }
        
        /// get the chromosomes in this file.
        void get_chromosomes(std::vector<chromosome> &result) {
            result = dat;
        }

        /// must be called after appending FASTA data.
        void make_index(size_t num_indexed_chars) {
            idx = index_type(str, num_indexed_chars);
        }

        /// number of base pairs in this file.
        size_t size() const {
            return str.size();
        }

        /// find a chomosome for a location.
        const chromosome &find_chromosome(size_t location) const {
            const chromosome *end = dat.data() + dat.size();
            const chromosome *i = std::lower_bound(dat.data(), end, location);
            if (i != end && location >= i[0].start && location < i[0].end) {
                return i[0];
            } else {
                return null_chr;
            }
        }
    private:
        // Note: for these data members, order matters because the map constructor requires this.

        /// empty chromosome for out-of-range queries.
        chromosome null_chr;

        /// vector of chromosomes.
        chromosome_type dat;

        /// usually an augmented_string type to store bases.
        string_type str;

        /// index on str
        index_type idx;
    };
    
    /// this container is a writable type for conversion from ASCII files.
    typedef basic_fasta_file<
        std::vector<chromosome>,
        augmented_string,
        two_stage_index
    > fasta_file;

    /// this container is read only for mapped files.
    typedef basic_fasta_file<
        std::vector<chromosome>,
        mapped_augmented_string,
        mapped_two_stage_index
    > mapped_fasta_file;
    
} }

#endif
