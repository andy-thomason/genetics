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
        virtual ~fasta_file_interface() {}
    };
    
    template <class ChromosomeType, class StringType, class IndexType>
    class basic_fasta_file : public fasta_file_interface { 
    public:
        typedef ChromosomeType chromosome_type;
        typedef StringType string_type;
        typedef IndexType index_type;

        /// create an empty FASTA reference
        basic_fasta_file() {
        }

        /// create a FASTA reference from a text file
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

        void append(const std::string &filename) {
            using namespace boost::interprocess;
            file_mapping fm(filename.c_str(), read_only);
            mapped_region region(fm, read_only);
            const char *p = (const char*)region.get_address();
            const char *end = p + region.get_size();
            append(p, end);
        }

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
                g.end = str.size();
                str.append(b, p);
                dat.push_back(g);
                g.start = g.end;
            }
        }
        
        void swap(basic_fasta_file &rhs) {
            std::swap(str, rhs.str);
            std::swap(dat, rhs.dat);
            std::swap(idx, rhs.idx);
        }
        
        void write_binary(writer &wr) const {
            wr.write(dat);
            wr.write(str);
            wr.write(idx);
        }
        
        void find_inexact(
            std::vector<fasta_result> &result,
            const std::string &str,
            size_t max_distance,
            size_t max_results
        ) /* overload */ {
            result.resize(0);
            dna_string dna_str = str;
            for (
                index_type::iterator i = idx.find_inexact(dna_str, 0, max_distance);
                i != idx.end();
                ++i
            ) {
                fasta_result r;
                r.location = (size_t)i;
                std::string substr = str.substr(r.location, str.size());
                r.distance = 0;
                for (size_t i = 0; i != str.size(); ++i) {
                    r.distance += (str[i] != substr[i]);
                }
                if (r.distance <= max_distance) {
                    result.push_back(r);
                    if (result.size() == max_results) {
                        return;
                    }
                }
            }
        }
        
        void get_chromosomes(std::vector<chromosome> &result) {
            result = dat;
        }

        void make_index(size_t num_indexed_chars) {
            idx = index_type(str, num_indexed_chars);
        }
    private:
        // Note: order matters
        chromosome_type dat;
        string_type str;
        index_type idx;
    };
    
    typedef basic_fasta_file<
        std::vector<chromosome>,
        augmented_string,
        two_stage_index
    > fasta_file;

    typedef basic_fasta_file<
        std::vector<chromosome>,
        mapped_augmented_string,
        mapped_two_stage_index
    > mapped_fasta_file;
    
} }

#endif
