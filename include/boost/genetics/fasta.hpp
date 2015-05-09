
#ifndef BOOST_GENETICS_FASTA_HPP
#define BOOST_GENETICS_FASTA_HPP

#include <boost/genetics/dna_string.hpp>
#include <boost/genetics/augmented_string.hpp>
#include <boost/genetics/two_stage_index.hpp>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

namespace boost { namespace genetics {
    struct genome_data {
        std::string name;
        std::string info;
        size_t start;
        size_t end;
    };

    template <class DataType, class StringType, class IndexType>
    class basic_fasta_file { 
    public:
        typedef DataType data_type;
        typedef StringType string_type;
        typedef IndexType index_type;

        /// create an empty FASTA reference
        basic_fasta_file() {
        }

        /// create a FASTA reference from a 
        basic_fasta_file(const std::string &filename) {
            append(filename);
        }
        
        basic_fasta_file &operator=(basic_fasta_file &&rhs) {
            str = std::move(rhs.str);
            dat = std::move(rhs.dat);
            idx = std::move(rhs.idx);
        }

        const string_type &string() const {
            return str;
        } 

        const data_type &data() const {
            return dat;
        }
        
        const index_type &index() const {
            return idx;
        }
        
        void make_index() {
            idx = two_stage_index(str, 6);
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
            genome_data g;
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

    private:
        string_type str;
        data_type dat;
        index_type idx;
    };
    
    typedef basic_fasta_file<
        std::vector<genome_data>,
        augmented_string,
        two_stage_index
    > fasta_file;

} }

#endif