// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#pragma warning(disable : 4273)
#include <boost/genetics/augmented_string.hpp>
#include <boost/genetics/fasta.hpp>
#include <boost/python.hpp>

using namespace boost::python;
using namespace boost::genetics;

class Fasta {
public:
    Fasta() {
    }

    /// load from single FASTA file
    /*Fasta(const std::string &filename) {
        fasta = new fasta_file();
        fasta->append(filename);
        int num_indexed_chars = (56-lzcnt(fasta->size())/2;
        num_indexed_chars = std::max(6, std::min(12, num_indexed_chars));
        fasta->make_index((size_t)num_indexed_chars);
    }*/

    ~Fasta() {
        delete fasta;
    }

    /// load from FASTA files (takes several seconds)
    Fasta(const list &filenames, int num_indexed_chars)
    {
        size_t len = boost::python::len(filenames);
        fasta = new fasta_file();
        for (size_t i = 0; i != len; ++i) {
            fasta->append(boost::python::extract<std::string>(filenames[i]));
        }
        fasta->make_index((size_t)num_indexed_chars);
    }

    /// map instantly from a binary file
    void map(const std::string &binary_filename)
    {
        fasta = new mapped_fasta_file(binary_filename);
    }

    /// return a list of matches to 
    list find_inexact(const std::string &str, int max_distance, int max_results) const {
        std::vector<fasta_result> result;
        fasta->find_inexact(result, str, (size_t)max_distance, (size_t)max_results);
        list py_result;
        for (size_t i = 0; i != result.size(); ++i) {
            fasta_result &r = result[i];
            const chromosome &c = fasta->find_chromosome(r.location);
            py_result.append(boost::python::make_tuple(r.location, r.location - c.start, c.name));
        }

        return py_result;
    }

    list get_chromosomes() const {
        list result;

        return result;
    }

    fasta_file_interface *fasta;
}; 

BOOST_PYTHON_MODULE(genetics)
{
    class_<Fasta>("Fasta", init<list, int>())
        .def("find_inexact", &Fasta::find_inexact)
        //.property("chromosomes", &Fasta::get_chromosomes)
    ;
}

