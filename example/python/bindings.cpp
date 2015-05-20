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

    /// load from FASTA files (takes several seconds)
    void load(const list &filenames, int num_indexed_chars)
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

    ~Fasta() {
    }

    /// return a list of matches to 
    list find_inexact(const std::string &str, int max_distance, int max_results) const {
        std::vector<fasta_result> result;
        fasta->find_inexact(result, str, (size_t)max_distance, (size_t)max_results);
        list py_result;
        for (size_t i = 0; i != result.size(); ++i) {
            py_result.append(result[i].location);
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
    class_<Fasta>("Fasta")
        .def("load", &Fasta::load)
        .def("map", &Fasta::map)
        .def("find_inexact", &Fasta::find_inexact)
        //.property("chromosomes", &Fasta::get_chromosomes)
    ;
}

