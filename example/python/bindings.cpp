// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifdef _MSC_VER
    #pragma warning(disable : 4273)
#endif

#include <boost/genetics/augmented_string.hpp>
#include <boost/genetics/fasta.hpp>
#include <boost/python.hpp>

#include <memory>
#include <iostream>
#include <fstream>

using namespace boost::python;
using namespace boost::genetics;

/// python class representing a collection of FASTA files in a reference.
class Reference {
public:
    /// create an empty Reference
    Reference() {
    }

    /// create a Reference from a binary file.
    Reference(const std::string &binary_filename) {
        using namespace boost::interprocess;

        fm = std::make_shared<file_mapping>(binary_filename.c_str(), read_only);
        region = std::make_shared<mapped_region>(*fm, read_only);
        char *p = (char*)region->get_address();
        char *end = p + region->get_size();
        mapper m(p, end);
        fasta = std::make_shared<mapped_fasta_file>(m);
    }

    /// load from FASTA files and build index (takes several seconds)
    Reference(const list &filenames, int num_indexed_chars)
    {
        size_t len = boost::python::len(filenames);
        fasta = std::make_shared<fasta_file>();
        for (size_t i = 0; i != len; ++i) {
            fasta->append(boost::python::extract<std::string>(filenames[i]));
        }
        fasta->make_index((size_t)num_indexed_chars);
    }

    /// clean up
    ~Reference() {
    }

    /// return a list of matches to 
    list find_inexact(const std::string &str, int max_distance, int max_gap, bool is_brute_force, int max_results) const {
        std::vector<fasta_result> result;
        fasta->find_inexact(result, str, (size_t)max_distance, (size_t)max_results, (size_t)max_gap, is_brute_force);
        list py_result;
        for (size_t i = 0; i != result.size(); ++i) {
            fasta_result &r = result[i];
            const chromosome &c = fasta->find_chromosome(r.location);
            py_result.append(boost::python::make_tuple(r.location, r.location - c.start, c.name));
        }

        return py_result;
    }

    /// write a binary file for use with a map
    void write_binary_file(const std::string &filename) const {
        using namespace boost::interprocess;

        writer sizer(nullptr, nullptr);
        fasta->write_binary(sizer);
        size_t size = (size_t)sizer.get_ptr();


        // there may be a better way of creating a pre-sized writable file!
        {
          std::ofstream os(filename);
          os.seekp(size-1);
          os.write("", 1);
        }
        file_mapping fm(filename.c_str(), read_write);
        mapped_region region(fm, read_write, 0, size);
        char *p = (char*)region.get_address();
        char *end = p + region.get_size();
        writer w(p, end);
        fasta->write_binary(w);
    }

    /// write an ASCII file for human and legacy use.
    void write_ascii_file(const std::string &filename) const {
        std::ofstream os(filename);
        fasta->write_ascii(os);
    }

    std::shared_ptr<boost::genetics::fasta_file_interface> fasta;
    std::shared_ptr<boost::interprocess::file_mapping> fm;
    std::shared_ptr<boost::interprocess::mapped_region> region;
}; 

BOOST_PYTHON_MODULE(genetics)
{
    class_<Reference>("Reference", init<list, int>())
        .def(init<const std::string &>())
        .def("find_inexact", &Reference::find_inexact, (arg("str"), arg("max_distance"), arg("max_gap"), arg("is_brute_force"), arg("max_results")), "find up to max_results hits with up to max_distance errors")
        .def("write_binary_file", &Reference::write_binary_file, (arg("filename")), "write the reference to a binary file")
        .def("write_ascii_file", &Reference::write_ascii_file, (arg("filename")), "write the reference to an ascii file")
    ;
}

