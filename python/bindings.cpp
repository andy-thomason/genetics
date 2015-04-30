
#include <boost/genetics/genetics.hpp>
#include <boost/python.hpp>

using namespace boost::python;

class Fasta {
    typedef boost::genetics::dna_string dna_string; 
    typedef boost::genetics::augmented_string augmented_string;
    boost::genetics::fasta_file fasta; 
public:
    Fasta(const std::string &filename) : fasta(filename)
    {
    }

    ~Fasta() {
    }

    list find_all(const std::string &key, int max_distance) const {
        list result;
        return result;
    }

    list get_chromosomes() const {
        list result;
        
        return result;
    }
}; 

BOOST_PYTHON_MODULE(genetics)
{
    class_<Fasta>("Fasta")
        .def("find_all", &Fasta::find_all)
        .def("get_chromosomes", &Fasta::get_chromosomes)
    ;
}

