// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include<string.h>
#include <boost/genetics/fasta.hpp>

#define BOOST_TEST_MODULE genetics
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( fasta_test )
{
    using namespace boost::genetics;
    
    fasta_file f("ensembl_chr21.fa");
    f.append("ensembl_chr21.fa");
    f.make_index(4);
    
    
    BOOST_CHECK(f.get_num_chromosomes() == 4);

    BOOST_CHECK(!strcmp(f.get_chromosome(0).name, "21"));
    BOOST_CHECK(!strcmp(f.get_chromosome(0).info, "21 fake chromosome 21 from ENSEMBL data"));

    BOOST_CHECK(!strcmp(f.get_chromosome(1).name, "22"));
    BOOST_CHECK(!strcmp(f.get_chromosome(1).info, "22 fake chromosome 22"));

    BOOST_CHECK(!strcmp(f.get_chromosome(2).name, "21"));
    BOOST_CHECK(!strcmp(f.get_chromosome(2).info, "21 fake chromosome 21 from ENSEMBL data"));

    BOOST_CHECK(!strcmp(f.get_chromosome(3).name, "22"));
    BOOST_CHECK(!strcmp(f.get_chromosome(3).info, "22 fake chromosome 22"));

    const fasta_file::string_type &str = f.get_string();
    BOOST_CHECK(str.substr((5-2)*60, 60) == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
    BOOST_CHECK(str.substr((6-2)*60, 60) == "GATCCACCCGCCTTGGCCTCCTAAAGTGCTGGGATTACAGGTGTTAGCCACCACGTCCAG");
}

