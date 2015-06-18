// Copyright Andy Thomason 2015
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

//! \file
//! \brief genetics library examples with snippets for Quickbook documentation.

// Caution: this file contains Quickbook markup as well as code
// and comments, don't change any of the special comment markups!

#include <boost/genetics/dna_string.hpp>

#include <memory>
#include <iostream>
#include <fstream>

//[simple_substring_search
void simple_substring_search()
{
    using namespace boost::genetics;

    // Make a DNA string. Only A, C, G or T allowed.
    dna_string str("ACAGTACGTAGGATACAGGTACA");

    // Search the string for a substring.
    dna_string GATACA = "GATACA";
    size_t gataca_pos = str.find(GATACA);
    if (gataca_pos != dna_string::npos) {
        std::cout << "'GATACA' found at location " << gataca_pos << "\n";
    }
}
//] [/simple_substring_search]


//[reverse_complement
void reverse_complement()
{
    using namespace boost::genetics;

    // in this example the string is reversed to AACCTTTT
    // and then complemented (T <-> A and C <-> G) to TTGGAAAA
    dna_string str("TTTTCCAA");

    std::cout << "reverse complement of " << str << " is " << rev_comp(str) << "\n";

    // this also works for std::string
    std::string stdstr("TTGGAAAA");

    std::cout << "reverse complement of " << stdstr << " is " << rev_comp(stdstr) << "\n";
}

//] [/reverse_complement]


//[brute_force_inexact_find
void inexact_find()
{
    using namespace boost::genetics;

    // make a DNA string. Only A, C, G or T allowed.
    dna_string str("ACAGAAACAGTACGTAGGATACAGGTACA");
    //                 GAAACA        GATACA

    // search the string for a substring with up to one error.
    dna_string GATACA = "GATACA";
    size_t gaaaca_pos = str.find_inexact(GATACA, 0, (size_t)-1, 1);
    if (gaaaca_pos != dna_string::npos) {
        // the first search finds GAAACA (distance one from GATACA)
        std::cout << str.substr(gaaaca_pos, 6) << " found at location " << gaaaca_pos << "\n";
        size_t gataca_pos = str.find_inexact(GATACA, gaaaca_pos+1, (size_t)-1, 1);
        if (gataca_pos != dna_string::npos) {
            // the second search finds GATACA itself.
            std::cout << str.substr(gataca_pos, 6) << " found at location " << gataca_pos << "\n";
        }
    }
}

//] [/brute_force_inexact_find]

int main()
{
    simple_substring_search();
    reverse_complement();
    inexact_find();
}
