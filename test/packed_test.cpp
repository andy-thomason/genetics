// Copyright Andy Thomason 2015
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <sstream>
#include <utility>

#define BOOST_TEST_MODULE Genetics
#include <boost/test/unit_test.hpp>

#include <boost/genetics/utils.hpp>
#include <boost/genetics/dna_string.hpp>
#include <boost/genetics/augmented_string.hpp>
#include <boost/genetics/two_stage_index.hpp>
#include <boost/genetics/fm_index.hpp>

static const char chr1[] =
  "TGTGATTAATGCCTGAGACTGTGTGAAGTAAGAGATGGATCAGAGGCCGGGCGCGGGGGC"
  "TCGCGCCTGTCATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGA"
  "TCGAGACCATCCTGGCTAACACGGGGAAACCCCGTCTCCACTAAAAATACAAAAAGTTAG" // 120
  "CCGGGCGCGGTGGCGGGCGCCTGCGGTCCCAGCTGCTGGGGAGGCCGAGGCGGGAGCATG"
  "GCGGGAACCGGGAGGCGGAGCCTGCAGTGAGCCGAGATGGCGCCACCGCACTCCAGCCTG"
  "GGCGACCCAGCGAGACTCCGCCTCAAAAAAAAAAAAAGAAGATTGATCAGAGAGTACCTC"
  "CCCTAAGGGTACATGCAGATAAATACAGTTAAGGCGATTAACATTTCAAATACGGTGACT"
  "GTTTCTTACGTGGACGACGTTGTGTTGAACATGGGTGAGTAAGACTGAAGCAGCCGTAAT"
  "TACTGCACGATGCGCATGGTAAAGAAGCACTCCGTTAGGGAAATTATATTCTTTGCCCCT"
  "CTAATCCTTCACTCCACCTGCCATATTCCCACATGATTTTTTTCTTTGCTGTTCTTGTCT"
  "AATTGTTATTAATAATTAATAAATAACTTATGATCTAATTGTTATTAATAATAACTTATC"
  "ATCACATGATTTATTAATAAATTAATAAATAACTTATTATCACCGCATTTCCCCAATTCA"
  "TTTATCTTTCTTTCATTTTCTCTCTTTGTGTGTTTTCTGTCTTCATATTTCAGCACTTGC"
  "CACATATTTCCCACAAAATCATTTATGGTCAAACAACACTTCAACGTGTAGCATTTGTAT"
  "TTCTCAATTCTTCCTCACTTTCTTCCTTCAGAATACTAAAGCTTCTTCTCTACTGACTGA"
  "GTCAATGGCCAATGGATAGAGTAAATAATTCTGCGGTATCTAAATTTGTATTGATTGGAC"
  "TTTCAAGCTCTTGGGAGATGCATCTTTTTCTTTTTTGGTTCTTCTCTGTGTTCTACATGG"
  "GAATTATCCTGGAAAATCTCTTCATTGTGTTCACAGTAATTATTGACTCTCATTTAAATT"
  "CCCCAGGTACTGCCTACTGGCCAACATTTATCTTCTTGATCTGGGTCTTCTCCTACAGTT"
  "CTGACTTTTTCACTAACTGCAGCATCATTTCTTTTCCAAGATGCATCATACAGATATTTT"
  "TCGAGACCATCCTGGCTAACACGGGGAAACCCCGTCTCCACTAAAAATACAAAAAGTTAG" // 1200
  "TCGAGACCATCCTGGCTAACACGCGGAAACCCCGTCTCCACTAAAAATACAAAAAGTTAG" // 1260
  "TCGAGACCATCCTGGCTAACACGCGGAAACCCCGTCTCCACTAATAATACAAAAAGTTAG" // 1320
  //"CCCCAGGTACTGCCTACTGGCCAACATTTATCTTCTTGATCTGGGTCTTCTCCTACAGTT"
;

BOOST_AUTO_TEST_CASE( dna_string_window )
{
    using namespace boost::genetics;

    const char test_string[] = "ACTGTGAC" "TTAACCGG" "TCTCAGAG" "ACGATTGG" "CCAATTGA";
    dna_string str(test_string);

    {
        dna_string::word_type w0 = str.window(0);
        dna_string::word_type w1 = str.window(32);
        dna_string::word_type w2 = str.window(16);
        dna_string::word_type w3 = str.window(-4);
        BOOST_CHECK(w2 == ((w0 << 32) | (w1 >> 32)));
        BOOST_CHECK(w3 == w0 >> 8);
    }
}

template <class Type>
void acgt_container_tests() {
    {
        std::string str("tTAACcGGTTACgtacgtacGTACGTACGTACGTACGT");
        Type b( str );
        BOOST_CHECK( b == "TTAACCGGTTACGTACGTACGTACGTACGTACGTACGT" );
    }

    {
        Type b;
        std::stringstream ss("TGCAACGTttttXXXX\nAACCGGTT");
        ss >> b;
        //BOOST_TEST_MESSAGE( b );
        BOOST_CHECK( b == "TGCAACGTTTTTAAAA" );
        ss >> b;
        BOOST_TEST_MESSAGE( b );
        BOOST_CHECK( b == "AACCGGTT" );
    }

}

template <class Type>
void generic_string_tests() {
    {
        Type b("ACGT");
        BOOST_TEST_MESSAGE( b );
        BOOST_CHECK( std::string(b) == "ACGT" );
        b.append("TTGA");
        BOOST_TEST_MESSAGE( b );
        BOOST_CHECK( std::string(b) == "ACGTTTGA" );
        b.append("ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG");
        BOOST_TEST_MESSAGE( b );
        BOOST_CHECK( std::string(b) == "ACGTTTGA" "ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG" );
        //BOOST_CHECK( b[-1] == 'N' );
        BOOST_CHECK( b[0] == 'A' );
        BOOST_CHECK( b[31] == 'G' );
        BOOST_CHECK( b[32] == 'A' );
        BOOST_CHECK( b[33] == 'C' );
        //BOOST_CHECK( b[999] == 'N' );
    }

    {
        Type b("TGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACAT");
        BOOST_TEST_MESSAGE( b );
        BOOST_CHECK( b == "TGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACAT");
    }

    {
        Type b("ACGTACGTACGT");
        char buf[64] = {0};
        std::stringstream ss(buf);
        ss << b;
        BOOST_TEST_MESSAGE( std::string(ss.str()) );
        BOOST_CHECK( std::string(ss.str()) == "ACGTACGTACGT" );
    }
}

template <class Type>
void find_test() {
    {
        Type a(chr1);
        Type key1("TCGCGCCTGTCATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGA");
        //std::cout << key1 << "\n";
        BOOST_CHECK( a.compare(60, key1.size(), key1) == 0);
    }
    {
        Type a(chr1);
        Type key1("T");
        BOOST_CHECK( a.compare(60, key1.size(), key1) == 0);
    }
    {
        Type a(chr1);
        Type key1("TCGCGCCTGTCATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGA");
        BOOST_CHECK( a.find(key1, 0) == (size_t)60);
    }
    {
        Type a("ACGTACGTACGTACGTACGTACGTACGTACGC" "ATGTACGTACGTACGTACGTACGTACGTACGT");
        Type key1("CA");
        BOOST_CHECK( a.find(key1, 0) == (size_t)31);
    }
    {
        Type a("ACGTACGTACGTACGTACGTACGTACGTACGC" "ATGTACGTACGTACGTACGTACGTACGTACGT");
        Type key1("CAT");
        BOOST_CHECK( a.find(key1, 0) == (size_t)31);
    }
    {
        Type a("ACGTACGTACGTACGTACGTACGTACGTACCA" "TTGTACGTACGTACGTACGTACGTACGTACGT");
        Type key1("CAT");
        BOOST_CHECK( a.find(key1, 0) == (size_t)30);
    }
}

BOOST_AUTO_TEST_CASE( generic_find_tests )
{
    using namespace boost::genetics;

    find_test<std::string>();
    find_test<dna_string>();
}

BOOST_AUTO_TEST_CASE( generic_string_test )
{
    using namespace boost::genetics;

    generic_string_tests<std::string>();
    generic_string_tests<dna_string>();
}

BOOST_AUTO_TEST_CASE( find_inexact )
{
    using namespace boost::genetics;

    {
        // wildcard searches
        dna_string a("ACGTACGTACGTACGTACGTACGTACGTACCA" "TTGTACGTACGTACGTACGTACGTACGTACGT");
        dna_string key3("ACA");
        BOOST_CHECK( a.find_inexact(key3, 0, 64, 1) == 0); // matches ACG
        BOOST_CHECK( a.find_inexact(key3, 1, 64, 1) == 4); // matches ACG
        dna_string key4("ACCGTTGT");
        BOOST_CHECK( a.find_inexact(key4, 0, 64, 1) == 28); // matches "ACCATTGT"
    }
    {
        // check boundaries and limits
        dna_string a("ACGTACGTACGTACGTACGTACGTACGTACCA" "TTGTACGTACGTACGTACGTACGTACGTACGT");
        dna_string key1("CAT");
        BOOST_CHECK( a.find_inexact(key1, 0, 31) == dna_string::npos);
        BOOST_CHECK( a.find_inexact(key1, 0, 32) == dna_string::npos);
        BOOST_CHECK( a.find_inexact(key1, 0, 33) == (size_t)30);
        BOOST_CHECK( a.find_inexact(key1, 30, 33) == (size_t)30);
        BOOST_CHECK( a.find_inexact(key1, 31, 33) == dna_string::npos);
        dna_string key2("TTG");
        BOOST_CHECK( a.find_inexact(key2, 0, 33) == dna_string::npos);
        BOOST_CHECK( a.find_inexact(key2, 0, 34) == dna_string::npos);
        BOOST_CHECK( a.find_inexact(key2, 0, 35) == (size_t)32);
    }
}

BOOST_AUTO_TEST_CASE( bandl_tests )
{
    using namespace boost::genetics;

    {
        static const char test_string[] = "ACGTACGTXXXXXACGT";
        augmented_string b(test_string);
        BOOST_CHECK(std::string(b) == test_string);
    }

    {
        static const char test_string[] = "NNNNACTGACTGACTG";
        augmented_string b(test_string);
        BOOST_CHECK(std::string(b) == test_string);
    }

    {
        std::string test_string(0xffff, 'N');
        test_string.append(256, 'A');
        augmented_string b(test_string);
        BOOST_CHECK(std::string(b) == test_string);
    }

    {
        std::string test_string(0x10000, 'N');
        test_string.append(256, 'A');
        augmented_string b(test_string);
        BOOST_CHECK(std::string(b) == test_string);
    }

}

BOOST_AUTO_TEST_CASE( rev_comp_test )
{
    using namespace boost::genetics;

    {
        dna_string b("ACGTACGTACGTACGTACGTACGTACGTACGT");
        BOOST_TEST_MESSAGE( rev_comp(b) );
        BOOST_CHECK( rev_comp(b) == "ACGTACGTACGTACGTACGTACGTACGTACGT" );
    }

    {
        dna_string b("TGCAACACACA");
        BOOST_TEST_MESSAGE( rev_comp(b) );
        BOOST_CHECK( rev_comp(b) == "TGTGTGTTGCA" );
    }

    {
        std::string b("TGNAACACACA");
        BOOST_CHECK( rev_comp(b) == "TGTGTGTTNCA" );
    }
}



BOOST_AUTO_TEST_CASE( substr )
{
    using namespace boost::genetics;

    const char test_string[] = "ACGTTTGA" "ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG";

    {
        dna_string b("ACGT");
        b.resize(3);
        b.resize(4);
        BOOST_TEST_MESSAGE( b );
        BOOST_CHECK( b == "ACGA" );
    }

    {
        dna_string b(test_string);
        dna_string c;
        c = b.substr(8, 8);
        BOOST_TEST_MESSAGE( c );
        BOOST_CHECK( c == "ACCGACCG" );
        c = b.substr(8, 40);
        BOOST_TEST_MESSAGE( c );
        dna_string d(test_string + 8);
        BOOST_TEST_MESSAGE( dna_string(test_string + 8) );
        BOOST_CHECK( c == test_string + 8);
    }

    {
        std::string b(test_string);
        std::string c;
        c = b.substr(8, 8);
        BOOST_TEST_MESSAGE( c );
        BOOST_CHECK( c == "ACCGACCG" );
        c = b.substr(8, 40);
        BOOST_TEST_MESSAGE( c );
        BOOST_CHECK( c == test_string + 8);
    }

}

BOOST_AUTO_TEST_CASE( two_stage_index_test )
{
    using namespace boost::genetics;

    augmented_string as(chr1);
    two_stage_index tsi(as, 4);
    //BOOST_TEST_MESSAGE(tsi);

    search_params params;
    search_stats stats;

    augmented_string key1("TCGAGACCATCCTGGCTAACACGGGGAAACCCCGTCTCCACTAAAAATACAAAAAGTTAG");

    {
        two_stage_index::iterator i = tsi.find_inexact(key1, 0, params, stats);
        BOOST_CHECK(i == 120);
        ++i;
        BOOST_CHECK(i == 1200);
        ++i;
        BOOST_CHECK(i == augmented_string::npos);
    }
    {
        params.max_distance = 1;
        two_stage_index::iterator i = tsi.find_inexact(key1, 121, params, stats);
        BOOST_CHECK(i == 1200);
        ++i;
        BOOST_CHECK(i == 1260);
        ++i;
        BOOST_CHECK(i == augmented_string::npos);
    }
    {
        params.max_distance = 2;
        two_stage_index::iterator i = tsi.find_inexact(key1, 1201, params, stats);
        BOOST_CHECK(i == 1260);
        ++i;
        BOOST_CHECK(i == 1320);
        ++i;
        BOOST_CHECK(i == augmented_string::npos);
    }
}

BOOST_AUTO_TEST_CASE( mapped_container_test )
{
    using namespace boost::genetics;

    boost::genetics::uint64_t buf[64];
    char *ptr = 0;

    {
        dna_string dna("TTTTTTTTGGGGGGGGCCCCCCCCAAAA");
        writer wr((char*)buf, (char*)(buf + 4));
        dna.write_binary(wr);
        ptr = wr.get_ptr();
        BOOST_CHECK(wr.is_end());
        BOOST_CHECK(buf[0] == 28);
        BOOST_CHECK(buf[1] == 8);
        BOOST_CHECK(buf[2] == 1);
        BOOST_CHECK(buf[3] == 0xffffaaaa55550000);
    }

    {
        mapper map((const char*)buf, ptr);
        mapped_dna_string dna(map);
        BOOST_CHECK(map.is_end());
        BOOST_CHECK(std::string(dna) == "TTTTTTTTGGGGGGGGCCCCCCCCAAAA");
    }

    {
        augmented_string dna("TTTTTTNNNNGGGGGGCCCCCCCCAAAA");
        writer wr((char*)buf, (char*)(buf + 32));
        dna.write_binary(wr);
        ptr = wr.get_ptr();
    }

    {
        mapper map((const char*)buf, ptr);
        mapped_augmented_string dna(map);
        BOOST_CHECK(map.is_end());
        BOOST_CHECK(std::string(dna) == "TTTTTTNNNNGGGGGGCCCCCCCCAAAA");
    }

    {
        augmented_string dna("TTTTTTNNNNGGGGGGCCCCCCCCAAAA");
        two_stage_index tsi(dna, 2);
        writer wr((char*)buf, (char*)(buf + 32));
        dna.write_binary(wr);
        tsi.write_binary(wr);
        ptr = wr.get_ptr();
    }

    {
        search_params params;
        search_stats stats;
        mapper map((const char*)buf, ptr);
        mapped_augmented_string dna(map);
        mapped_two_stage_index tsi(dna, map);
        BOOST_CHECK(map.is_end());
        BOOST_CHECK(std::string(dna) == "TTTTTTNNNNGGGGGGCCCCCCCCAAAA");
        augmented_string key("GG");
        mapped_two_stage_index::iterator i = tsi.find_inexact(key, 0, params, stats);
        BOOST_CHECK(i == 10);
        i++;
        BOOST_CHECK(i == 11);
    }
}

BOOST_AUTO_TEST_CASE( fm_index_test )
{
    using namespace boost::genetics;
    {
        dna_string dna(chr1);
        fm_index fm(dna);
        BOOST_CHECK(fm.verify());
    }
}

BOOST_AUTO_TEST_CASE( bwt_test )
{
    using namespace boost::genetics;

    dna_string bwt;
    dna_string ibwt;
    size_t inverse_sa0 = 0;

    {
        dna_string dna("ACGT");
        dna.bwt(bwt, inverse_sa0);
        bwt.ibwt(ibwt, inverse_sa0);
        BOOST_CHECK(bwt == dna_string("T$ACG"));
        BOOST_CHECK(inverse_sa0 == 1);
        BOOST_CHECK(ibwt == dna);
    }

    {
        dna_string dna("TGCA");
        dna.bwt(bwt, inverse_sa0);
        bwt.ibwt(ibwt, inverse_sa0);
        BOOST_CHECK(bwt == dna_string("ACGT$"));
        BOOST_CHECK(inverse_sa0 == 4);
        BOOST_CHECK(ibwt == dna);
    }

    {
        dna_string dna("TCGCGCCTGTCATCCCTTCGCGCCTGTCATCCC");
        dna.bwt(bwt, inverse_sa0);
        bwt.ibwt(ibwt, inverse_sa0);
        BOOST_CHECK(bwt == dna_string("CCCCTTCTTGGCGGTTCCCCCCCTTGGAAT$CCC"));
        BOOST_CHECK(inverse_sa0 == 30);
        BOOST_CHECK(ibwt == dna);
    }

    {
        dna_string dna("AAAAAAAAAAAAAAAAAAAA");
        dna.bwt(bwt, inverse_sa0);
        bwt.ibwt(ibwt, inverse_sa0);
        BOOST_CHECK(bwt == dna_string("AAAAAAAAAAAAAAAAAAAA$"));
        BOOST_CHECK(inverse_sa0 == 20);
        BOOST_CHECK(ibwt == dna);
    }

    {
        dna_string dna;
        dna.bwt(bwt, inverse_sa0);
        bwt.ibwt(ibwt, inverse_sa0);
        BOOST_CHECK(bwt == dna_string("$"));
        BOOST_CHECK(inverse_sa0 == 0);
        BOOST_CHECK(ibwt == dna);
    }

    {
        dna_string dna(chr1);
        //auto t0 = std::chrono::high_resolution_clock::now();
        dna.bwt(bwt, inverse_sa0);
        //auto t1 = std::chrono::high_resolution_clock::now();
        //std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() << "us\n";
        bwt.ibwt(ibwt, inverse_sa0);
        //auto t2 = std::chrono::high_resolution_clock::now();
        //std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us\n";
        BOOST_CHECK(ibwt == dna);
    }
}

BOOST_AUTO_TEST_CASE( occurance_test )
{
    using namespace boost::genetics;
    
    auto count = [](const char *b, const char *e) {
        std::array<size_t, 4> res;
        res[0] = res[1] = res[2] = res[3] = 0;
        while (b != e) {
            if (*b == 'C') {
                res[1]++;
            } else if (*b == 'G') {
                res[2]++;
            } else if (*b == 'T') {
                res[3]++;
            } else {
                res[0]++;
            }
            ++b;
        }
        return res;
    };
    
    {
        const char *str = "AAAACCCGGGGGTTTTTT";
        dna_string dna(str);
        auto dna_occ = dna.occurance(0, dna.size());
        auto str_occ = count(str, str+dna.size());
        BOOST_CHECK(dna_occ == str_occ);
        //printf("%d %d %d %d\n", (int)dna_occ[0], (int)dna_occ[1], (int)dna_occ[2], (int)dna_occ[3]);
        //printf("%d %d %d %d\n", (int)str_occ[0], (int)str_occ[1], (int)str_occ[2], (int)str_occ[3]);
    }
    {
        const char *str = "TTACGATTATTAATTACGATTATTAATTACGATTATTAA";
        dna_string dna(str);
        {
            size_t b = 3, e = 39;
            auto dna_occ = dna.occurance(b, e);
            auto str_occ = count(str+b, str+e);
            BOOST_CHECK(dna_occ == str_occ);
        }
        {
            size_t b = 3, e = 4;
            auto dna_occ = dna.occurance(b, e);
            auto str_occ = count(str+b, str+e);
            BOOST_CHECK(dna_occ == str_occ);
        }
        {
            size_t b = 3, e = 31;
            auto dna_occ = dna.occurance(b, e);
            auto str_occ = count(str+b, str+e);
            BOOST_CHECK(dna_occ == str_occ);
        }
        {
            size_t b = 3, e = 32;
            auto dna_occ = dna.occurance(b, e);
            auto str_occ = count(str+b, str+e);
            BOOST_CHECK(dna_occ == str_occ);
        }
        {
            size_t b = 3, e = 33;
            auto dna_occ = dna.occurance(b, e);
            auto str_occ = count(str+b, str+e);
            BOOST_CHECK(dna_occ == str_occ);
        }
    }
    {
        const char *str = chr1;
        dna_string dna(str);
        auto dna_occ = dna.occurance(0, dna.size());
        auto str_occ = count(str, str+dna.size());
        BOOST_CHECK(dna_occ == str_occ);
    }
}
