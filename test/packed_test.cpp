
#include <fstream>
#include <sstream>

#define BOOST_TEST_MODULE Genetics
#include <boost/test/unit_test.hpp>

#include <boost/genetics/genetics.hpp>

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
        //BOOST_MESSAGE( b );
        BOOST_CHECK( b == "TGCAACGTTTTTAAAA" );
        ss >> b;
        BOOST_MESSAGE( b );
        BOOST_CHECK( b == "AACCGGTT" );
    }

}

template <class Type>
void generic_container_tests() {
    {
        Type b("ACGT");
        BOOST_MESSAGE( b );
        BOOST_CHECK( std::string(b) == "ACGT" );
        b.append("TTGA");
        BOOST_MESSAGE( b );
        BOOST_CHECK( std::string(b) == "ACGTTTGA" );
        b.append("ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG");
        BOOST_MESSAGE( b );
        BOOST_CHECK( std::string(b) == "ACGTTTGA" "ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG" );
        //BOOST_CHECK( b[-1] == 'N' );
        BOOST_CHECK( b[0] == 'A' );
        BOOST_CHECK( b[31] == 'G' );
        BOOST_CHECK( b[32] == 'A' );
        BOOST_CHECK( b[33] == 'C' );
        //BOOST_CHECK( b[999] == 'N' );
    }
    {
        Type b("ACGT");
        b.resize(3);
        b.resize(4);
        BOOST_MESSAGE( b );
        BOOST_CHECK( b == "ACGA" );
    }

    {
        Type b("TGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACAT");
        BOOST_MESSAGE( b );
        BOOST_CHECK( b == "TGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACAT");
    }

    {
        Type b("ACGTACGTACGT");
        char buf[64] = {0};
        std::stringstream ss(buf);
        ss << b;
        BOOST_MESSAGE( std::string(ss.str()) );
        BOOST_CHECK( std::string(ss.str()) == "ACGTACGTACGT" );
    }
}

BOOST_AUTO_TEST_CASE( acgt_tests )
{
    using namespace boost::genetics;

    //acgt_container_tests<bases>();
    //acgt_container_tests<bases_and_letters>();
    //acgt_container_tests<std::string>();

    //generic_container_tests<bases_and_letters>();
    //generic_container_tests<std::string>();
}

BOOST_AUTO_TEST_CASE( generic_tests )
{
    using namespace boost::genetics;

    //acgt_container_tests<bases>();
    //acgt_container_tests<bases_and_letters>();
    //acgt_container_tests<std::string>();

    //generic_container_tests<bases_and_letters>();
    //generic_container_tests<std::string>();
}

BOOST_AUTO_TEST_CASE( bandl_tests )
{
    using namespace boost::genetics;

    {
        static const char test_string[] = "ACGTACGTXXXXACGT";
        bases_and_letters b(test_string);
        BOOST_CHECK(std::string(b) == test_string);
    }

    {
        static const char test_string[] = "NNNNACTGACTGACTG";
        bases_and_letters b(test_string);
        BOOST_CHECK(std::string(b) == test_string);
    }

    {
        std::string test_string(0xffff, 'N');
        test_string.append(256, 'A');
        bases_and_letters b(test_string);
        BOOST_CHECK(std::string(b) == test_string);
    }

    {
        std::string test_string(0x10000, 'N');
        test_string.append(256, 'A');
        bases_and_letters b(test_string);
        BOOST_CHECK(std::string(b) == test_string);
    }

}

BOOST_AUTO_TEST_CASE( rev_comp_test )
{
    using namespace boost::genetics;

    {
        bases b("TGCAACACACA");
        BOOST_MESSAGE( rev_comp(b) );
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
        bases b(test_string);
        bases c;
        c = b.substr(8, 8);
        BOOST_MESSAGE( c );
        BOOST_CHECK( c == "ACCGACCG" );
        c = b.substr(8, 40);
        BOOST_MESSAGE( c );
        bases d(test_string + 8);
        BOOST_MESSAGE( bases(test_string + 8) );
        BOOST_CHECK( c == test_string + 8);
    }

    {
        std::string b(test_string);
        std::string c;
        c = b.substr(8, 8);
        BOOST_MESSAGE( c );
        BOOST_CHECK( c == "ACCGACCG" );
        c = b.substr(8, 40);
        BOOST_MESSAGE( c );
        BOOST_CHECK( c == test_string + 8);
    }

}

