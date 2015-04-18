
#include <fstream>
#include <sstream>

#include <boost/genetics/genetics.hpp>

#define BOOST_TEST_MODULE Genetics
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test_bases )
{
    using namespace boost::genetics;

    {
        bases b("ACGT");
        BOOST_MESSAGE( b );
        BOOST_CHECK( b == "ACGT" );
        b.append("TTGA");
        BOOST_MESSAGE( b );
        BOOST_CHECK( b == "ACGTTTGA" );
        b.append("ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG");
        BOOST_MESSAGE( b );
        BOOST_CHECK( b == "ACGTTTGA" "ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG" );
        BOOST_CHECK( b[-1] == 'N' );
        BOOST_CHECK( b[0] == 'A' );
        BOOST_CHECK( b[31] == 'G' );
        BOOST_CHECK( b[32] == 'A' );
        BOOST_CHECK( b[33] == 'C' );
        BOOST_CHECK( b[999] == 'N' );
    }

    {
        bases b("ACGT");
        b.resize(3);
        b.resize(4);
        BOOST_MESSAGE( b );
        BOOST_CHECK( b == "ACGA" );
    }

    {
        bases b("TGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACAT");
        BOOST_MESSAGE( b );
        BOOST_CHECK( b == "TGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACAT");
    }

    {
        std::string str("tTAACcGGTTACgtacgtacGTACGTACGTACGTACGT");
        bases b( str );
        BOOST_CHECK( b == "TTAACCGGTTACGTACGTACGTACGTACGTACGTACGT" );
    }

    {
        bases b;
        std::stringstream ss("TGCAACGTttttXXXX\nAACCGGTT");
        ss >> b;
        //BOOST_MESSAGE( b );
        BOOST_CHECK( b == "TGCAACGTTTTTAAAA" );
        ss >> b;
        BOOST_MESSAGE( b );
        BOOST_CHECK( b == "AACCGGTT" );
    }

    {
        bases b("ACGTACGTACGT");
        char buf[64] = {0};
        std::stringstream ss(buf);
        ss << b;
        BOOST_MESSAGE( std::string(ss.str()) );
        BOOST_CHECK( std::string(ss.str()) == "ACGTACGTACGT" );
    }

}

/*BOOST_AUTO_TEST_CASE( rev_comp )
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
}*/

BOOST_AUTO_TEST_CASE( substr )
{
    using namespace boost::genetics;

    {
        bases b("ACGTTTGA" "ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG");
        bases c;
        c = b.substr(8, 8);
        BOOST_MESSAGE( c );
        BOOST_CHECK( c == "ACCGACCG" );
        c = b.substr(8, 40);
        BOOST_MESSAGE( c );
        BOOST_CHECK( c == "ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG");
    }

    {
        std::string b("ACGTTTGA" "ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG");
        std::string c;
        c = b.substr(8, 8);
        BOOST_MESSAGE( c );
        BOOST_CHECK( c == "ACCGACCG" );
        c = b.substr(8, 40);
        BOOST_MESSAGE( c );
        BOOST_CHECK( c == "ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG");
    }

}

