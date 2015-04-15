
#include <fstream>
#include <sstream>

#include <boost/genetics/genetics.hpp>

#define BOOST_TEST_MODULE Genetics
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( packed )
{
    using namespace boost::genetics;

    {
        bases<8> b("ACGT");
        BOOST_CHECK( b.to_ullong() == 0x1b00000000000000ull );
        BOOST_CHECK( b.to_string() == "ACGTAAAA" );
    }

    {
        bases<4> b("TGCAACACACAC");
        BOOST_CHECK( b.to_ullong() == 0xe400000000000000ull );
        BOOST_CHECK( b.to_string() == "TGCA" );
    }

    {
        bases<100> b("TGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACACCCCCCCCCCCCC");
        //BOOST_MESSAGE( b );
        BOOST_CHECK( b.to_string() == "TGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACATGCAACACACAT");
    }

    {
        std::string str("tTAACcGGTTACGTACGTACGTACGTACGTACGTACGT");
        bases<64> b( str );
        BOOST_CHECK( b.to_ullong() == 0xf05af1b1b1b1b1b1ull );
    }

    {
        bases<16> b;
        std::stringstream ss("TGCAACGTttttXXXX\nAACCGGTT");
        ss >> b;
        //BOOST_MESSAGE( b );
        BOOST_CHECK( b.to_string() == "TGCAACGTTTTTAAAA" );
        ss >> b;
        BOOST_MESSAGE( b );
        BOOST_CHECK( b.to_string() == "AACCGGTTAAAAAAAA" );
    }

    {
        bases<32> b("ACGTACGTACGT");
        char buf[64];
        std::stringstream ss(buf);
        ss << b;
        BOOST_CHECK( b.to_string() == "ACGTACGTACGTAAAAAAAAAAAAAAAAAAAA" );
    }

    BOOST_CHECK( false );
}

