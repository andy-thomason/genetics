
#include <fstream>
#include <sstream>
#include <utility>

#define BOOST_TEST_MODULE Genetics
#include <boost/test/unit_test.hpp>

#include <boost/genetics/genetics.hpp>

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
        //printf("%016llx!!\n", (long long)w0);
        //printf("%016llx\n", (long long)w1);
        //printf("%016llx\n", (long long)w2);
        //printf("%016llx\n", (long long)w3);
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
        //BOOST_MESSAGE( b );
        BOOST_CHECK( b == "TGCAACGTTTTTAAAA" );
        ss >> b;
        BOOST_MESSAGE( b );
        BOOST_CHECK( b == "AACCGGTT" );
    }

}

template <class Type>
void generic_string_tests() {
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

BOOST_AUTO_TEST_CASE( acgt_tests )
{
    using namespace boost::genetics;

    find_test<std::string>();
    find_test<dna_string>();
}

BOOST_AUTO_TEST_CASE( generic_tests )
{
    using namespace boost::genetics;

    //acgt_container_tests<dna_string>();
    //acgt_container_tests<augmented_string>();
    //acgt_container_tests<std::string>();

    //generic_string_tests<augmented_string>();
    //generic_string_tests<std::string>();

    {
        dna_string a("ACGTACGTACGTACGTACGTACGTACGTACCA" "TTGTACGTACGTACGTACGTACGTACGTACGT");
        dna_string key1("CAT");
        BOOST_CHECK( a.find(key1, 0, 31) == dna_string::npos);
        BOOST_CHECK( a.find(key1, 0, 32) == dna_string::npos);
        BOOST_CHECK( a.find(key1, 0, 33) == (size_t)30);
        BOOST_CHECK( a.find(key1, 30, 33) == (size_t)30);
        BOOST_CHECK( a.find(key1, 31, 33) == dna_string::npos);
        dna_string key2("TTG");
        BOOST_CHECK( a.find(key2, 0, 33) == dna_string::npos);
        BOOST_CHECK( a.find(key2, 0, 34) == dna_string::npos);
        BOOST_CHECK( a.find(key2, 0, 35) == (size_t)32);
        dna_string key3("ACA");
        BOOST_CHECK( a.find(key3, 0, 64, 1) == 0); // matches ACG
        BOOST_CHECK( a.find(key3, 1, 64, 1) == 4); // matches ACG
        dna_string key4("ACCGTTGT");
        BOOST_CHECK( a.find(key4, 0, 64, 1) == 28); // matches "ACCATTGT"
    }
}

BOOST_AUTO_TEST_CASE( bandl_tests )
{
    using namespace boost::genetics;

    {
        static const char test_string[] = "ACGTACGTXXXXACGT";
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
        BOOST_MESSAGE( rev_comp(b) );
        BOOST_CHECK( rev_comp(b) == "ACGTACGTACGTACGTACGTACGTACGTACGT" );
    }

    {
        dna_string b("TGCAACACACA");
        BOOST_MESSAGE( rev_comp(b) );
        BOOST_CHECK( rev_comp(b) == "TGTGTGTTGCA" );
    }

    {
        std::string b("TGNAACACACA");
        BOOST_CHECK( rev_comp(b) == "TGTGTGTTNCA" );
    }
    {
        augmented_string b("TGNAACACACA");
        BOOST_MESSAGE( rev_comp(b) );
        BOOST_CHECK( rev_comp(b) == "TGTGTGTTNCA" );
    }
}



BOOST_AUTO_TEST_CASE( substr )
{
    using namespace boost::genetics;

    const char test_string[] = "ACGTTTGA" "ACCGACCG" "ACCGACCG" "TGCGACCG" "ACCGACCG";

    {
        dna_string b(test_string);
        dna_string c;
        c = b.substr(8, 8);
        BOOST_MESSAGE( c );
        BOOST_CHECK( c == "ACCGACCG" );
        c = b.substr(8, 40);
        BOOST_MESSAGE( c );
        dna_string d(test_string + 8);
        BOOST_MESSAGE( dna_string(test_string + 8) );
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

/*BOOST_AUTO_TEST_CASE( two_stage_index_test )
{
    using namespace boost::genetics;
    
    augmented_string as(chr1);
    two_stage_index<augmented_string> tsi(as, 4);
    //BOOST_MESSAGE(tsi);

    augmented_string key1("TCGAGACCATCCTGGCTAACACGGGGAAACCCCGTCTCCACTAAAAATACAAAAAGTTAG");
    if (0) {
      two_stage_index<augmented_string>::iterator i = tsi.find(key1, 0, 0);
      BOOST_CHECK(i == 120);
      ++i;
      BOOST_CHECK(i == 1200);
      ++i;
      BOOST_CHECK(i == augmented_string::npos);
    }
    if (0) {
      two_stage_index<augmented_string>::iterator i = tsi.find(key1, 121, 1);
      BOOST_CHECK(i == 1200);
      ++i;
      BOOST_CHECK(i == 1260);
      ++i;
      BOOST_CHECK(i == augmented_string::npos);
    }
    {
      two_stage_index<augmented_string>::iterator i = tsi.find(key1, 1201, 2);
      BOOST_CHECK(i == 1260);
      ++i;
      BOOST_CHECK(i == 1320);
      ++i;
      BOOST_CHECK(i == augmented_string::npos);
    }
}*/
