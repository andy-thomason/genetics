
#include <fstream>
#include <sstream>
#include <utility>

#define BOOST_TEST_MODULE Genetics
#include <boost/test/unit_test.hpp>

#include <boost/genetics/genetics.hpp>

static const char chr1[] =
  "TGTGATTAATGCCTGAGACTGTGTGAAGTAAGAGATGGATCAGAGGCCGGGCGCGGGGGC"
  "TCGCGCCTGTCATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGA"
  "TCGAGACCATCCTGGCTAACACGGGGAAACCCCGTCTCCACTAAAAATACAAAAAGTTAG"
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
  "TCGAGACCATCCTGGCTAACACGGGGAAACCCCGTCTCCACTAAAAATACAAAAAGTTAG"
  "TCGAGACCATCCTGGCTAACACGCGGAAACCCCGTCTCCACTAAAAATACAAAAAGTTAG"
  "TCGAGACCATCCTGGCTAACACGCGGAAACCCCGTCTCCACTAATAATACAAAAAGTTAG"
;

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

    //acgt_container_tests<dna_string>();
    //acgt_container_tests<augmented_string>();
    //acgt_container_tests<std::string>();

    //generic_string_tests<augmented_string>();
    //generic_string_tests<std::string>();
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

BOOST_AUTO_TEST_CASE( two_stage_index_test )
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
    {
      two_stage_index<augmented_string>::iterator i = tsi.find(key1, 121, 1);
      BOOST_CHECK(i == 1200);
      ++i;
      BOOST_CHECK(i == 1260);
      ++i;
      BOOST_CHECK(i == 1320);
      ++i;
      BOOST_CHECK(i == augmented_string::npos);
    }
}


#if 0
template<class Type>
class alloc : public std::char_traits<char> {
public:
	typedef alloc<Type> other;

	typedef typename Type value_type;

	typedef value_type *pointer;
	typedef const value_type *const_pointer;
	typedef void *void_pointer;
	typedef const void *const_void_pointer;

	typedef value_type& reference;
	typedef const value_type& const_reference;

	typedef size_t size_type;
	typedef ptrdiff_t difference_type;

	typedef std::false_type propagate_on_container_copy_assignment;
	typedef std::false_type propagate_on_container_move_assignment;
	typedef std::false_type propagate_on_container_swap;

  /*static size_t copy(Type *d, const Type *a, size_t b) {
      while (b) *d++ = *a++;
      return 0;
  }

  static void assign(Type &d, const Type &a) {
      d = a;
  }*/
};
#endif

char ptr[10];

template<class Type>
class zalloc {
public:
  //zalloc();
	typedef typename Type value_type;

	typedef value_type *pointer;
	typedef const value_type *const_pointer;
	typedef void *void_pointer;
	typedef const void *const_void_pointer;

	typedef value_type& reference;
	typedef const value_type& const_reference;

	typedef size_t size_type;
	typedef ptrdiff_t difference_type;

  template <class U> struct rebind { typedef zalloc<U> other; };

  pointer allocate(size_t n) {
      return (pointer)ptr;
  }

  void deallocate(pointer p, size_t n) {
  }

  template <class... Args>
  void construct(pointer p, Args... a) {
      new ((void*)p)value_type(std::forward<Args>(a)...);
  }

  void destroy(pointer p) {
      p->~value_type();
  }
};

BOOST_AUTO_TEST_CASE( allocators )
{
    char x[10];
    zalloc<char> a;
    std::basic_string<char, std::char_traits<char>, zalloc<char>> str(a);
    str.resize(10);
}
