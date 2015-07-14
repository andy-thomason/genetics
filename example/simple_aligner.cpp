

#include <fstream>
#include <random>

#include <boost/genetics/fasta.hpp>
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

const char *binary_filename = "index.bin";
const char *fastq_filename = "random.fq";
const char *sam_filename = "random.sam";
const size_t key_size = 101;
const size_t max_distance = 2;
const size_t max_results = 100;
const size_t max_gap = 0;

void generate_index() {
    using namespace boost::genetics;

    fasta_file builder;
    builder.append_random("chr1", 1000000, 0x9bac7615+0);
    builder.append_random("chr2",  900000, 0x9bac7615+1);
    builder.append_random("chr3",  800000, 0x9bac7615+2);
    builder.append_random("chr4",  700000, 0x9bac7615+3);
    builder.append_random("chr5",  600000, 0x9bac7615+4);
    builder.append_random("chr6",  500000, 0x9bac7615+5);
    builder.append_random("chr7",  400000, 0x9bac7615+6);
    builder.append_random("chr8",  300000, 0x9bac7615+7);
    builder.append_random("chr9",  200000, 0x9bac7615+0);
    builder.make_index(10);

    writer sizer(nullptr, nullptr);
    builder.write_binary(sizer);

    std::vector<char> buffer(sizer.size());
    writer wr(buffer.data(), buffer.data() + sizer.size());
    builder.write_binary(wr);

    std::ofstream out_file(binary_filename, std::ios_base::binary);
    out_file.write(buffer.data(), sizer.size());
}

void generate_fastq(boost::genetics::mapped_fasta_file &ref) {
    using namespace boost::genetics;

    std::default_random_engine generator;
    std::uniform_int_distribution<size_t> loc_distribution(0, ref.size()-key_size);
    std::uniform_int_distribution<size_t> error_loc_distribution(0, key_size-1);
    std::uniform_int_distribution<size_t> error_prob_distribution(0, 99);
    std::uniform_int_distribution<int> error_value_distribution(0, 4);
    std::uniform_int_distribution<int> rev_comp_value_distribution(0, 1);
    std::uniform_int_distribution<int> phred_distribution('#', 'J');
    std::vector<fasta_result> result;

    std::ofstream fastq_file(fastq_filename);
    if (!fastq_file.good()) {
        std::cerr << "unable to create " << fastq_filename << "\n";
        return;
    }

    std::string phred(key_size, ' ');
    for (int i = 0; i != 1000000; ) {
        size_t locus = loc_distribution(generator);
        size_t error_prob = error_prob_distribution(generator);
        int rev_comp = rev_comp_value_distribution(generator);
        std::string key = ref.get_string().substr(locus, key_size, rev_comp != 0);
        if (std::count(key.begin(), key.end(), 'N') < 3) {
            // Inject up to two random errors or an 'N'
            if (error_prob < 10) {
                size_t error_loc = error_loc_distribution(generator);
                int error_value = error_value_distribution(generator);
                key[error_loc] = "ACGTN"[error_value];
                if (error_prob < 1) {
                    size_t error_loc = error_loc_distribution(generator);
                    int error_value = error_value_distribution(generator);
                    key[error_loc] = "ACGTN"[error_value];
                }
            }
            for (int j = 0; j != key_size; ++j) {
                phred[j] = phred_distribution(generator);
            }
            fastq_file << "@r" << i << "." << locus << "." << rev_comp << "\n" << key << "\n+\n" << phred << "\n";
            ++i;
        }
    }
}

//! Very basic aligner.
void align(std::ifstream &file, boost::genetics::mapped_fasta_file &ref) {
    using namespace boost::genetics;

    std::ofstream sam_file(sam_filename, std::ios_base::binary);

    char name[80];
    char key[key_size*2+2];
    char plus[8];
    char phred[key_size*2+2];
    std::vector<fasta_result> result;
    std::string rev;
    size_t num_multiple = 0;
    size_t num_unmatched = 0;
    size_t num_reads = 0;
    while(!file.eof()) {
        file.getline(name, sizeof(name));
        if (file.eof()) break;
        file.getline(key, sizeof(key));
        file.getline(plus, sizeof(plus));
        file.getline(phred, sizeof(phred));
        ref.find_inexact(result, key, max_distance, max_results, max_gap, false, true);

        num_unmatched += result.size() == 0 ? 1 : 0;
        num_multiple += result.size() >= 2 ? 1 : 0;
        num_reads++;

        if (result.size() == 0) {
            std::cerr << name << ", " << key << ", " << file.eof() << " unmatched\n";
        }

        for (size_t i = 0; i != result.size(); ++i) {
            fasta_result &r = result[i];
            const chromosome &c = ref.find_chromosome(r.location);

            // SAM flags
            // https://ppotato.files.wordpress.com/2010/08/sam_output.pdf
            int flags = 0;
            flags |= r.reverse_complement ? 0x10 : 0x00;
            flags |= i != 0 ? 0x100 : 0x000;
            sam_file
                << name + 1 << '\t'
                << flags << '\t'
                << c.name << '\t'
                << r.location - c.start + 1 << '\t'
                << '\n'
            ;
        }
    }
    std::cerr << num_unmatched << " unmatched\n";
    std::cerr << num_multiple << " multiple\n";
    std::cerr << num_reads << " reads\n";
}


//! This example shows how to write a simple single-threaded aligner.
//! We also generate some 
int main() {
    using namespace boost::genetics;
    using namespace boost::interprocess;

    // If the index does not exist, use traditional storage and files
    // to create it.
    {
        std::ifstream test_file(binary_filename);
        if (!test_file.good()) {
            generate_index();
        }
    }

    // after we have built the index, map it in instantly.
    file_mapping fm(binary_filename, read_only);
    mapped_region region(fm, read_only);
    char *p = (char*)region.get_address();
    char *end = p + region.get_size();
    mapper m(p, end);
    mapped_fasta_file ref(m);

    // If the FASTQ file does not exist, generate it.
    {
        std::ifstream test_file(fastq_filename);
        if (!test_file.good()) {
            generate_fastq(ref);
        }
    }

    {
        std::ifstream test_file(fastq_filename);
        if (test_file.good()) {
            align(test_file, ref);
        }
    }
}

