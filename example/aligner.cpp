

#include <fstream>
#include <random>
#include <thread>
#include <future>
#include <atomic>
#include <cstdint>
#include <chrono>

#include <boost/genetics/fasta.hpp>
#include <boost/genetics/utils.hpp>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>



class aligner {
public:
    // implement the "aligner index" mode
    void index(int argc, char **argv) {
        using namespace boost::program_options;
        using namespace boost::genetics;
        using namespace boost::interprocess;

        // These are the options for the aligner.
        options_description desc("aligner index <file1.fa> <file2.fa> ... <-o index.bin>");
        desc.add_options()
            ("help", "produce help message")
            ("fasta-files", value<std::vector<std::string> >(), "input filename(s)")
            ("output-file,o", value<std::string>()->default_value("index.bin"), "output filename")
            ("num-index-chars,n", value<int>()->default_value(12), "number of chars in first stage index")
        ;

        positional_options_description pod;
        pod.add("fasta-files", -1);

        // Parse the command line.
        variables_map vm;
        command_line_parser clp(argc, argv);
        store(
            clp.options(desc).positional(pod).run(),
            vm
        );
        notify(vm);

        if (vm.count("help") || argc == 1) {
            std::cout << desc << "\n";
            return;
        }

        // Append the fasta files to the reference.
        fasta_file builder;
        auto fa_files = vm["fasta-files"].as<std::vector<std::string> >();
        for (auto &f : fa_files) {
            std::cerr << f << "\n";
            builder.append(f);
        }

        if (builder.get_string().size() == 0) {
            throw std::runtime_error("no fa files or empty reference");
        }

        // Build the index. This will take some time.
        builder.make_index(vm["num-index-chars"].as<int>());

        // Find the size of the reference file.
        writer sizer(nullptr, nullptr);
        builder.write_binary(sizer);
        size_t size = (size_t)sizer.get_ptr();

        // Get the file name to write to (defaults to index.bin).
        std::string filename = vm["output-file"].as<std::string>();

        // there may be a better way of creating a pre-sized writable file!
        {
            std::ofstream os(filename, std::ios_base::binary);
            if (os.bad()) {
                throw std::runtime_error("unable to write index file");
            }
            os.seekp(size-1);
            os.write("", 1);
        }

        // Use file mapping to write the index file directly to buffers.
        try {
            file_mapping fm(filename.c_str(), read_write);
            mapped_region region(fm, read_write, 0, size);
            char *p = (char*)region.get_address();
            char *end = p + region.get_size();
            writer w(p, end);
            builder.write_binary(w);
        } catch(interprocess_exception &) {
            throw std::runtime_error("unable to write to the index file");
        }
    }

    // Implement the "aligner align" function
    void align(int argc, char **argv) {
        using namespace boost::program_options;
        using namespace boost::genetics;
        using namespace boost::interprocess;

        // These are the options.
        options_description desc("aligner align <index.bin> <file1.fq> <file2.fq> <-o out.sam>");
        desc.add_options()
            ("help", "produce help message")
            ("index,i", value<std::string>()->default_value("index.bin"), "index file")
            ("fastq-files,q", value<std::vector<std::string> >(), "fastq files (2 max)")
            ("output-file,o", value<std::string>()->default_value("out.sam"), "output filename")
            ("num-threads,t", value<int>()->default_value(1), "number of threads to use.")
        ;

        positional_options_description pod;
        pod.add("fastq-files", 2);

        // Parse the command line.
        variables_map vm;
        command_line_parser clp(argc, argv);
        store(
            clp.options(desc).positional(pod).run(),
            vm
        );
        notify(vm);

        if (vm.count("help") || !vm.count("index") || argc == 1) {
            std::cout << desc << "\n";
            return;
        }

        // Create a conventional file output stream for the results.
        sam_file = std::ofstream(vm["output-file"].as<std::string>(), std::ios_base::binary);

        // Map in the index (usually index.bin)
        file_mapping fm(vm["index"].as<std::string>().c_str(), read_only);
        mapped_region region(fm, read_only);
        char *p = (char*)region.get_address();
        char *end = p + region.get_size();
        mapper m(p, end);
        mapped_fasta_file ref(m);

        // Get the fastq filenames (max of 2).
        // These contain sequences we want to align.
        auto fq_filenames = vm["fastq-files"].as<std::vector<std::string> >();

        if (fq_filenames.size() != 1 && fq_filenames.size() != 2) {
            throw std::runtime_error("expected one or two FASTQ files");
        }

        // Build a vector of input files.
        std::vector<std::ifstream> read_files;
        for (size_t i = 0; i != fq_filenames.size(); ++i) {
            auto &f = fq_filenames[i];
            std::cerr << f << "\n";
            read_files.emplace_back(f, std::ios_base::binary);
        }

        // Run on multiple threads.
        int num_threads = vm["num-threads"].as<int>();
        std::vector<std::thread> align_threads;
        std::mutex read_mutex;
        std::atomic<size_t> num_reads;
        std::atomic<size_t> num_merges;
        std::atomic<size_t> num_compares;
        std::atomic<size_t> num_matches;

        search_params params;
        params.max_distance = 5;
        params.max_gap = 0;
        params.max_results = 100;
        params.always_brute_force = false;
        params.never_brute_force = true;
        params.search_rev_comp = true;

        for (size_t i = 0; i != ref.get_num_chromosomes(); ++i) {
            const chromosome &c = ref.get_chromosome(i);
            sam_file << "@SQ	SN:" << c.name << "\tLN:" << c.num_leading_N + (c.end - c.start) + c.num_trailing_N << "\n";
        }
        sam_file << "@PG	ID:boost\tPN:boost\tVN:1.0\n";

        auto start_time = std::chrono::system_clock::now();
        for (int tid = 0; tid != num_threads; ++tid) {
            align_threads.emplace_back(
                [&](int tid) {
                    // Each thread has a buffer for reading a section of the fastq file.
                    // And an index to find names and keys.
                    aligner_thread at(read_files.size());

                    // Keep reading sections of the input file until we have
                    // processed them All.
                    for(;;) {
                        size_t max_reads = 0;
                        {
                            // Only one thread can read at a time.
                            std::unique_lock<std::mutex> lock(read_mutex);
                            at.read_components(read_files);

                            max_reads = at.namess[0].size()-1;
                            for (size_t i = 1; i < read_files.size(); ++i) {
                                max_reads = std::min(max_reads, at.namess[0].size()-1);
                            }

                            if (max_reads == 0) break;

                            // If we have read too much, seek back to allow the next thread
                            // to read the whole entry.
                            for (size_t i = 0; i != read_files.size(); ++i)  {
                                auto &buffer = at.buffers[i];
                                auto &stream = read_files[i];
                                auto &names = at.namess[i];
                                size_t bytes = (size_t)stream.gcount();
                                const char *end = buffer.data() + bytes;
                                if (names[max_reads] != end) {
                                    stream.seekg(names[max_reads] - end, std::ios_base::cur);
                                }
                            }
                        }

                        for (size_t read_idx = 0; read_idx != max_reads; ++read_idx) {
                            // read and align a pair of reads or a single read
                            for (size_t i = 0; i != read_files.size(); ++i)  {
                                at.align_read(i, read_idx, params, ref);
                                num_merges += at.stats.merges_done;
                                num_compares += at.stats.compares_done;
                                num_matches += at.resultss[i].size() != 0;
                            }

                            // Todo: add pair matching.

                            // Write the results for each input read.
                            for (size_t i = 0; i != read_files.size(); ++i)  {
                                at.write_sam(i, sam_file, ref);
                            }
                        }
                        num_reads += max_reads;
                    } // for(;;)
                },
                tid
            );
        }

        for (int i = 0; i != num_threads; ++i) {
            align_threads[i].join();
        }

        auto end_time = std::chrono::system_clock::now();
        std::cerr << std::chrono::nanoseconds(end_time - start_time).count() * 1e-9 << "s\n";
        std::cerr << num_reads << " pairs\n";
        std::cerr << (double)num_reads * 1000000000.0 / std::chrono::nanoseconds(end_time - start_time).count() << "pairs/s\n";
        if (num_reads) {
            std::cerr << (double)num_matches / num_reads << " matches/read\n";
            std::cerr << (double)num_merges / num_reads << " merges/read\n";
            std::cerr << (double)num_compares / num_reads << " compares/read\n";
        }
    }

    void generate_reference(int argc, char **argv) {
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

        const char *binary_filename = "index.bin";

        writer sizer(nullptr, nullptr);
        builder.write_binary(sizer);

        std::vector<char> buffer(sizer.get_size());
        writer wr(buffer.data(), buffer.data() + sizer.get_size());
        builder.write_binary(wr);

        std::ofstream out_file(binary_filename, std::ios_base::binary);
        out_file.write(buffer.data(), sizer.get_size());
    }

    void generate_fastq(int argc, char **argv) {
        using namespace boost::genetics;
        /*size_t key_size = 101;

        std::default_random_engine generator;
        std::uniform_int_distribution<size_t> loc_distribution(0, ref.size()-key_size);
        std::uniform_int_distribution<size_t> error_loc_distribution(0, key_size-1);
        std::uniform_int_distribution<size_t> error_prob_distribution(0, 99);
        std::uniform_int_distribution<int> error_value_distribution(0, 4);
        std::uniform_int_distribution<int> rev_comp_value_distribution(0, 1);
        std::uniform_int_distribution<int> phred_distribution('#', 'J');
        std::vector<fasta_result> result;

        const char *fastq_filename = "fastq.fq";
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
        */
    }
private:
    static const size_t buffer_size = 0x100000;
    size_t num_multiple = 0;
    size_t num_unmatched = 0;
    size_t num_reads = 0;
    std::ofstream sam_file;

    struct aligner_thread {
        std::vector<std::array<char, buffer_size> > buffers;
        std::vector<std::vector<const char*> > keyss;
        std::vector<std::vector<const char*> > namess;

        std::vector<std::string> name_strs;
        std::vector<std::string> key_strs;
        std::vector<std::string> phred_strs;
        std::vector< std::vector<boost::genetics::fasta_result> > resultss;

        std::string out_buf;

        boost::genetics::search_stats stats;

        aligner_thread(size_t num_files) {
            buffers.resize(num_files);
            keyss.resize(num_files);
            namess.resize(num_files);
            name_strs.resize(num_files);
            key_strs.resize(num_files);
            phred_strs.resize(num_files);
            resultss.resize(num_files);
        }

        void read_components(std::vector<std::ifstream> &read_files) {
            for (size_t i = 0; i != read_files.size(); ++i)  {
                auto &buffer = buffers[i];
                auto &stream = read_files[i];
                auto &keys = keyss[i];
                auto &names = namess[i];
                stream.read(buffer.data(), buffer_size);
                size_t bytes = (size_t)stream.gcount();
                bool is_last_block = stream.eof();
                //std::cerr << bytes << " " << is_last_block << " " << stream.tellg() << "\n";
                const char *p = buffer.data();
                const char *end = buffer.data() + bytes;
                keys.resize(0);
                names.resize(0);
                
                //printf("here\n");
                while (p != end) {
                    if (*p != '@') {
                        throw std::runtime_error("Not a FASTQ file: no @ found for name");
                    }
                    names.push_back(p);
                    while (p != end && *p != '\n') ++p;
                    p += p != end;
                    keys.push_back(p);
                    while (p != end && *p != '\n') ++p;
                    p += p != end;
                    while (p != end && *p != '\n') ++p;
                    p += p != end;
                    while (p != end && *p != '\n') ++p;
                    p += p != end;
                }
                if (is_last_block) names.push_back(end);
                size_t num_reads = names.size() - 1;
            }
        }

        void align_read(size_t file_idx, size_t read_idx, boost::genetics::search_params &params, boost::genetics::mapped_fasta_file &ref) {
            using namespace boost::genetics;

            auto &results = resultss[file_idx];
            auto &keys = keyss[file_idx];
            auto &names = namess[file_idx];
            std::string &name_str = name_strs[file_idx];
            std::string &key_str = key_strs[file_idx];
            std::string &phred_str = phred_strs[file_idx];

            // FASTQ reads have the form:
            // @name          name of this read
            // ACGGATTGAG     probable sequence from machine
            // +
            // JJC#DDGGGH     log-scale quality (Phred).

            const char *p = names[read_idx]+1, *q = p;
            while (*p != ' ' && *p != '\n') ++p;
            name_str.assign(q, p);

            p = keys[read_idx], q = p;
            while (*p != '\n') ++p;
            key_str.assign(q, p);

            ++p;
            while (*p != '\n') ++p;

            ++p;
            q = p;
            while (*p != '\n') ++p;
            phred_str.assign(q, p);

            stats.merges_done = 0;
            stats.compares_done = 0;
            ref.find_inexact(results, key_str, params, stats);
        }

        void write_sam(size_t file_idx, std::ofstream &sam_file, boost::genetics::mapped_fasta_file &ref) {
            using namespace boost::genetics;

            auto &results = resultss[file_idx];
            auto &keys = keyss[file_idx];
            auto &names = namess[file_idx];
            std::string &name_str = name_strs[file_idx];
            std::string &key_str = key_strs[file_idx];
            std::string &phred_str = phred_strs[file_idx];

            // No alignment, write a null record.
            if (results.size() == 0) {
                out_buf.resize(1024);
                auto dest = out_buf.begin();
                dest = make_str(dest, name_str.c_str());
                dest = make_str(dest, "\t4\t*\t0\t0\t*\t*\t0\t0\t");
                dest = make_str(dest, key_str.c_str());
                dest = make_str(dest, "\t");
                dest = make_str(dest, phred_str.c_str());
                dest = make_str(dest, "\n");
                sam_file.write(out_buf.c_str(), dest - out_buf.begin());
            } else {
                for (size_t res_idx = 0; res_idx != results.size(); ++res_idx) {
                    fasta_result &r = results[res_idx];
                    const chromosome &c = ref.find_chromosome(r.location);
                    size_t key_len = key_str.size();

                    // SAM flags
                    // https://ppotato.files.wordpress.com/2010/08/sam_output.pdf

                    int flags = 0;
                    flags |= r.reverse_complement ? 0x10 : 0x00;
                    flags |= res_idx != 0 ? 0x100 : 0x000;
                    int qual = results.size() > 2 ? 0 : 37;

                    out_buf.resize(0x10000);
                    auto dest = out_buf.begin();
                    dest = make_str(dest, name_str.c_str());
                    dest = make_str(dest, "\t");
                    dest = make_int(dest, flags);
                    dest = make_str(dest, "\t");
                    dest = make_str(dest, c.name);
                    dest = make_str(dest, "\t");
                    dest = make_int(dest, r.location - c.start + c.num_leading_N + 1);
                    dest = make_str(dest, "\t");
                    dest = make_int(dest, qual);
                    dest = make_str(dest, "\t");
                    dest = make_int(dest, key_str.size());
                    dest = make_str(dest, "M\t*\t0\t0\t");
                    if (r.reverse_complement) {
                        dest = make_rev_comp(dest, key_str);
                        dest = make_str(dest, "\t");
                        dest = make_rev_str(dest, phred_str.c_str());
                    } else {
                        dest = make_str(dest, key_str.c_str());
                        dest = make_str(dest, "\t");
                        dest = make_str(dest, phred_str.c_str());
                    }
                    dest = make_str(dest, "\tXT:A:U\tNM:i:");
                    dest = make_int(dest, r.distance);
                    dest = make_str(dest, "\tX0:i:1\tX1:i:0\tXM:i:");
                    dest = make_int(dest, r.distance);
                    dest = make_str(dest, "\tXO:i:0\tXG:i:0\tMD:Z:");
                    std::string ref_str = ref.get_string().substr(r.location, key_str.length(), r.reverse_complement);
                    dest = make_MD_field(dest, ref_str, key_str);
                    dest = make_str(dest, "\n");
                    sam_file.write(out_buf.c_str(), dest - out_buf.begin());
                }
            }
        }
    };
};


int main(int argc, char **argv) {
    aligner al;

    try {
        if (argc >= 2) {
            if (!strcmp(argv[1], "index")) {
                al.index(argc-1, argv+1);
                return 0;
            } else if (!strcmp(argv[1], "align")) {
                al.align(argc-1, argv+1);
                return 0;
            } else if (!strcmp(argv[1], "genreads")) {
                //al.generate(argc-1, argv+1);
                return 0;
            } else if (!strcmp(argv[1], "genref")) {
                //al.generate(argc-1, argv+1);
                return 0;
            } else {
                std::cerr << "unknown function " << argv[1] << "\n";
                return 1;
            }
        } else {
            std::cerr << "Usage:\n";
            std::cerr << "  aligner index <file1.fa> <file2.fa> ... <-o index.bin>       (Generate Index)\n";
            std::cerr << "  aligner align <index.bin> <file1.fq> <file2.fq> <-o out.sam> (Align against index) \n";
            //std::cerr << "  aligner genreads ...\n";
            //std::cerr << "  aligner genref ...\n";
            std::cerr << "  aligner <index|align>                                        (Get help for each function)\n";
            return 1;
        }
    } catch (boost::program_options::unknown_option u) {
        std::cerr << "error: " << u.what() << "\n";
        return 1;
    } catch(std::runtime_error e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

