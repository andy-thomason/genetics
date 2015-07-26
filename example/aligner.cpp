

#include <fstream>
#include <random>
#include <thread>
#include <future>
#include <atomic>
#include <cstdint>
#include <chrono>

#include <boost/genetics/fasta.hpp>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>



class aligner {
public:
    void index(int argc, char **argv) {
        using namespace boost::program_options;
        using namespace boost::genetics;
        using namespace boost::interprocess;


        options_description desc("aligner index <file1.fa> <file2.fa> ... <-o index.bin>");
        desc.add_options()
            ("help", "produce help message")
            ("fasta-files", value<std::vector<std::string> >(), "input filename(s)")
            ("output-file,o", value<std::string>()->default_value("index.bin"), "output filename")
            ("num-index-chars,n", value<int>()->default_value(12), "number of chars in first stage index")
        ;

        positional_options_description pod;
        pod.add("fasta-files", -1);

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

        fasta_file builder;
        auto fa_files = vm["fasta-files"].as<std::vector<std::string> >();
        for (auto &f : fa_files) {
            std::cerr << f << "\n";
            builder.append(f);
        }

        if (builder.get_string().size() == 0) {
            throw std::runtime_error("no fa files or empty reference");
        }

        builder.make_index(vm["num-index-chars"].as<int>());


        writer sizer(nullptr, nullptr);
        builder.write_binary(sizer);
        size_t size = (size_t)sizer.get_ptr();


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
        file_mapping fm(filename.c_str(), read_write);
        mapped_region region(fm, read_write, 0, size);
        char *p = (char*)region.get_address();
        char *end = p + region.get_size();
        writer w(p, end);
        builder.write_binary(w);
    }

    void align(int argc, char **argv) {
        using namespace boost::program_options;
        using namespace boost::genetics;
        using namespace boost::interprocess;

        options_description desc("aligner align <index.bin> <file1.fq> <file2.fq> <-o out.sam>");
        desc.add_options()
            ("help", "produce help message")
            ("index,i", value<std::string>(), "index file")
            ("fastq-files,q", value<std::vector<std::string> >(), "fastq files (2 max)")
            ("output-file,o", value<std::string>()->default_value("out.sam"), "output filename")
        ;

        positional_options_description pod;
        pod.add("index", 1);
        pod.add("fastq-files", 2);

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

        sam_file = std::ofstream(vm["output-file"].as<std::string>(), std::ios_base::binary);

        file_mapping fm(vm["index"].as<std::string>().c_str(), read_only);
        mapped_region region(fm, read_only);
        char *p = (char*)region.get_address();
        char *end = p + region.get_size();
        mapper m(p, end);
        mapped_fasta_file ref(m);

        auto fq_filenames = vm["fastq-files"].as<std::vector<std::string> >();

        if (fq_filenames.size() != 1 && fq_filenames.size() != 2) {
            throw std::runtime_error("expected one or two FASTQ files");
        }

        std::vector<std::ifstream> read_files;
        for (size_t i = 0; i != fq_filenames.size(); ++i) {
            auto &f = fq_filenames[i];
            std::cerr << f << "\n";
            read_files.emplace_back(f, std::ios_base::binary);
        }

        int num_threads = 8;
        std::vector<std::thread> align_threads;
        std::mutex read_mutex;
        std::atomic<size_t> num_reads;
        std::atomic<size_t> num_merges;
        std::atomic<size_t> num_compares;
        std::atomic<size_t> num_matches;
        auto start_time = std::chrono::system_clock::now();
        for (int tid = 0; tid != num_threads; ++tid) {
            align_threads.emplace_back(
                [&](int tid) {
                    buffers.resize(read_files.size());
                    keyss.resize(read_files.size());
                    namess.resize(read_files.size());

                    //std::vector<size_t> ends(read_files.size());
                    for(;;) {
                        std::unique_lock<std::mutex> lock(read_mutex);
                        //std::cerr << "t" << tid << "\n";
                        read_components(read_files);

                        size_t max_reads = namess[0].size()-1;
                        for (size_t i = 1; i < read_files.size(); ++i) {
                            max_reads = std::min(max_reads, namess[0].size()-1);
                        }

                        if (max_reads == 0) break;

                        for (size_t i = 0; i != read_files.size(); ++i)  {
                            auto &buffer = buffers[i];
                            auto &stream = read_files[i];
                            auto &names = namess[i];
                            size_t bytes = (size_t)stream.gcount();
                            const char *end = buffer.data() + bytes;
                            //std::cerr << names[max_reads] - end << "\n";
                            if (names[max_reads] != end) {
                                stream.seekg(names[max_reads] - end, std::ios_base::cur);
                            }
                            //std::cerr << std::string(names[0], std::find(names[0], names[0]+0x1000, '\n')) << "\n";
                        }

                        //std::cerr << "max_reads = " << max_reads << "\n";
                        //lock.release();

                        std::array< std::vector<fasta_result>, 2 > resultss;

                        std::vector<std::string> key_strs(read_files.size());
                        std::string name_str;

                        search_params params;
                        params.max_distance = 2;
                        params.max_gap = 0;
                        params.max_results = 100;
                        params.always_brute_force = false;
                        params.never_brute_force = true;
                        params.search_rev_comp = true;
                        for (size_t read_idx = 0; read_idx != max_reads; ++read_idx) {
                            for (size_t i = 0; i != read_files.size(); ++i)  {
                                auto &results = resultss[i];
                                auto &keys = keyss[i];
                                auto &names = namess[i];
                                key_strs[i].assign(keys[read_idx], std::find(keys[read_idx], names[read_idx+1], '\n'));

                                params.stats.merges_done = 0;
                                params.stats.compares_done = 0;
                                ref.find_inexact(results, key_strs[i], params);

                                //if (params.stats.merges_done > 10000) {
                                    //std::cerr << read_idx << ": " << key_strs[i] << " " << params.stats.merges_done << "\n";
                                //}
                                num_merges += params.stats.merges_done;
                                num_compares += params.stats.compares_done;

                                num_matches += results.size() != 0;

                                /*name_str.assign(names[read_idx]+1, std::find(names[read_idx], names[read_idx+1], '\n'));

                                for (size_t res_idx = 0; res_idx != results.size(); ++res_idx) {
                                    fasta_result &r = results[res_idx];
                                    const chromosome &c = ref.find_chromosome(r.location);

                                    // SAM flags
                                    // https://ppotato.files.wordpress.com/2010/08/sam_output.pdf
                                    int flags = 0;
                                    flags |= r.reverse_complement ? 0x10 : 0x00;
                                    flags |= res_idx != 0 ? 0x100 : 0x000;
                                    sam_file
                                        << name_str << '\t'
                                        << flags << '\t'
                                        << c.name << '\t'
                                        << r.location - c.start + c.num_leading_N + 1 << '\t'
                                        << '\n'
                                    ;
                                }*/
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
        std::cerr << std::chrono::nanoseconds(end_time - start_time).count() << "ns\n";
        std::cerr << num_reads << " pairs\n";
        std::cerr << (double)num_reads * 1000000000.0 / std::chrono::nanoseconds(end_time - start_time).count() << "pairs/s\n";
        std::cerr << (double)num_matches / num_reads << " matches/read\n";
        std::cerr << (double)num_merges / num_reads << " merges/read\n";
        std::cerr << (double)num_compares / num_reads << " compares/read\n";
    }

    void generate_reference((int argc, char **argv) {
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

    void generate_fastq((int argc, char **argv) {
        using namespace boost::genetics;

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
    }
private:
    static const size_t buffer_size = 0x100000;
    std::vector<std::array<char, buffer_size> > buffers;
    std::vector<std::vector<const char*> > keyss;
    std::vector<std::vector<const char*> > namess;
    size_t num_multiple = 0;
    size_t num_unmatched = 0;
    size_t num_reads = 0;
    std::ofstream sam_file;

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
            while (p != end) {
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

