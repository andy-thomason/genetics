clang++ -fpic -std=c++11 -mpopcnt -mlzcnt -O3 -I ../include packed_test.cpp -l boost_test_exec_monitor -o packed_test
clang++ -fpic -std=c++11 -mpopcnt -mlzcnt -O3 -I ../include fasta_test.cpp -l boost_test_exec_monitor -o fasta_test
