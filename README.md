# Boost.Genetics
A boost library for genetic searches

The Boost genetics library provides data structures and algorithms for searching very large databases of genetic sequences.

The library is primarily designed to provide a reliable, maintainable service to find 100% of available matches. Readability is a primary focus and unit testing a priority through the Travis system.

Genetic data usually consists of the letters 'A', 'C', 'G' and 'T' with other optional characters denoting
special cases.

**Reference data** is usually in FASTA format which is a series of chromomosome names and bases.

**Read data** from the Illumina machines is in FASTQ format with large numbers of shorter sequences of around
30-200 bases. Currently this is around 100, but the trend is for longer sequences.

**Aligned data** the result of an aligner such as Tophat or Star is usually in the SAM or BAM (compressed SAM) format and is similar to the FASTA format but has additional metadata.

We provide both packed (2 bits/base) and unpacked data (eg. std::string) algorithms.

Most genetic searches revolve around the creation of Suffix Arrays which are a sort of all the substrings of a sequence to the end.

eg.

"hello" consists of the substrings "hello", "ello", "llo", "lo", "p" and ""

The suffix array is:

5 ""
1 "ello"
0 "hello"
3 "lo"
2 "llo"
4 "o"

See Wikipedia for a more comprehensive discussion.

For simplicity, we calculate a hybrid indexed suffix array which gives you all the substrings
starting with a fixed number of letters. For this we add an index of 2^(2N) entries for N
letters the suffixes sorted by address in the array instead of value. We hope to extend this
to optionally use a variable length and hash methods in the future as well as to support reduced
size indices for low memory machines.

This means that we can use two very expensive memory accesses to find a set of addresses
for a sequence and then a small number of very inexpensive accesses and a merge operation
to find exact matches of any length (N x M).

For inexact matches, we can filter the M matches we have to a much smaller set where we can
then do expensive tests on the reference itself.

For a very large number of potential errors, such as 6 errors in a 20 character string, it is
fastest to brute-force search the genome. Here we can make use of leading zero and population
count operations to rapidly screen large number of bases.

Many thanks to Paul Bristow, Riley Doyle, Andre and others for their contributions.
