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

## Classes

Class                  | Description
--------------------   | -------------------------------------------------------------------------------
```dna_string```       | represents genetic data as two bits per base.
```augmented_string``` | adds auxililiary characters such as long runs on 'N' to the dna_string
```two_stage_index```  | provides very large (Gigabyte sized) index of a dna_string
```fm_index```         | is an implementation of an FM-Index used by the BWA and Bowtie aligners
```fasta```            | is an implementation of a genetic reference such as the ENSEMBL human genome

## Examples

```
#include <boost/genetics/fasta.hpp>
using namespace boost::genetics;

  // load an EBSEMBL reference
  // from ftp://ftp.ensembl.org/pub/release-83/fasta/homo_sapiens/dna/
  fasta my_reference("Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz");
  my_reference.make_index(12);

  // search for a single string.
  search_params params;
  search_stats stats;
  std::vector<fasta_result> result;
  my_reference.find_inexact(result, "GATACAGATACAGATACA", params, stats);

  // now results contains locations with strings close to "GATACAGATACAGATACA"

```


## Performance

Several measures improve performance significantly over traditional implementations.

**Memory mapped I/O** allows us to load an index instantly and share them between processes.
**Use of special instructions** Popcnt and lzcnt allow us to search faster
**Careful use of memory** Avoiding indexing when possible to improve cache latency

