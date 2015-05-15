"""
Test module for example Boost.Genetics Python wrapper
"""

import genetics  # must match name of module generated in bindings.cpp


def test_genetics():

    s = genetics.Fasta()
    print(s)

    #s.load('/home/andy/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa')
    s.load('1.fa')
    result = s.find_all('ACTGACTGACTG', 0)
    print(len(result))
    print(result)


if __name__ == "__main__":
    test_genetics()
