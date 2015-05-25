"""
Test module for example Boost.Genetics Python wrapper
"""

import genetics  # must match name of module generated in bindings.cpp

print("loading...")

s = genetics.Fasta(['/projects/deskgen-rosalind/22.fa'], 10)

s.write_binary("")

print("finding...")

result = s.find_inexact('ACTGACTGACTG', 0, 100)
print(len(result))
print(result)
