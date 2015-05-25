"""
Test module for example Boost.Genetics Python wrapper
"""

import genetics  # must match name of module generated in bindings.cpp


if False:
  print("loading...")
  s = genetics.Fasta(['/projects/deskgen-rosalind/22.fa'], 10)
  print("writing...")
  s.write_binary_file("1.bin")
else:
  print("loading...")
  s = genetics.Fasta("1.bin")

  print("finding...")

  result = s.find_inexact('ACTGACTGACTG', 0, 100)
  print(len(result))
  print(result)
