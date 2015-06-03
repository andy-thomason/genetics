"""
Test module for example Boost.Genetics Python wrapper
"""

import genetics  # must match name of module generated in bindings.cpp
import sys


if sys.argv[1] == "index":
  print("loading...")
  s = genetics.Reference(['Homo_sapiens.GRCh38.dna.primary_assembly.fa'], 14)
  print("writing...")
  s.write_binary_file("GRCh38.bin")
else:
  print("loading...")
  s = genetics.Reference(sys.argv[1])

  #s.write_ascii_file("1.fa")

  print("finding...")

  is_brute_force = int(sys.argv[3]) != 0
  print(is_brute_force)
  result = s.find_inexact('ACTGACTGACTGACTGACTGACTG', int(sys.argv[2]), 0, is_brute_force, 1000)
  print(len(result))
  print(result)
