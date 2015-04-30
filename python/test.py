import genetics

s = genetics.AugmentedString()
print(s)

s.load_fasta('~/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa')
result = s.find_all('ACGT', 0)
print(len(result))


