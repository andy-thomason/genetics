import genetics

s = genetics.Fasta()
print("loading...")

s.load(['/projects/deskgen-rosalind/22.fa'], 10)

print("finding...")

result = s.find_inexact('ACTGACTGACTG', 0, 100)


