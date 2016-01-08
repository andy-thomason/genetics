import sys

bits = 8

for i in range(1<<bits):
    chars = [(i >> j*2) & 3 for j in range(bits/2)]
    occ = [0 for j in range(4)]
    for c in chars:
        occ[c] += 1
    value = occ[0] + (occ[1] << 6) + (occ[2] << 12)
    sys.stdout.write("0x%04x," % value)
    if i % 32 == 31:
        sys.stdout.write("\n")

for i in range(256):
    res = 8
    for j in range(8):
        if i & 128>>j:
           res = j
           break
    sys.stdout.write("%d," % res)
    if i % 32 == 31:
        sys.stdout.write("\n")
    
