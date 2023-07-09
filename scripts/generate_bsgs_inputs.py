import random
from Crypto.Util.number import getPrime

test_nbits = [40, 48, 52]

for nbits in test_nbits:
    for i in range(5):
        p = getPrime(nbits)
        g = random.randint(2, p)
        x = random.randint(2, p)
        y = pow(g, x, p)
        assert pow(g, x, p) == y

        f = open(f'input/bsgs_{nbits}_{i}.txt', 'w')
        f.write(f"{g} {y} {p}")
        f.close()