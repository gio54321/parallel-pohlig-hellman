import random
from Crypto.Util.number import getPrime


for i in range(8, 49, 4):
    for _ in range(5):
        p = getPrime(i)
        g = random.randint(2, p)
        x = random.randint(2, p)
        y = pow(g, x, p)
        print(f"{g} {y} {p}")
        assert pow(g, x, p) == y