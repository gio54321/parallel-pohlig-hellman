import random
from Crypto.Util.number import getPrime, isPrime


"""
Generate a prime p such that p-1 is k(bit)-smooth
"""
def gen_k_smooth_prime(k, nbits, max_e = 3):
    while True:
        # 2 is the only even prime
        # it is necessary to generate a prime p such that p-1 is k-smooth
        # because if p-1 is not even, then p cannot possibly be a prime
        result = 2
        factors = [2]
        es = [1]
        
        while result.bit_length() < nbits - k*max_e:
            e = random.randint(1, max_e)
            p = getPrime(k)
            result *= p ** e
            factors.append(p)
            es.append(e)

        while result.bit_length() < nbits - k*2:
            p = getPrime(k)
            result *= p
            factors.append(p)
            es.append(1)

        i = 0
        while True:
            i += 1
            if i > 1000:
                # failed to find a suitable prime, retry
                break

            # make sure that the bit length of the largest prime factor
            # is exactly k
            if nbits - result.bit_length() < k:
                p1 = getPrime(nbits - result.bit_length())
                candidate = result * p1
                p2 = -1
            else:
                p1 = getPrime((nbits - result.bit_length())//2)
                p2 = getPrime((nbits - result.bit_length())//2)
                candidate = result * p1 * p2

            if candidate.bit_length() != nbits:
                continue

            
            if isPrime(candidate + 1):
                factors.append(p1)
                es.append(1)
                if p2 != -1:
                    factors.append(p2)
                    es.append(1)
                return candidate + 1, factors, es


test_bits = [40]
max_es = [1, 3]

for e in max_es:
    for nbits in test_bits:
        for i in range(10):
            print(f"Generating {nbits} bits smooth input {i}, max_e = {e}")
            p, factors, es =  gen_k_smooth_prime(nbits, 512, max_e=e)
            g = random.randint(2, p)
            x = random.randint(2, p)
            y = pow(g, x, p)
            assert pow(g, x, p) == y

            f = open(f"input/pohlig_hellman_{nbits}_{e}_{i}.txt", "w")
            f.write(f"{g} {y} {p}\n")
            f.write(str(len(factors)) + ' ' + ' '.join([f"{x} {y}" for x, y in zip(factors, es)]))
            f.close()