import random
from Crypto.Util.number import getPrime, isPrime


"""
Generate a prime p such that p-1 is k(bit)-smooth
"""
def gen_k_smooth_prime(k, nbits):
    while True:
        # 2 is the only even prime
        # it is necessary to generate a prime p such that p-1 is k-smooth
        # because if p-1 is not even, then p cannot possibly be a prime
        result = 2
        factors = [2]
        es = [1]
        
        while result.bit_length() < nbits - k*3:
            e = random.randint(1, 3)
            p = getPrime(k)
            result *= p ** e
            factors.append(p)
            es.append(e)

        while result.bit_length() < nbits - k*2:
            p = getPrime(k)
            result *= p
            factors.append(p)
            es.append(1)

        #print(result.bit_length())
        #print(nbits - result.bit_length())
        
        if nbits - result.bit_length() < k:
            continue
        i = 0
        while True:
            i += 1
            if i > 1000:
                # print("failed")
                break
            # make sure that the bit length of the largest prime factor
            # is exactly k
            p1 = getPrime((nbits - result.bit_length())//2)
            p2 = getPrime((nbits - result.bit_length())//2)
            candidate = result * p1 * p2

            if candidate.bit_length() != nbits:
                continue

            
            if isPrime(candidate + 1):
                factors.append(p1)
                factors.append(p2)
                es.append(1)
                es.append(1)
                return candidate + 1, factors, es

    
for i in range(16, 49, 4):
    p, factors, es =  gen_k_smooth_prime(i, 512)
    g = random.randint(2, p)
    x = random.randint(2, p)
    y = pow(g, x, p)
    print(f"# largest prime factor has {i} bits")
    print(f"{g} {y} {p}")
    print(str(len(factors)) + ' ' + ' '.join([f"{x} {y}" for x, y in zip(factors, es)]))
    print()
    assert pow(g, x, p) == y