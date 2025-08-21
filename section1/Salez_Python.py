# This script produces R_k.

from math import gcd
from sympy import isprime, divisors

def chrem(r, m):
    r0, r1, r2, m0, m1, m2, k = 0, 0, 0, 0, 0, 0, 0
    if m[0] > m[1]:
        m0 = m[1]
        r0 = r[1]
        m1 = m[0]
        r1 = r[0]
    else:
        m0 = m[0]
        r0 = r[0]
        m1 = m[1]
        r1 = r[1]
    m2 = pow(m1, -1, m0)
    if m2 == 0:
        return 0
    r2 = m0 - r1 % m0
    k = (r0 + r2) * m2
    k %= m0
    return (r1 + k * m1)


def modsqrt(s, n):
    Q, q = [], 1
    s %= n
    while q < (n + 1) // 2:
        if (q * q) % n == s: 
            Q.extend([q, n - q])
        q += 1
    return Q


def filter1a(m, S):
    m1 = None
    if m % 4 == 3:
        m1 = m
    elif m % 3 != 0:
        m1 = 3 * m
    else:
        return
    v = gcd(m1, 3)
    q = m1 + 1
    D = divisors(q // 4)

    for b in D:
        for c in D:
            if gcd(b, c) > 1:
                continue
            r = q // c
            r *= b
            if (r + 1) % v != 0:
                continue
            sf = m - r % m
            if sf in S:
                continue
            S.add(sf)


def filter1b(m, S):
    T = [1, 2, 3, 6]
    u = gcd(m, 3)
    if u == 3:
        T[2] = 0
    for t in T:
        v = 4 * t * u
        D1 = divisors(t * m)
        for a in D1:
            b = (t * m) // a
            if b < a:
                continue
            if gcd(a, b) > 1:
                continue
            D2 = divisors(a + b)
            for E in D2:
                if (E + 1) % v != 0:
                    continue
                sf = m - E % m
                if sf in S:
                    continue
                S.add(sf)


def filter1c(m, S):
    T = [1, 2, 3, 6]
    u = gcd(m, 3)
    if u == 3:
        T[2] = 0
    for t in T:
        v = 4 * t * u
        D = divisors(t * m)
        for E in D:
            if (E + 1) % 4 != 0:
                continue
            bd = (t * m) // E
            if gcd(bd, E) > 1:
                continue
            for b in D:
                if bd % b != 0:
                    continue
                s = E + 4 * b * bd
                if (s + 1) % v != 0:
                    continue
                sf = m - s % m
                if sf in S:
                    continue
                S.add(sf)


def filter2a(m, S):
    T = [1, 2, 3, 6]
    u = gcd(m, 3)
    if u == 3:
        T[2] = 0
    for t in T:
        v = 4 * t * u
        D1 = divisors(t * m)
        for a in D1:
            b = (t * m) // a
            if b < a:
                continue
            if gcd(a, b) > 1:
                continue
            D2 = divisors(a + b)
            for E in D2:
                if (E + 1) % v != 0:
                    continue
                sf = m - pow(E, -1, m)
                if sf in S:
                    continue
                S.add(sf)


def filter2b(m, S):
    T = [1, 2, 3, 6]
    R, M = [0, 0], [0, 0]
    u = gcd(m, 3)
    if u == 3:
        T[2] = 0
    for t in T:
        v = 4 * t * u
        D = divisors(t * m)
        for F in D:
            if (F + 1) % 4 != 0:
                continue
            bc = (t * m) // F
            if gcd(bc, F) > 1:
                continue
            w = gcd(4 * bc, 24)
            if (F + 1) % w != 0:
                continue
            for b in D:
                if bc % b != 0:
                    continue
                c = bc // b
                R[0] = c * pow(b, -1, F)
                R[0] %= F
                M[0] = F
                R[1] = F
                M[1] = 4 * bc
                s = chrem(R, M)
                if (s + 1) % v != 0:
                    continue
                sf = m - s % m
                if sf in S:
                    continue
                S.add(sf)


def filter2c(m, S):
    T = [1, 2, 3, 6]
    u = gcd(m, 3)
    if u == 3:
        T[2] = 0
    for t in T:
        v = 4 * t * u
        m1 = t * m
        D1 = divisors(m1)
        for b in D1:
            d = m1 // b
            D2 = divisors(4 * b * b * d + 1)
            for F in D2:
                if (F + 1) % v != 0:
                    continue
                sf = m - F % m
                if sf in S:
                    continue
                S.add(sf)


def filter2d(m, S):
    T = [1, 2, 3, 6]
    R, M = [0, 0], [0, 0]
    u = gcd(m, 3)
    if u == 3:
        T[2] = 0
    for t in T:
        v = 4 * t * u
        D = divisors(t * m)
        for F in D:
            if (F + 1) % 4 != 0:
                continue
            cd = (t * m) // F
            if gcd(cd, F) > 1:
                continue
            w1 = gcd(4 * cd, 24)
            if (F + 1) % w1 != 0:
                continue
            w2 = gcd(F, 24)
            cdf = cd % F
            for c in D:
                if cd % c != 0:
                    continue
                s2 = F - (4 * c * cdf) % F
                if (s2 - 1) % w2 != 0:
                    continue
                Q = modsqrt(s2, F)
                for q in Q:
                    R[0] = q
                    M[0] = F
                    R[1] = 4 * cd - F % (4 * cd)
                    M[1] = 4 * cd
                    s = chrem(R, M)
                    if (s - 1) % v != 0:
                        continue
                    sf = s % m
                    if sf in S:
                        continue
                    S.add(sf)



def filter(m):
    S = set()
    for filter_fn in [ filter1a, filter1b, filter1c, filter2a, filter2b, filter2c, filter2d ]:
        filter_fn(m, S)
    if isprime(m):
        S.add(0)
    return list(S)


def sieve(base, mod):
    filters = { m: filter(m) for m in base + mod }
    gap = 24
    for base_m in base:
        gap *= base_m
    for m in mod:
        if isprime(m):
            continue
        for d in mod:
            if m % d != 0:
                continue
            for a in range(m):
                if a % d in filters[d] and a not in filters[m]:
                    filters[m].append(a)
        q = gcd(m, gap)
        if q == 1:
            continue
        for d in base + mod:
            if q % d != 0:
                continue
            for a in range(m):
                if a % d in filters[d] and a in filters[m]:
                    filters[m].remove(a)
    return filters


# Add more primes to base to create R_8, R_9, ..., from them.
gap, base = [24], [5, 7, 11, 13, 17, 19, 23]

for m in base:
    gap.append(gap[-1] * m)

GAP = gap[-1]
DG = divisors(GAP // 24)

if 0 in DG: DG.remove(0)
if 1 in DG: DG.remove(1)
DG.sort()

kfs = { d: [] for d in DG }


def set_killing_filter(x):
    S = filter(x)
    C = []
    for r in S:
        old = False
        for m in DG:
            if x % m != 0:
                continue
            if r % m in kfs[m]:
                old = True
                break
        if not old:
            C.append(r)
    return C


def set_killing():
    for idx, x in enumerate(DG):
        print(f"KF done ({idx}/{len(DG)})")
        kfs[x] = set_killing_filter(x)


def gcdex(a, b, x):
    x0 = 1
    y0 = 0
    x1 = 0
    y1 = 1
    while b != 0:
        Q = a // b
        x0 -= Q * x1
        y0 -= Q * y1
        a -= Q * b
        if a != 0:
            Q = b // a
            x1 -= Q * x0
            y1 -= Q * y0
            b -= Q * a
        else:
            x[0] = x1
            x[1] = y1
            return b
    x[0] = x0
    x[1] = y0
    return a


def first_residues(row, R1):
    x = [0, 0]
    R = []
    m0 = base[row]
    m1 = gap[row]
    m = m0 * m1
    gcdex(m0, m1, x)
    if x[0] < 0: x[0] += m1
    if x[1] < 0: x[1] += m0
    t0 = m1 * x[1]
    t1 = m0 * x[0]
    for a in range(base[row]):
        if isprime(a) and a in filter(base[row]):
            continue
        r0 = (t0 * a) % m
        for b in R1:
            r1 = (t1 * b) % m
            r = (r0 + r1) % m
            R.append(r)
    return R


def reduction(row, R2):
    m = base[row]
    g = gap[row + 1] // 24
    R = []
    for r in R2:
        broke_early = False
        for x in DG:
            if g % x > 0 or x % m > 0:
                continue
            C = kfs[x]
            if r % x in C:
                broke_early = True
                break
        if not broke_early:
            R.append(r)
    return R
        

def set_residues():
    set_killing()
    R1 = [1]
    R2 = []
    for row in range(len(base)):
        R2 = first_residues(row, R1)
        R1 = reduction(row, R2)
        with open(f"R{row + 1}.txt", "w+") as file:
            print(f"R{row + 1} has size {len(R1)}")
            file.write(" ".join(list(map(str, R1))))


set_residues()
