import hashlib
import os
from typing import List, Tuple
from toy_key_gen import MLDSA

# Toy signing steps
#   w = A*y
#   w1 = HighBits(w, gamma2)
#   c = H(m || w1)   (toy: c is a single bit 0/1)
#   z = y + c*s1
#   reject if norms are too large (toy rejection sampling)

########### highbits and lowbits functions

# function to compute high and lowbits of a polynomial expressed as list of its coeffs
def _highbits_poly(mldsa: MLDSA, poly: List[int], gamma2: int) -> Tuple[List[int], List[int]]:
    highs, lows = [], []
    for c in poly:
        cc = mldsa._center(c,mldsa.q)
        # low  = coeff mod gamma2
        low = mldsa._centered_mod(cc, gamma2)       # symmetric remainder
        # high = (coeff - low)/gamma2
        high = (cc - low) // gamma2           # integer exact
        highs.append(high)
        lows.append(low)
    return highs, lows

# function to compute the high and lowbits of a vector of polys
def _highbits_vec(mldsa: MLDSA, vec: List[List[int]], gamma2: int) -> Tuple[List[List[int]], List[List[int]]]:
    highs, lows = [], []
    for p in vec:
        h, l = _highbits_poly(mldsa, p, gamma2)
        highs.append(h)
        lows.append(l)
    return highs, lows

def _lowbits_vec(mldsa: MLDSA, vec: List[List[int]], gamma2: int) -> List[List[int]]:
    return [_highbits_poly(mldsa, p, gamma2)[1] for p in vec]



########################3
# this function returns either 0 or 1 so it is in the B_tau although -1 is not used for simplicity in the toy
def _hash_to_bit(message: bytes, w1: List[List[int]]) -> int:
    # Toy challenge: 1 bit from SHA-256.
    h = hashlib.sha256()
    h.update(message)   # this is H(M)
    # Serialize w1 as signed-ish small ints, we mod to 16-bit for stable bytes
    for poly in w1:
        for v in poly:
            h.update(int(v & 0xFFFF).to_bytes(2, "little")) # little endian
    return h.digest()[0] & 1  

#####   functions for algebraic ops
# function to do: poly1 + scalar*poly2 of two polys expressed as lists of their polynomials
def _poly_scalar_add(mod: int, poly1: List[int], scalar: int, poly2: List[int]) -> List[int]:
    if scalar == 0:
        return [c1 % mod for c1 in poly1]
    return [(c1 + scalar*c2) % mod for c1, c2 in zip(poly1, poly2)]

# the same operation expanded to vectors of polys
# vec1+scalar*vec2
def _vec_scalar_add(mod: int, vec1: List[List[int]], scalar: int, vec2: List[List[int]]) -> List[List[int]]:
    return [_poly_scalar_add(mod, v, scalar, a) for v, a in zip(vec1, vec2)]

# function to do matrix*vec of polynomials
def _mat_vec_mul(mldsa: MLDSA, A: List[List[List[int]]], v: List[List[int]]) -> List[List[int]]:
    # Returns k polynomials.
    out = []
    for i in range(mldsa.k):
        acc = [0] * mldsa.n  # vector of n 0s
        for j in range(mldsa.l):
            acc = mldsa._poly_add(acc, mldsa._poly_mul_negacyclic(A[i][j], v[j]))
        out.append(acc)
    return out



##############################33
# signature generation
def sign(mldsa: MLDSA, message: bytes, private_key: tuple, max_tries: int = 1000) -> Tuple[int, List[List[int]]]:
    rho, s1, s2, _t = private_key
    A = mldsa.sample_uniform_matrix(rho)

    for _ in range(max_tries):
        # 1) sample y with small coefficients in [-gamma1, gamma1]
        y = [mldsa.sample_bounded(os.urandom(32), i, mldsa.gamma1, mldsa.n) for i in range(mldsa.l)]

        # 2) w = A*y
        w = _mat_vec_mul(mldsa, A, y)

        # 3) w1 = HighBits(w, gamma2)
        w1, _w0 = _highbits_vec(mldsa, w, mldsa.gamma2)

        # 4) c = H(m || w1)  (toy: c in {0,1})
        c = _hash_to_bit(message, w1)

        # 5) z = y + c*s1
        z = _vec_scalar_add(mldsa.q, y, c, s1)

        # 6) rejection sampling 
        if mldsa._vec_norm_inf(z) > (mldsa.gamma1 - mldsa.beta):
            continue

        # w0' = LowBits(w - c*s2, gamma2)
        if c == 0:
            w_minus = w
        else:
            # subtract s2 from each component (mod q)
            w_minus = [mldsa._poly_sub(w[i], [c_ % mldsa.q for c_ in s2[i]]) for i in range(mldsa.k)]
        w0_prime = _lowbits_vec(mldsa, w_minus, mldsa.gamma2)

        # ||w0'||_inf <= gamma2 - beta
        if max(max(abs(x) for x in poly) for poly in w0_prime) > (mldsa.gamma2 - mldsa.beta):
            continue

        return c, z

    raise RuntimeError(f"Signing failed after {max_tries} attempts (try increasing gamma1/gamma2 or lowering beta).")

def verify(mldsa: MLDSA, message: bytes, signature: Tuple[int, List[List[int]]], public_key: tuple) -> bool:
    c, z = signature
    rho, t = public_key
    A = mldsa.sample_uniform_matrix(rho)

    if c not in (0, 1):
        return False

    if mldsa._vec_norm_inf(z) > (mldsa.gamma1 - mldsa.beta):
        return False

    # w1' = HighBits(Az - c t)
    Az = _mat_vec_mul(mldsa, A, z)
    if c == 0:
        Az_minus_ct = Az
    else:
        Az_minus_ct = [mldsa._poly_sub(Az[i], t[i]) for i in range(mldsa.k)]

    w1_prime, _ = _highbits_vec(mldsa, Az_minus_ct, mldsa.gamma2)
    c_prime = _hash_to_bit(message, w1_prime)

    return c_prime == c

if __name__ == "__main__":
    m = MLDSA(security_level=2)
    pk, sk = m.generate_keypair()

    msg = b"Hello, world!"
    sig = sign(m, msg, sk)

    print("Signature c:", sig[0])
    print("Signature z[0][:8]:", sig[1][0][:8])
    print("Verify:", verify(m, msg, sig, pk))
