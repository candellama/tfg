import hashlib
import os
from typing import Tuple, List

class MLDSA:
# 1- Initialization and parameters
    # function to initialize parameters based on security level
    def __init__(self, security_level: int = 2):
        self.security_level = security_level
        self.set_parameters()
    
    # function to set parameters based on security level
    def set_parameters(self):
        params = {
            2: {"k": 4, "l": 4, "d": 13, "eta": 2, "gamma1": 16, "gamma2": 8, "beta": 1},
            3: {"k": 6, "l": 5, "d": 13, "eta": 4, "gamma1": 16, "gamma2": 8, "beta": 1},
            5: {"k": 8, "l": 7, "d": 13, "eta": 2, "gamma1": 16, "gamma2": 8, "beta": 1}
        }
        
        p = params.get(self.security_level, params[2])
        self.k = p["k"]  # number of rows in A matrix
        self.l = p["l"]  # number of columns in A matrix
        self.d = p["d"]  # number of dropped bits
        self.eta = p["eta"]  # eta parameter for key sampling
        self.n = 16 # real is 256 ring dimension
        self.q =  257 #8380417 modulus for coefficients
    
        # Toy signing bounds
        self.gamma1 = p["gamma1"]
        self.gamma2 = p["gamma2"]
        self.beta = p["beta"]

# 2- Seed expansion
    # function to generate random seed for key generation
    # default seed size is 32 bytes (256 bits)
    def generate_seed(self) -> bytes:
        return os.urandom(32)
    
    # function to expand seed using SHAKE256 XOF into a pseudorandom byte stream of some length
    def expand_seed(self, seed: bytes, nonce: int, length: int) -> bytes:
        h = hashlib.shake_256()
        h.update(seed)
        h.update(nonce.to_bytes(2, "little"))  # robust (0..65535)
        return h.digest(length)
    
# 3- Sampling

    # function to obtain uniform matrix A k x l from seed rho
    def sample_uniform_matrix(self, rho: bytes) -> list:
        A = []
        for i in range(self.k):
            row = []
            for j in range(self.l):
                nonce = i * self.l + j
                # sample polynomial coefficients
                # chain of bytes representing the polynomial
                poly_bytes = self.expand_seed(rho, nonce, 4 * self.n)
                # Convert the bytes to polynomial coefficients mod q
                poly = self._bytes_to_coefficients_uniform_mod_q(poly_bytes)
                row.append(poly)
            A.append(row)
        return A
    
    # function to sample secret small polynomials s1 and s2, coefficients in [-eta, eta]
    def sample_bounded(self, sigma: bytes, nonce: int, eta: int, n: int) -> list:
        poly_bytes = self.expand_seed(sigma, nonce, n)
        coeffs = []
        for b in poly_bytes:
            coeffs.append((b % (2 * eta + 1)) - eta)
        return coeffs

    def _bytes_to_coefficients_uniform_mod_q(self, data: bytes) -> list:
        coeffs = []
        needed = 4 * self.n
        if len(data) < needed:
            raise ValueError(f"Need at least {needed} bytes, got {len(data)}")

        for i in range(self.n):
            chunk = data[4*i : 4*i + 4]
            val = int.from_bytes(chunk, "little") % self.q
            coeffs.append(val)
        return coeffs
    
# 4- Algebraic and ring arithmetic operations 
    def _center(self, x: int) -> int:
        x %= self.q
        if x > self.q // 2:
            return x - self.q
        return x

    def _centered_mod(self, x: int, m: int) -> int:
        r = x % m
        return r - m if r > m // 2 else r

    def _poly_norm_inf(self, a: List[int]) -> int:
        return max(abs(self._center(x)) for x in a)

    def _vec_norm_inf(self, v: List[List[int]]) -> int:
        return max(self._poly_norm_inf(p) for p in v)
    

    # function to do: poly1 + scalar*poly2 of two polys expressed as lists of their polynomials
    def _poly_scalar_add(self, poly1: List[int], scalar: int, poly2: List[int]) -> List[int]:
        if scalar == 0:
            return [c1 % self.q for c1 in poly1]
        return [(c1 + scalar*c2) % self.q for c1, c2 in zip(poly1, poly2)]
    
    def _poly_add(self, a: List[int], b: List[int]) -> List[int]:
        return self._poly_scalar_add(a, 1, b)

    def _poly_sub(self, a: List[int], b: List[int]) -> List[int]:
        return self._poly_scalar_add(a, -1, b)

    # the same function expanded to vectors of polys
    # vec1+scalar*vec2
    def _vec_scalar_add(self, vec1: List[List[int]], scalar: int, vec2: List[List[int]]) -> List[List[int]]:
        return [self._poly_scalar_add(v, scalar, a) for v, a in zip(vec1, vec2)]

    # function to multiply two polynomials in R_q = Z_q[x]/(x^n + 1)
    def _poly_mul_negacyclic(self, a: List[int], b: List[int]) -> List[int]:
        """
        Multiply in R_q = Z_q[x]/(x^n + 1)
        """
        n = self.n
        q = self.q
        res = [0] * n
        for i in range(n):
            for j in range(n):
                prod = a[i] * b[j]
                idx = i + j
                # If idx >= n, wrap around and subtract because x^n == -1.
                if idx < n:
                    res[idx] = (res[idx] + prod) % q
                else:
                    res[idx - n] = (res[idx - n] - prod) % q
        return res


    
    def _matrix_mult_add(self, A: List[List[List[int]]], s1: List[List[int]], s2: List[List[int]]) -> List[List[int]]:
        """
        Compute t = A*s1 + s2 over the ring R_q (polynomial arithmetic).
        A: k x l matrix of polynomials (each polinomial is length-n list mod q)
        s1: length-l vector of polynomials
        s2: length-k vector of polynomials (small coefficients)
        """
        t = []
        for i in range(self.k):
            acc = [0] * self.n
            for j in range(self.l):
                acc = self._poly_add(acc, self._poly_mul_negacyclic(A[i][j], s1[j]))
            # add s2 (convert negatives to mod q if needed)
            acc = [(acc[k] + (s2[i][k] % self.q)) % self.q for k in range(self.n)]
            t.append(acc)
        return t
    
    # function to do matrix*vec of polynomials
    def _mat_vec_mul(self, A, v):
        zero_vec = [[0] * self.n for _ in range(self.k)]
        return self._matrix_mult_add(A, v, zero_vec)
    
    
# 5- highbits and lowbits

    # function to compute high and lowbits of a polynomial expressed as list of its coeffs
    def _highbits_poly(self, poly: List[int], gamma2: int) -> Tuple[List[int], List[int]]:
        highs, lows = [], []
        for c in poly:
            cc = self._center(c)
            # low  = coeff mod gamma2
            low = self._centered_mod(cc, gamma2)       # symmetric remainder
            # high = (coeff - low)/gamma2
            high = (cc - low) // gamma2           # integer exact
            highs.append(high)
            lows.append(low)
        return highs, lows

    # extended function to compute the high and lowbits of a vector of polys
    def _highbits_vec(self, vec: List[List[int]], gamma2: int) -> Tuple[List[List[int]], List[List[int]]]:
        highs, lows = [], []
        for p in vec:
            h, l = self._highbits_poly(p, gamma2)
            highs.append(h)
            lows.append(l)
        return highs, lows

    def _lowbits_vec(self, vec: List[List[int]], gamma2: int) -> List[List[int]]:
        return [self._highbits_poly(p, gamma2)[1] for p in vec]

    # this function returns either 0 or 1 so it is in the B_tau although -1 is not used for simplicity in the toy
    def _hash_to_bit(self, message: bytes, w1: List[List[int]]) -> int:
        # Toy challenge: 1 bit from SHA-256.
        h = hashlib.sha256()
        h.update(message)   # this is H(M)
        # Serialize w1 as signed-ish small ints, we mod to 16-bit for stable bytes
        for poly in w1:
            for v in poly:
                h.update(int(v & 0xFFFF).to_bytes(2, "little")) # little endian
        return h.digest()[0] & 1  

#########################################################################################################
#########################################################################################################
# Toy version of Dilithium algorithms 
    # 1 KEYPAIR GENERATION FUNCTION 
    def generate_keypair(self, seed: bytes = None) -> Tuple[Tuple, Tuple]:

        if seed is None:
            seed = self.generate_seed()
        
        # Seed expansion
        h = hashlib.shake_256()
        h.update(b"KeyGen" + seed)
        expanded = h.digest(128)
        
        rho = expanded[:32]  # seed for A
        sigma = expanded[32:64]  # seed for s1 and s2
        
        # Sample matrix A
        A = self.sample_uniform_matrix(rho)
        
        # Sample secretes s1 and s2
        s1 = []
        s2 = []
        
        for i in range(self.l):
            nonce = i
            s1_poly = self.sample_bounded(sigma, nonce, self.eta, self.n)
            s1.append(s1_poly)
        
        for i in range(self.k):
            nonce = self.l + i
            s2_poly = self.sample_bounded(sigma, nonce, self.eta, self.n)
            s2.append(s2_poly)
        
        # Compute t = A * s1 + s2 (simplified matrix multiplication)
        t = self._matrix_mult_add(A, s1, s2)
        
        # Return raw components without encoding
        public_key = (rho, t)
        private_key = (rho, s1, s2, t)
        
        return public_key, private_key
    

    ############################################
    # 2 SIGNATURE GENERATION FUNCTION
    def sign(self, message: bytes, private_key: tuple, max_tries: int = 10000) -> Tuple[int, List[List[int]]]:
        rho, s1, s2, _t = private_key
        A = self.sample_uniform_matrix(rho)

        for attempt in range(max_tries):
            # 1) sample y with small coefficients in [-gamma1, gamma1]
            y = [self.sample_bounded(os.urandom(32), i, self.gamma1, self.n) for i in range(self.l)]

            # 2) w = A*y
            w = self._mat_vec_mul(A, y)

            # 3) w1 = HighBits(w, gamma2)
            w1, _w0 = self._highbits_vec(w, self.gamma2)

            # 4) c = H(m || w1)  (toy: c in {0,1})
            c = self._hash_to_bit(message, w1)

            # 5) z = y + c*s1
            z = self._vec_scalar_add( y, c, s1)

            # 6) rejection sampling 
            if self._vec_norm_inf(z) > (self.gamma1 - self.beta):
                continue

            # w0' = LowBits(w - c*s2, gamma2)
            if c == 0:
                w_minus = w
            else:
                # subtract s2 from each component (mod q)
                w_minus = [self._poly_sub(w[i], [c_ % self.q for c_ in s2[i]]) for i in range(self.k)]
            w0_prime = self._lowbits_vec(w_minus, self.gamma2)

            # ||w0'||_inf <= gamma2 - beta
            if max(max(abs(x) for x in poly) for poly in w0_prime) > (self.gamma2 - self.beta):
                continue

            return c, z, attempt

        raise RuntimeError(f"Signing failed after {max_tries} attempts (try increasing gamma1/gamma2 or lowering beta).")


    #############################################
    # 3 SIGNATURE VERIFICATION FUNCTION
    def verify(self, message: bytes, signature: Tuple[int, List[List[int]]], public_key: tuple) -> bool:
        c, z, a = signature
        rho, t = public_key
        A = self.sample_uniform_matrix(rho)

        if c not in (0, 1):
            return False

        if self._vec_norm_inf(z) > (self.gamma1 - self.beta):
            return False

        # w1' = HighBits(Az - c t)
        Az = self._mat_vec_mul(A, z)
        if c == 0:
            Az_minus_ct = Az
        else:
            Az_minus_ct = [self._poly_sub(Az[i], t[i]) for i in range(self.k)]

        w1_prime, _ = self._highbits_vec(Az_minus_ct, self.gamma2)
        c_prime = self._hash_to_bit(message, w1_prime)

        return c_prime == c
