"""
Are Chebyshev polynomials irreducible over Q?
"""

from sympy import chebyshevt, Symbol, factor, Poly, QQ

x = Symbol('x')

print("FACTORISATION OF CHEBYSHEV T_n(x) OVER Q")
print("=" * 65)

for n in range(1, 21):
    T = chebyshevt(n, x)
    f = factor(T, domain=QQ)
    poly = Poly(T, x, domain=QQ)
    factors = poly.factor_list()
    # factors = (content, [(factor1, mult1), (factor2, mult2), ...])
    content, factor_pairs = factors
    num_factors = len(factor_pairs)
    degrees = [Poly(fi, x).degree() for fi, _ in factor_pairs]

    irr = "IRREDUCIBLE" if num_factors == 1 else f"REDUCIBLE → {len(factor_pairs)} factors"
    print(f"\n  T_{n:>2}(x) = {T}")
    print(f"         {irr}  (degrees: {degrees})")
    if num_factors > 1:
        print(f"         = {f}")
