# pyrforest

This package is a simple wrapper for the [rforest](https://github.com/jjgarzella/rforest) library code into SageMath.

## Install

The native `rforest` source is pinned at
`9a0a53f6b1a12cab4a845d181a4207905307b6fa`, the exact head of [rforest PR
#6](https://github.com/jjgarzella/rforest/pull/6) on
`features/ring-remainder-forest`. The submodule source repository remains
`https://github.com/jjgarzella/rforest.git`.

This native revision retains the integer-polynomial remainder forest and adds
the C APIs `rforest_p2`, `rforest_pn`, and `rforest_pnq` for matrix forests over
`Z[P]/(P^2)`, `Z[P]/(P^n)`, and `Z[P,Q]/(P^N,Q^N)`. See the pinned native
[`rforest.h`](https://github.com/jjgarzella/rforest/blob/9a0a53f6b1a12cab4a845d181a4207905307b6fa/rforest.h)
for their signatures. This wrapper exposes the univariate `P²` and `Pⁿ`
facades; it does not expose the bivariate `P^N,Q^N` facade.

Clone this wrapper change and initialize its pinned native submodule:

```
git clone --branch main https://github.com/jjgarzella/pyrforest.git
cd pyrforest
git submodule sync --recursive
git submodule update --init --recursive
git -C pyrforest/lib rev-parse HEAD
```

The final command prints the native source commit listed above. Install from the
checkout with the Python environment that has SageMath available:

```
python -m pip install --no-build-isolation --upgrade .
```

To install for your user account, add `--user` to the install command.

## Development

Use the pinned checkout steps above, then install from the Python environment
that has SageMath available:

```
python -m pip install --no-build-isolation -e .
```

Import Sage before importing `pyrforest`:

```
from sage.all import Matrix, PolynomialRing, ZZ
from pyrforest import remainder_forest_p2, remainder_forest_pn
```

Run the wrapper doctests, independent ring-binding references, and native
checks from the checkout:

```
python -m sage.doctest --force-lib pyrforest
python tests/test_ring_forest_bindings.py
make -C pyrforest/lib check
```

## Functions provided
### Interface with rforest
#### remainder_forest(M, m, k, kbase=0, indices=None, V=None, ans=None, kappa=None, projective=False)
#### remainder_forest_p2(M, m, k, kbase=0, indices=None, V=None, kappa=None, return_state=False, *, z=None)
#### remainder_forest_pn(M, m, k, nP, kbase=0, indices=None, V=None, kappa=None, return_state=False, *, z=None)
#### remainder_forest_generic_prime(M, d, e, k, indices=None, m=None, kbase=0, V=None, ans=None, kappa=None)

The ring functions accept an exact Sage matrix over a two-generator polynomial
ring over `ZZ`. Its first generator is formal `P` and its second is transition
variable `x`, regardless of the names chosen for those generators. `x` is
evaluated at each transition index; `P` remains formal.
`remainder_forest_p2` retains the coefficients of `1` and `P`.
`remainder_forest_pn` takes an explicit positive `nP` and retains powers below
`P^nP` (the powers `0` through `nP-1`). Coefficients at higher P powers are
truncated before the native call.

The initial `V` may be rectangular. It can use exact integer constants, a
one-generator integer polynomial ring (whose generator is formal `P`), or
the same two-generator ring as `M`; its entries must be independent of `x`.
If omitted, it is the identity. The native arrays use
`M[row][column][P exponent][ascending x degree]`,
`V[P exponent][row][column]`, and output
`[endpoint][P exponent][row][column]`. Results are a dict keyed by `indices`
(or `range(len(m))`); each value is a Sage matrix over
`(ZZ/m[i]ZZ)[P]` with degree below `nP`. Endpoints are exclusive, and
`kbase`, `indices`, and `kappa` follow the integer `remainder_forest`
conventions. Moduli must be positive integers and endpoints must be
nondecreasing and at least `kbase`.

With `return_state=True`, a ring function returns `(results, state)`. The state
keys are `z` and `final_V`: `z` is the exact residual native modulus, and
`final_V` is the final rows-by-dim coefficient-ring matrix represented over
`(ZZ/zZZ)[P]` with powers below `nP`. The keyword-only initial `z` must be a
positive multiple of all endpoint moduli; it is copied and never modified.
When omitted, it is their product, so the residual `z` is one and `final_V`
is zero in the modulus-one ring. An explicit cofactor lets callers inspect a
nontrivial final state. An empty batch returns the initial `V` reduced modulo
the initial `z`. Without `return_state`, the return value is just the result
dict, matching the usual `remainder_forest` result shape.

### Examples of usage

- `batch_factorial(n, e, gamma)`

  Return a dict whose value at a prime $p$ equals $(\lceil \gamma p \rceil-1)! \pmod p^{e}$.

- `batch_harmonic(n, e, gamma, j, proj=False)`

  Return a dict whose value at a prime $p$ is the truncated harmonic sum

$$ \sum_{k=1}^{\lceil \gamma p \rceil-1} k^{-j} mod p^e. $$

  If `proj` is True, instead return a $1 \times 2$ matrix $[x, y]$ representing $x/y \pmod{p^e}$.
