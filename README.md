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
from pyrforest import (
    remainder_forest_p2,
    remainder_forest_p2_packed,
    remainder_forest_pn,
    remainder_forest_pn_native,
)
```

Run the wrapper doctests, independent ring-binding references, and native
checks from the checkout:

```
python -m sage.doctest --force-lib pyrforest
python tests/test_ring_forest_bindings.py
python tests/test_build_native_pn_consumer.py
make -C pyrforest/lib check
```

To verify the separately compiled PN consumer against the files actually
included in a wheel, build a wheel and pass it to the retained consumer test.
The helper extracts the wheel, compiles against its header and Cython
declaration, then runs the producer and consumer in an isolated child process.
It prints the loaded package, compiled extension, consumer, header, and
declaration paths with SHA-256 hashes, and exercises a conflicting ambient
`PYTHONPATH` sentinel so a source checkout or older installed producer cannot
make the check pass accidentally:

```
python setup.py bdist_wheel --dist-dir build/dist
sage -python tests/build_native_pn_consumer.py build/dist/pyrforest-*.whl
```

## Functions provided
### Interface with rforest
#### remainder_forest(M, m, k, kbase=0, indices=None, V=None, ans=None, kappa=None, projective=False)
#### remainder_forest_p2(M, m, k, kbase=0, indices=None, V=None, kappa=None, return_state=False, *, z=None)
#### remainder_forest_p2_packed(M, primes, kbase=0, indices=None, V=None, kappa=None)
#### remainder_forest_pn(M, m, k, nP, kbase=0, indices=None, V=None, kappa=None, return_state=False, *, z=None)
#### remainder_forest_pn_native(M, m, k, nP, kbase=0, indices=None, V=None, kappa=None, return_state=False, *, z=None)
#### remainder_forest_generic_prime(M, d, e, k, indices=None, m=None, kbase=0, V=None, ans=None, kappa=None)

The ring functions accept an exact Sage matrix over a two-generator polynomial
ring over `ZZ`. Its first generator is formal `P` and its second is transition
variable `x`, regardless of the names chosen for those generators. `x` is
evaluated at each transition index; `P` remains formal.
`remainder_forest_p2` retains the coefficients of `1` and `P`.
`remainder_forest_pn` takes an explicit positive `nP` and retains powers below
`P^nP` (the powers `0` through `nP-1`). Coefficients at higher P powers are
truncated before the native call.

`remainder_forest_p2_packed` is an explicit opt-in output path for prime
moduli `p^2`, with endpoint `k=p`. It accepts primes in nondecreasing order
from `2` through `min(2^32-1, LONG_MAX)`, where `LONG_MAX` is the largest
signed value of the host's native C `long`. A 64-bit `long` supports the full
uint32 prime range; a narrower `long` caps the packed endpoint domain at its
representable maximum. The function rejects composites and unsupported
endpoints with `PackedP2UnsupportedDomain` before passing any endpoint to C.
Catch that exception when a consumer needs to fall back to the
arbitrary-precision `remainder_forest_p2` interface. The ordinary Sage matrix
functions keep their existing result type and behavior.

The packed result has logical shape `(batch, 2, rows, columns)`: for each
prime, plane 0 contains `P^0`, plane 1 contains `P^1`, and each matrix plane is
row-major. `result.primes`, `result.endpoints`, `result.moduli`, and
`result.indices` retain the exact input order. `result.coefficients` is a
read-only flat memoryview with format `Q` (native-endian uint64); the linear
index is `(((batch_position*2 + plane)*rows + row)*columns + column)`.
`result.matrix(batch_position, plane)` returns one shaped read-only matrix
view, and `result.buffer` exposes the immutable bytes that own the copy.
Every coefficient is checked to lie in `[0, p^2)` before export.
`PackedP2ForestResult` snapshots metadata sequences as tuples. Direct
construction validates aligned metadata, the supported prime/endpoint/modulus
relationship, positive matrix dimensions, canonical coefficient words, and
shape/storage limits before it exposes a view. Total packed words are bounded
by native indexing and Python buffer limits; each batch plane remains C-int
bounded.

`remainder_forest_pn_native` is a separate opt-in arbitrary-precision output
path. It retains the same general PN inputs as `remainder_forest_pn`, including
independent formal precision `nP`, arbitrary positive GMP moduli, caller
indices, and native C `long` endpoints. Its immutable
`NativePnForestResult` reports `shape == (batch, nP, rows, columns)` plus
separate tuples for `indices`, `endpoints`, and `moduli`. Endpoints are the
actual requested forest endpoints; they are not inferred from the moduli and
need not equal any prime. Coefficients stay in producer-owned GMP storage and
are not converted to per-coefficient Python or Sage integers. `return_state`
is unsupported on this result path; use `remainder_forest_pn` when the final
`V` state is needed.

Compiled consumers call `result.as_capsule()` and keep both the result and
capsule alive for the complete native call. Each call returns a fresh capsule;
the capsule holds a strong reference to the result, so its pointers remain
valid even if the caller drops its direct result reference. Pointers must not
be retained after the capsule/owner lifetime or used to modify the source
GMP values. The fixed name is
`pyrforest.PnGmpForestResult/v1`; `pyrforest/native_pn.h` and
`pyrforest/native_pn.pxd` define the matching versioned ABI. The C-contiguous
coefficient order is `(batch, degree, row, column)`, with offset
`(((batch*nP + degree)*rows + row)*columns + column)`. Each coefficient is
canonical in `[0, moduli[batch])`. The payload records exact GMP version and
limb sizes so consumers can reject incompatible builds before reading data.

The public feature markers are `PN_GMP_RESULT_ABI_VERSION == 1`,
`PN_GMP_RESULT_CAPABILITY == "pn-gmp-result-v1"`, and
`PN_GMP_RESULT_CAPSULE_NAME == "pyrforest.PnGmpForestResult/v1"`. Wheels and
source distributions include the ABI header and Cython declaration. A
compiled client should cimport `pyrforest.native_pn` and add the installed
`pyrforest` package directory to its C compiler include path.

For Trace inputs, the same function can be used for each matrix and the scalar
call, with `V` omitted for the identity matrix:

```python
packed = remainder_forest_p2_packed(M, primes, kbase=1)
trace_plane_0 = packed.matrix(0, 0)
trace_plane_1 = packed.matrix(0, 1)
```

The buffer remains valid after the native forest temporaries are released and
after later calls. The export traverses native GMP coefficients in compiled
code and does not construct intermediate Sage result matrices or polynomial
coefficients.

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
