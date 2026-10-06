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
for their signatures. This pyrforest change keeps its existing Python
interface; it adds no Python bindings for those ring APIs and makes no wrapper
algorithm changes.

Clone this wrapper change and initialize its pinned native submodule:

```
git clone --branch deps/rforest-pr3 https://github.com/jjgarzella/pyrforest.git
cd pyrforest
git submodule sync --recursive
git submodule update --init --recursive
git -C pyrforest/lib rev-parse HEAD
```

The final command prints the native source commit listed above. Install from the
checkout with Sage:

```
sage -pip install --no-build-isolation --upgrade .
```

To install for your user account, add `--user` to the install command.

## Development

Use the pinned checkout steps above, then install in editable mode:

```
make install
```

Run tests:

```
make test
```

## Functions provided
### Interface with rforest
#### remainder_forest(M, m, k, kbase=0, indices=None, V=None, ans=None, kappa=None, projective=False)
#### remainder_forest_generic_prime(M, d, e, k, indices=None, m=None, kbase=0, V=None, ans=None, kappa=None)

### Examples of usage

- `batch_factorial(n, e, gamma)`

  Return a dict whose value at a prime $p$ equals $(\lceil \gamma p \rceil-1)! \pmod p^{e}$.

- `batch_harmonic(n, e, gamma, j, proj=False)`

  Return a dict whose value at a prime $p$ is the truncated harmonic sum

$$ \sum_{k=1}^{\lceil \gamma p \rceil-1} k^{-j} mod p^e. $$

  If `proj` is True, instead return a $1 \times 2$ matrix $[x, y]$ representing $x/y \pmod{p^e}$.
