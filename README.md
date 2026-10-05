# pyrforest

This package is a simple wrapper for the [rforest](https://github.com/jjgarzella/rforest) library code into SageMath.

## Install

The tested wrapper source is pinned to commit
`08a4ebbbca9008ae7478d4ec541fa026034d76b6`. It pins the native `rforest`
source to `4bad17781bb3e0842753b1d16fd8a2bf36af99cd`.

Clone that wrapper commit and initialize its pinned native submodule:

```
git clone https://github.com/jjgarzella/pyrforest.git
cd pyrforest
git checkout 08a4ebbbca9008ae7478d4ec541fa026034d76b6
git submodule sync --recursive
git submodule update --init --recursive
git -C lib rev-parse HEAD
```

The final command prints the native source commit listed above. Install from
the checkout with Sage:

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
