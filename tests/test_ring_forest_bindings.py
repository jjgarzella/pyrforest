"""Independent Sage references for the native P² and P^n wrappers."""

import gc
from pathlib import Path

from sage.all import Integers, LaurentPolynomialRing, Matrix, PolynomialRing, QQ, ZZ, prod

from pyrforest import (
    PACKED_P2_UINT64_CAPABILITY,
    PackedP2ForestResult,
    PackedP2UnsupportedDomain,
    remainder_forest,
    remainder_forest_p2,
    remainder_forest_p2_packed,
    remainder_forest_pn,
)
from pyrforest.rforest import deflate_matrix, inflate_matrix
import pyrforest.rforest as binding_module


def _terms(entry):
    if entry.parent() == ZZ:
        return {(): entry}
    result = {}
    for exponents, coefficient in entry.dict().items():
        try:
            exponents = tuple(exponents)
        except TypeError:
            exponents = (exponents,)
        result[exponents] = coefficient
    return result


def _ring_entry(entry, x_value, nP, ring):
    result = ring.zero()
    p = ring.gen()
    base = ring.base_ring()
    variable_count = len(entry.parent().variable_names()) if entry.parent() != ZZ else 0
    for exponents, coefficient in _terms(entry).items():
        p_exponent = exponents[0] if variable_count else 0
        x_degree = exponents[1] if variable_count == 2 else 0
        if p_exponent < nP:
            result += base(coefficient) * p**p_exponent * base(x_value)**x_degree
    return result


def direct_ring_forest(M, moduli, endpoints, nP, kbase=0, V=None):
    """Sequential exact Sage coefficient-ring products, independent of rforest."""
    result = {}
    dim = M.nrows()
    rows = dim if V is None else V.nrows()
    p_name = M.base_ring().variable_names()[0]
    for index, modulus, endpoint in zip(range(len(moduli)), moduli, endpoints):
        ring = PolynomialRing(Integers(modulus), p_name)
        if V is None:
            current = Matrix.identity(ring, dim)
        else:
            current = Matrix(
                ring,
                rows,
                dim,
                [_ring_entry(V[row, col], 0, nP, ring)
                 for row in range(rows) for col in range(dim)],
            )
        for transition in range(kbase, endpoint):
            step = Matrix(
                ring,
                dim,
                dim,
                [_ring_entry(M[row, col], transition, nP, ring)
                 for row in range(dim) for col in range(dim)],
            )
            product = current * step
            current = Matrix(
                ring,
                rows,
                dim,
                [ring([entry[p] for p in range(nP)]) for entry in product.list()],
            )
        result[index] = current
    return result


def specialize_p_zero(M, V):
    """Build the legacy integer-polynomial matrices by taking P coefficient 0."""
    x_ring = PolynomialRing(ZZ, "x")
    x = x_ring.gen()

    def specialize(entry):
        result = x_ring.zero()
        variable_count = len(entry.parent().variable_names()) if entry.parent() != ZZ else 0
        for exponents, coefficient in _terms(entry).items():
            if not variable_count or exponents[0] == 0:
                x_degree = exponents[1] if variable_count == 2 else 0
                result += coefficient * x**x_degree
        return result

    matrix = Matrix(x_ring, M.nrows(), M.ncols(), [specialize(v) for v in M.list()])
    if V is None:
        return matrix, None
    initial = Matrix(ZZ, V.nrows(), V.ncols(), [specialize(v) for v in V.list()])
    return matrix, initial


def assert_raises(error_type, function):
    try:
        function()
    except error_type:
        return
    raise AssertionError("expected %s" % error_type.__name__)


def assert_packed_equals(packed, expected, primes, indices=None):
    if indices is None:
        indices = tuple(range(len(primes)))
    else:
        indices = tuple(indices)
    rows = next(iter(expected.values())).nrows()
    columns = next(iter(expected.values())).ncols()
    assert isinstance(packed, PackedP2ForestResult)
    assert packed.shape == (len(primes), 2, rows, columns)
    assert packed.indices == indices
    assert packed.primes == tuple(primes)
    assert packed.endpoints == tuple(primes)
    assert packed.moduli == tuple(p * p for p in primes)
    assert packed.plane_names == ("P^0", "P^1")
    flat = packed.coefficients
    assert flat.format == "Q" and flat.itemsize == 8
    assert flat.readonly and flat.c_contiguous
    assert len(flat) == len(primes) * 2 * rows * columns
    for batch_position, index in enumerate(indices):
        for plane in range(2):
            matrix = packed.matrix(batch_position, plane)
            assert tuple(matrix.shape) == (rows, columns)
            assert matrix.format == "Q" and matrix.readonly
            for row in range(rows):
                for column in range(columns):
                    expected_value = int(expected[index][row, column][plane])
                    assert matrix[row, column] == expected_value
                    offset = (((batch_position * 2 + plane) * rows + row)
                              * columns + column)
                    assert flat[offset] == expected_value


def main():
    assert binding_module.__file__.endswith((".so", ".pyd"))
    assert Path(binding_module.__file__).resolve().parent == \
        Path(__file__).resolve().parents[1] / "pyrforest"
    # The generator names are intentionally arbitrary. Their positions define
    # P first and x second, as specified by the bridge contract.
    S = PolynomialRing(ZZ, names=("u", "t"))
    P, x = S.gens()
    M = Matrix(
        S,
        [
            [1 + 2 * P - 3 * x + P * x**2 + 2 * P**3 + 7 * P**4 * x + 10**35,
             -7 + P + 4 * x - 2 * P**2 * x + 5 * P**4],
            [5 * P - 9 * x + P**2 * (3 + x) - 10**30,
             2 - 4 * P * x + 3 * x**2 + 11 * P**3 - P**4],
        ],
    )
    V = Matrix(S, [[2 - 3 * P + 5 * P**2 + 2 * P**4,
                    -4 + P + 6 * P**3]])
    V_ring = PolynomialRing(ZZ, "q")
    q = V_ring.gen()
    V_one_generator = Matrix(V_ring, [[2 - 3 * q + 5 * q**2 + 2 * q**4,
                                       -4 + q + 6 * q**3]])
    moduli = [5, 7, 11, 13]
    endpoints = [1, 3, 5, 7]
    kbase = 1
    M_before = M.__copy__()
    V_before = V.__copy__()
    moduli_before = list(moduli)
    endpoints_before = list(endpoints)

    # x and P remain distinct, and two successive transition matrices do not
    # commute for this fixture.
    R101 = PolynomialRing(Integers(101), "u")
    T0 = Matrix(R101, 2, 2, [_ring_entry(entry, 1, 5, R101) for entry in M.list()])
    T1 = Matrix(R101, 2, 2, [_ring_entry(entry, 2, 5, R101) for entry in M.list()])
    assert T0 * T1 != T1 * T0

    for nP in (1, 2, 3, 5):
        expected = direct_ring_forest(M, moduli, endpoints, nP, kbase, V)
        for kappa in (0, 1, 3):
            actual = remainder_forest_pn(
                M, moduli, endpoints, nP, kbase=kbase, V=V, kappa=kappa
            )
            assert actual == expected, (nP, kappa, actual, expected)
            assert list(actual) == list(range(len(moduli)))
            assert actual[0].base_ring().variable_names() == ("u",)
            assert actual[0] is not M and actual[0] is not V

    assert remainder_forest_pn(
        M, moduli, endpoints, 3, kbase=kbase, V=V_one_generator, kappa=1
    ) == remainder_forest_pn(M, moduli, endpoints, 3, kbase=kbase, V=V, kappa=1)

    keyed = remainder_forest_pn(
        M,
        {"first": 5, "second": 7},
        {"first": 3, "second": 5},
        2,
        kbase=kbase,
        indices=(key for key in ("first", "second")),
        V=V,
    )
    keyed_expected = direct_ring_forest(M, [5, 7], [3, 5], 2, kbase, V)
    assert keyed == {"first": keyed_expected[0], "second": keyed_expected[1]}
    callable_keyed = remainder_forest_pn(
        M,
        lambda key: {"first": 5, "second": 7}[key],
        lambda key: {"first": 3, "second": 5}[key],
        2,
        kbase=kbase,
        indices=(key for key in ("first", "second")),
        V=V,
    )
    assert callable_keyed == keyed

    expected_p2 = direct_ring_forest(M, moduli, endpoints, 2, kbase, V)
    for kappa in (0, 1, 3):
        p2 = remainder_forest_p2(
            M, moduli, endpoints, kbase=kbase, V=V, kappa=kappa
        )
        pn2 = remainder_forest_pn(
            M, moduli, endpoints, 2, kbase=kbase, V=V, kappa=kappa
        )
        assert p2 == expected_p2
        assert pn2 == p2

        # The historical P² block embedding remains a secondary reference.
        block_M = inflate_matrix(M, [P], 2)
        block_V = inflate_matrix(V, [P], 2)
        block_outputs = remainder_forest(
            block_M, moduli, endpoints, kbase=kbase, V=block_V, kappa=kappa
        )
        for index in range(len(moduli)):
            assert deflate_matrix(block_outputs[index], [P], 2) == p2[index]

    # PN precision one agrees with the old integer forest.
    integer_M, integer_V = specialize_p_zero(M, V)
    integer_outputs = remainder_forest(
        integer_M, moduli, endpoints, kbase=kbase, V=integer_V, kappa=1
    )
    pn1 = remainder_forest_pn(
        M, moduli, endpoints, 1, kbase=kbase, V=V, kappa=1
    )
    for index, modulus in enumerate(moduli):
        base = Integers(modulus)
        constant = Matrix(
            base,
            1,
            2,
            [entry[0] for entry in pn1[index].list()],
        )
        assert constant == integer_outputs[index].change_ring(base)

    # Verify exact final native V/z with a nontrivial residual modulus.
    z0 = ZZ(13) * prod(ZZ(m) for m in moduli)
    z_before = ZZ(z0)
    with_state = remainder_forest_pn(
        M, moduli, endpoints, 5, kbase=kbase, V=V, return_state=True, z=z0
    )
    state_results, state = with_state
    assert state_results == direct_ring_forest(M, moduli, endpoints, 5, kbase, V)
    assert state["z"] == 13
    expected_final = direct_ring_forest(M, [13], [endpoints[-1]], 5, kbase, V)[0]
    assert state["final_V"] == expected_final
    assert z0 == z_before
    default_state = remainder_forest_p2(
        M, moduli, endpoints, kbase=kbase, V=V, return_state=True
    )[1]
    assert default_state["z"] == 1
    assert all(entry == 0 for entry in default_state["final_V"].list())
    assert with_state[0][0] is not state["final_V"]

    modulus_one = remainder_forest_pn(M, [1], [3], 2, kbase=kbase, V=V)[0]
    assert all(entry == 0 for entry in modulus_one[0].list())

    # Empty batches avoid native calls and preserve the initial state modulo z.
    empty = remainder_forest_pn(M, [], [], 5, V=V, return_state=True, z=13)
    assert empty[0] == {}
    assert empty[1]["z"] == 13
    assert empty[1]["final_V"] == direct_ring_forest(M, [13], [0], 5, 0, V)[0]
    assert remainder_forest_p2(M, [], []) == {}

    # Signed offsets, repeated endpoints, empty prefixes, and ordered products.
    signed_kbase = -6
    signed_endpoints = [-6, -4, -4, -2]
    signed_moduli = [5, 7, 7, 13]
    signed_expected = direct_ring_forest(
        M, signed_moduli, signed_endpoints, 3, signed_kbase, V
    )
    signed_actual = remainder_forest_pn(
        M, signed_moduli, signed_endpoints, 3, kbase=signed_kbase, V=V, kappa=0
    )
    assert signed_actual == signed_expected
    assert signed_actual[0] == direct_ring_forest(M, [5], [-6], 3, -6, V)[0]
    assert signed_actual[1] == signed_actual[2]  # repeated endpoint

    shifted_start, shifted_end = -6, -2
    shifted_M = Matrix(
        S,
        M.nrows(),
        M.ncols(),
        [entry.subs({x: x + shifted_start}) for entry in M.list()],
    )
    original_interval = remainder_forest_pn(
        M, [17], [shifted_end], 3, kbase=shifted_start, V=V
    )[0]
    translated_interval = remainder_forest_pn(
        shifted_M, [17], [shifted_end - shifted_start], 3, kbase=0, V=V
    )[0]
    assert original_interval == translated_interval

    # Nonzero positive kbase and repeated endpoints are also covered.
    positive = remainder_forest_pn(
        M, [5, 7, 7], [3, 5, 5], 2, kbase=3, V=V
    )
    assert positive == direct_ring_forest(M, [5, 7, 7], [3, 5, 5], 2, 3, V)
    assert positive[1] == positive[2]

    # Input arrays are immutable and repeated calls are independent.
    assert M == M_before and V == V_before
    assert moduli == moduli_before and endpoints == endpoints_before
    assert remainder_forest_p2(M, moduli, endpoints, kbase=kbase, V=V) == \
        remainder_forest_p2(M, moduli, endpoints, kbase=kbase, V=V)

    # Scalar matrices use the native dimension-one path directly.
    scalar_M = Matrix(S, [[1 + P + x]])
    scalar_moduli = [5, 7]
    scalar_endpoints = [1, 2]
    for nP in (1, 2, 3):
        assert remainder_forest_pn(
            scalar_M, scalar_moduli, scalar_endpoints, nP
        ) == direct_ring_forest(scalar_M, scalar_moduli, scalar_endpoints, nP)
    assert remainder_forest_p2(
        scalar_M, scalar_moduli, scalar_endpoints
    ) == direct_ring_forest(scalar_M, scalar_moduli, scalar_endpoints, 2)

    zero_M = Matrix(S, 2, 2, [0, 0, 0, 0])
    assert remainder_forest_pn(
        zero_M, moduli, endpoints, 5, kbase=kbase, V=V
    ) == direct_ring_forest(zero_M, moduli, endpoints, 5, kbase, V)

    # Invalid calls fail before native access, and an output-marshalling error
    # after a native call still releases all owned GMP storage.
    assert_raises(ValueError, lambda: remainder_forest_pn(M, [5], [0], 2, kbase=1))
    assert_raises(ValueError, lambda: remainder_forest_pn(M, [5, 7], [3, 2], 2, kbase=1))
    assert_raises(ValueError, lambda: remainder_forest_pn(M, [0], [1], 2))
    assert_raises(ValueError, lambda: remainder_forest_pn(M, [5], [1], 0))
    assert_raises(ValueError, lambda: remainder_forest_pn(M, [5], [1], 2, z=7))
    assert_raises(ValueError, lambda: remainder_forest_pn(M, [5], [1], 2, z=-5))
    long_min = -(ZZ(2) ** 63)
    long_max = ZZ(2) ** 63 - 1
    assert_raises(
        OverflowError,
        lambda: remainder_forest_pn(M, [5], [long_max], 2, kbase=long_min),
    )
    assert_raises(
        OverflowError,
        lambda: remainder_forest_pn(M, [5], [1], 2, kbase=long_min - 1),
    )
    assert_raises(ValueError, lambda: remainder_forest_pn(Matrix(S, 1, 2, [1, 2]), [5], [1], 2))
    assert_raises(ValueError, lambda: remainder_forest_pn(M, [5], [1], 2, V=Matrix(ZZ, 1, 1, [1])))
    assert_raises(ValueError, lambda: remainder_forest_pn(M, [5], [1], 2, V=Matrix(S, 1, 2, [x, 0])))
    Q = PolynomialRing(QQ, names=("u", "t"))
    assert_raises(TypeError, lambda: remainder_forest_pn(Matrix(Q, 1, 1, [1]), [5], [1], 2))
    assert_raises(TypeError, lambda: remainder_forest_pn(M, [5], [1], 2, indices=[[]]))

    # Laurent inputs are rejected before marshalling: the ring bridge only
    # accepts ordinary polynomial exponents and must not drop negative terms.
    L = LaurentPolynomialRing(ZZ, names=("u", "t"))
    P_laurent, x_laurent = L.gens()
    negative_P = Matrix(L, 1, 1, [P_laurent**-1 + 1])
    negative_x = Matrix(L, 1, 1, [x_laurent**-1 + 1])
    assert_raises(TypeError, lambda: remainder_forest_pn(negative_P, [5], [1], 2))
    assert_raises(TypeError, lambda: remainder_forest_pn(negative_x, [5], [1], 2))
    negative_V_P = Matrix(L, 1, 2, [P_laurent**-1, 1])
    negative_V_x = Matrix(L, 1, 2, [x_laurent**-1, 0])
    assert_raises(TypeError, lambda: remainder_forest_pn(M, [5], [1], 2, V=negative_V_P))
    assert_raises(TypeError, lambda: remainder_forest_pn(M, [5], [1], 2, V=negative_V_x))

    assert remainder_forest_pn(M, [5], [1], 2) == direct_ring_forest(M, [5], [1], 2)

    # The explicit packed P² route preserves prime order and avoids Sage
    # result-ring/matrix construction while retaining the legacy output API.
    assert PACKED_P2_UINT64_CAPABILITY == 1
    packed_primes = [5, 7, 11, 13]
    packed_moduli = [prime**2 for prime in packed_primes]
    packed_expected = direct_ring_forest(
        M, packed_moduli, packed_primes, 2, kbase, V
    )
    for packed_kappa in (0, 1, 3):
        legacy = remainder_forest_p2(
            M, packed_moduli, packed_primes,
            kbase=kbase, V=V, kappa=packed_kappa,
        )
        packed = remainder_forest_p2_packed(
            M, packed_primes, kbase=kbase, V=V, kappa=packed_kappa
        )
        assert legacy == packed_expected
        assert_packed_equals(packed, packed_expected, packed_primes)

    keyed_primes = {"g2": 5, "g3": 7, "g4": 11}
    keyed_order = ("g2", "g3", "g4")
    keyed_reference = direct_ring_forest(
        M, [25, 49, 121], [5, 7, 11], 2, kbase, V
    )
    keyed_expected_packed = {
        key: keyed_reference[position]
        for position, key in enumerate(keyed_order)
    }
    keyed_packed = remainder_forest_p2_packed(
        M, keyed_primes, kbase=kbase, indices=keyed_order, V=V, kappa=1
    )
    assert_packed_equals(
        keyed_packed, keyed_expected_packed, [5, 7, 11], keyed_order
    )

    # A nonsymmetric 3x3 matrix exercises row-major planes across a second
    # size and uses independent sequential products as the reference.
    M3 = Matrix(
        S,
        [
            [x + P - 2, 2*x - 3*P, -3 + P*x],
            [1 + 2*P, x**2 - 4, 5 - P],
            [x - P, 2 + P*x, 3*x + 1],
        ],
    )
    M3_before = M3.__copy__()
    primes3 = [3, 5]
    expected3 = direct_ring_forest(
        M3, [p*p for p in primes3], primes3, 2, kbase=1
    )
    for packed_kappa in (0, 2):
        actual3 = remainder_forest_p2_packed(
            M3, primes3, kbase=1, kappa=packed_kappa
        )
        assert_packed_equals(actual3, expected3, primes3)
    assert M3 == M3_before

    scalar_M = Matrix(S, [[-7 + 3*P + x]])
    scalar_primes = [2, 3, 5]
    scalar_expected = direct_ring_forest(
        scalar_M, [p*p for p in scalar_primes], scalar_primes, 2, kbase=1
    )
    scalar_packed = remainder_forest_p2_packed(
        scalar_M, scalar_primes, kbase=1, kappa=0
    )
    assert_packed_equals(scalar_packed, scalar_expected, scalar_primes)

    empty_packed = remainder_forest_p2_packed(M, [], V=V)
    assert empty_packed.shape == (0, 2, 1, 2)
    assert empty_packed.primes == empty_packed.endpoints == empty_packed.moduli == ()
    assert empty_packed.indices == ()
    assert empty_packed.coefficients.format == "Q"
    assert len(empty_packed.coefficients) == 0

    # The largest 32-bit prime produces coefficients above signed-long range.
    # One transition is enough to verify checked uint64 export at the boundary.
    largest_u32_prime = 4294967291
    largest_u32_square = largest_u32_prime**2
    boundary_M = Matrix(
        S, [[largest_u32_square - 1 + (largest_u32_square - 2)*P]]
    )
    boundary_expected = direct_ring_forest(
        boundary_M,
        [largest_u32_square],
        [largest_u32_prime],
        2,
        kbase=largest_u32_prime - 1,
    )
    boundary_packed = remainder_forest_p2_packed(
        boundary_M,
        [largest_u32_prime],
        kbase=largest_u32_prime - 1,
        kappa=0,
    )
    assert_packed_equals(
        boundary_packed, boundary_expected, [largest_u32_prime]
    )
    assert boundary_packed.matrix(0, 0)[0, 0] == largest_u32_square - 1
    assert boundary_packed.matrix(0, 1)[0, 0] == largest_u32_square - 2
    assert boundary_packed.matrix(0, 0)[0, 0] > 2**63
    assert boundary_packed.matrix(0, 1)[0, 0] > 2**63

    # Memoryviews retain the immutable byte owner after the result object and
    # native arrays are gone, and a later call cannot mutate older results.
    held = remainder_forest_p2_packed(M, [5], kbase=kbase, V=V, kappa=1)
    held_buffer = held.buffer
    held_view = held.matrix(0, 1)
    held_values = tuple(held_view[r, c] for r in range(held.rows)
                        for c in range(held.columns))
    del held
    gc.collect()
    later = remainder_forest_p2_packed(M, [7], kbase=kbase, V=V, kappa=1)
    assert later.buffer is not held_buffer
    assert tuple(held_view[r, c] for r in range(1) for c in range(2)) == held_values

    # Positive control: the legacy path calls Sage result constructors; the
    # packed path calls none of them after its Sage inputs already exist.
    constructor_counts = {"Matrix": 0, "Integers": 0, "PolynomialRing": 0}
    constructors = {
        name: getattr(binding_module, name) for name in constructor_counts
    }

    def observe_constructor(name):
        original = constructors[name]

        def observed(*args, **kwargs):
            constructor_counts[name] += 1
            return original(*args, **kwargs)

        return observed

    for name in constructor_counts:
        setattr(binding_module, name, observe_constructor(name))
    try:
        legacy_probe = remainder_forest_p2(
            M, [25], [5], kbase=kbase, V=V, kappa=1
        )
        assert legacy_probe
        assert constructor_counts["Matrix"] > 0
        assert constructor_counts["Integers"] > 0
        assert constructor_counts["PolynomialRing"] > 0
        for name in constructor_counts:
            constructor_counts[name] = 0
        packed_probe = remainder_forest_p2_packed(
            M, [5], kbase=kbase, V=V, kappa=1
        )
        assert all(count == 0 for count in constructor_counts.values())
    finally:
        for name, constructor in constructors.items():
            setattr(binding_module, name, constructor)
    assert_packed_equals(packed_probe, {0: packed_expected[0]}, [5])

    # A post-native result-construction exception exercises cleanup of GMP and
    # packed allocations; the next call must remain independent and usable.
    original_result_factory = binding_module._construct_packed_p2_result

    def reject_result_construction(*args, **kwargs):
        raise RuntimeError("injected packed result construction failure")

    binding_module._construct_packed_p2_result = reject_result_construction
    try:
        assert_raises(
            RuntimeError,
            lambda: remainder_forest_p2_packed(M, [5], kbase=kbase, V=V),
        )
    finally:
        binding_module._construct_packed_p2_result = original_result_factory
    after_failure = remainder_forest_p2_packed(
        M, [5], kbase=kbase, V=V, kappa=1
    )
    assert_packed_equals(after_failure, {0: packed_expected[0]}, [5])

    for invalid_prime in (1, 9, 2**32 - 1, 2**32):
        assert_raises(
            PackedP2UnsupportedDomain,
            lambda invalid_prime=invalid_prime: remainder_forest_p2_packed(
                M, [invalid_prime], kbase=kbase, V=V
            ),
        )
    assert_raises(
        TypeError,
        lambda: remainder_forest_p2_packed(M, [True], kbase=kbase, V=V),
    )
    assert_raises(
        ValueError,
        lambda: remainder_forest_p2_packed(M, [7, 5], kbase=kbase, V=V),
    )

    assert M == M_before and V == V_before
    print(
        "PASS exact P²/P^n references and packed uint64 planes, ownership, "
        "bypass, and validation"
    )


if __name__ == "__main__":
    main()
