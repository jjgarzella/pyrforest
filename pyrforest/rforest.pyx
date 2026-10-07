# distutils: language=c
# clang c
# Copyright 2023 Edgar Costa, Kiran Kedlaya

r"""
Wrapper for remainder forests.

The remainder forest construction is an efficient mechanism for computing
a sequence of quantities of the form ``V M(0) \cdots M(k_i-1) \pmod{m_i}``,
where `M` is a square matrix over `\ZZ[x]`, `V` is a matrix over `\ZZ`,
and `k_i` and `m_i` are sequences of integers.
"""

from cython cimport sizeof
from cysignals.signals cimport sig_on, sig_off
from libc.stddef cimport size_t
from libc.stdlib cimport malloc, free
from sage.libs.gmp.mpz cimport (
    mpz_clear,
    mpz_init,
    mpz_init_set,
    mpz_init_set_ui,
    mpz_set,
)
from sage.combinat.integer_vector import IntegerVectors
from sage.functions.log import log
from sage.functions.other import ceil
from sage.matrix.constructor import Matrix
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing



cpdef remainder_forest(M, m, k, kbase=0, indices=None, V=None, ans=None, kappa=None,  cutoff=None, projective=False):
    r"""
    Compute modular reductions of matrix products using a remainder forest.

    INPUT:

     - ``M``: a matrix of polynomials with integer coefficients.
     - ``m``: a list or dict of integers, or a function (see below)
     - ``k``: a list or dict of integers, or a function (see below). This list must be strictly monotone.
     - ``kbase``: an integer (defaults to 0).
     - ``indices``: a list or generator arbitrary values (optional). If included,
       we treat ``m`` and ``k`` as lambda functions to be evaluated on these indices.
     - ``V``: a matrix of integers (optional). If omitted, use the identity matrix.
     - ``ans``: a dict of matrices (optional).
     - ``kappa``: a tuning parameter (optional). This controls the number of trees in the forest.
     - `cutoff`: an integer (optional). If specified, answers are truncated to this many columns (counting from the left).
     - ``projective``: a boolean (optional). If True, the answer is allowed to be off by a scalar multiple.

    OUTPUT:

     - if ``ans`` is omitted, a dict ``l`` indexed by ``indices`` (or by default
       ``range(len(m))``) in which
       ``l[i] == V*prod(M.apply_map(lambda x: x(j)) for j in range(kbase, k[i])) % m[i]``.
       If ``ans`` is included, we return ``None`` and update ``ans[i]`` to
       ``ans[i]*(V*prod(M.apply_map(lambda x: x(j)) for j in range(kbase, k[i])) % m[i])``.
       Entries of ``ans`` whose keys do not appear in ``indices`` are unaffected.

    EXAMPLES:

    For every prime `p` up to 1000, we confirm Wilson's theorem that
    `(p-1)! \equiv -1 \pmod{p}`, and identify the primes for which the
    congruence holds modulo `p^2` (these are called Wilson primes)::

        sage: from pyrforest import remainder_forest
        sage: P.<x> = ZZ[]
        sage: M = Matrix([x])
        sage: m = [p^2 for p in prime_range(600)]
        sage: k = [p for p in prime_range(600)]
        sage: l = remainder_forest(M, m, k, kbase=1)
        sage: all(l[i]%k[i] == k[i]-1 for i in range(len(l)))
        True
        sage: [k[i] for i in range(len(l)) if l[i] == m[i]-1]
        [5, 13, 563]

    Redo the previous example using indices and ``return_dict``::

        sage: P.<x> = ZZ[]
        sage: M = Matrix([x])
        sage: indices = prime_range(600)
        sage: m = lambda p: p^2
        sage: k = lambda p: p
        sage: d = remainder_forest(M, m, k, kbase=1, indices=indices)
        sage: all(d[p]%p == p-1 for p in indices)
        True
        sage: [p for p in indices if d[p] == p^2-1]
        [5, 13, 563]

    A simple example where the matrix product computes a sum::

        sage: P.<x> = ZZ[]
        sage: M = Matrix([[1,0],[x,1]])
        sage: m = [n for n in range(1, 100)]
        sage: k = [n for n in range(1, 100)]
        sage: V = Matrix([[0,1]])
        sage: l = remainder_forest(M, m, k, V=V)
        sage: all(l[i][0,1] == 1 for i in range(1, len(l)))
        True
        sage: all(l[i][0,0] == (0 if i%2==0 else (i+1)//2) for i in range(len(l)))
        True
    """
    cdef:
        int rows, deg, dim, kappa1, numcols
        bint mdict, kdict, ansdict, errorflag, proj
        Integer tmp
        long *k1
        long n, i, j, j2, t, kbase1
        mpz_t *A1, *V1, *M1, *m1
        mpz_t z

    # Sanitize input.

    if not M.is_square():
        raise ValueError("Matrix must be square")
    dim = M.dimensions()[0]
    if V is None:
        rows = dim
    else:
        rows = V.dimensions()[0]
        if V.dimensions()[1] != dim:
            raise ValueError("Matrix dimension mismatch")
    if indices is None:
        n = len(m)
        if len(k) != n:
            raise ValueError("m and k must have the same length")
    else:
        try:
            n = len(indices)
        except TypeError:
            n = 0
            for _ in indices:
                n += 1

    # determine the maximum degree
    deg = 0
    for i in range(dim):
        for j in range(dim):
            if M[i,j] in ZZ:
                continue
            deg = max(deg, M[i,j].degree())

    # Translate other inputs into C variables
    kbase1 = kbase
    proj = projective
    ansdict = (ans is not None)
    numcols = dim if cutoff is None else cutoff

    # Allocate and set input variables.
    m1 = <mpz_t *>malloc(n*sizeof(mpz_t))
    k1 = <long *>malloc(n*sizeof(long))
    errorflag = 0
    if indices is None:
        for t in range(n):
            tmp = Integer(m[t])
            mpz_init_set(m1[t], tmp.value)
            k1[t] = k[t]
            if (t == 0 and k1[t] < kbase1) or (t > 0 and k1[t] < k1[t-1]):
                errorflag = 1
    else:
        mdict = isinstance(m, dict)
        kdict = isinstance(k, dict)
        t = 0
        for x in indices:
            tmp = Integer(m[x] if mdict else m(x))
            mpz_init_set(m1[t], tmp.value)
            k1[t] = k[x] if kdict else k(x)
            if (t == 0 and k1[t] < kbase1) or (t > 0 and k1[t] < k1[t-1]):
                errorflag = 1
            t += 1
    if errorflag:
        for i in range(n):
            mpz_clear(m1[i])
        free(m1)
        free(k1)
        raise ValueError("k must be a monotone sequence of values not less than kbase")

    M1 = <mpz_t *>malloc(dim*dim*(deg+1)*sizeof(mpz_t))
    t = 0
    for i in range(dim):
        for j in range(dim):
            for j2 in range(deg+1):
                tmp = Integer(M[i,j][j2])
                mpz_init_set(M1[t], tmp.value)
                t += 1

    V1 = <mpz_t *>malloc(rows*dim*sizeof(mpz_t))
    t = 0
    for i in range(rows):
        for j in range(dim):
            if V is None:
                mpz_init_set_ui(V1[t], 1 if i==j else 0)
            else:
                tmp = Integer(V[i,j])
                mpz_init_set(V1[t], tmp.value)
            t += 1

    if kappa is None:
        kappa1 = 1 if n <= 1 else ceil(log(log(n,2),2)) + 1
    else:
        kappa1 = kappa

    # Allocate output variables.

    A1 = <mpz_t *>malloc(rows*dim*n*sizeof(mpz_t))
    for t in range(rows*dim*n):
        mpz_init(A1[t])
    mpz_init(z)

    try:
        # Call rforest.
        sig_on()
        mproduct(z, m1, n)
        rforest(A1, V1, rows, M1, deg, dim, m1, kbase1, k1, n, z, kappa1)
        sig_off()

        # Retrieve answers. If ans is specified, we assume it is a dict of integer matrices
        # and fill in the entries rather than creating matrices from scratch.

        if indices is None:
            indices = range(n)
        if not ansdict:
            ans = []
        tmp = Integer(0)
        tmp_mat = Matrix(ZZ, rows, numcols)
        t = 0
        for i in range(n):
            for j in range(rows):
                for j1 in range(numcols):
                    tmp.set_from_mpz(A1[t])
                    tmp_mat[j,j1] = tmp
                    t += 1
                t += dim - numcols
            if ansdict:
                ans[indices[i]] *= tmp_mat
            else:
                ans.append(tmp_mat.__copy__())
        if not ansdict:
            return dict(zip(indices, ans))

    # Free malloc-ed memory, even if an exception was raised.

    finally:
        for i in range(dim*dim*(deg+1)):
            mpz_clear(M1[i])
        free(M1)
        for i in range(rows*dim):
            mpz_clear(V1[i])
        free(V1)
        for i in range(n):
            mpz_clear(m1[i])
        free(m1)
        free(k1)
        for i in range(rows*dim*n):
            mpz_clear(A1[i])
        free(A1)
        mpz_clear(z)


def _ring_parent_positions(parent, label, allow_x, require_bivariate=False):
    """Return variable positions, with roles defined by generator order."""
    if parent == ZZ:
        if require_bivariate:
            raise TypeError("%s must use a two-generator integer polynomial ring" % label)
        return {}
    try:
        names = tuple(parent.variable_names())
        coefficient_ring = parent.base_ring()
    except (AttributeError, TypeError):
        raise TypeError("%s must have entries in ZZ or an integer polynomial ring" % label)
    if coefficient_ring != ZZ:
        raise TypeError("%s polynomial coefficients must be exact integers" % label)
    if require_bivariate:
        if len(names) != 2:
            raise TypeError("%s must use a two-generator integer polynomial ring" % label)
        return {"P": 0, "x": 1}
    if len(names) == 1:
        return {"P": 0}
    if len(names) == 2:
        return {"P": 0, "x": 1}
    raise TypeError("%s must use ZZ, ZZ[P], or the same two-generator ring as M" % label)


def _ring_entry_terms(entry, positions):
    """Return an entry's exact terms with exponent tuples in parent order."""
    if not positions:
        return {(): entry}
    terms = entry.dict()
    normalized = {}
    for exponents, coefficient in terms.items():
        try:
            exponents = tuple(exponents)
        except TypeError:
            exponents = (exponents,)
        normalized[exponents] = coefficient
    return normalized


def _ring_integer(value, label):
    if isinstance(value, bool) or value not in ZZ:
        raise TypeError("%s must be an exact integer" % label)
    return Integer(value)


def _ring_batch_value(values, position, index):
    if isinstance(values, dict):
        return values[index]
    if callable(values):
        return values(index)
    return values[position]


def _check_mpz_array_size(count, label):
    size_max = Integer(2) ** (8 * sizeof(size_t)) - 1
    if count < 0 or count > size_max // sizeof(mpz_t):
        raise OverflowError("%s is too large for a native GMP array" % label)


def _ring_matrix_terms(M, nP, label, allow_x, require_bivariate=False):
    """Pack sparse Sage polynomial entries as (P exponent, x degree) maps."""
    positions = _ring_parent_positions(M.base_ring(), label, allow_x,
                                       require_bivariate)
    p_position = positions.get("P")
    x_position = positions.get("x")
    term_maps = []
    degree = 0
    int_max = Integer(2) ** (8 * sizeof(int) - 1) - 1
    for row in range(M.dimensions()[0]):
        for column in range(M.dimensions()[1]):
            coefficients = {}
            for exponents, coefficient in _ring_entry_terms(M[row, column], positions).items():
                p_exponent = exponents[p_position] if p_position is not None else 0
                x_degree = exponents[x_position] if x_position is not None else 0
                if not allow_x and x_degree:
                    raise ValueError("V entries must be independent of x")
                if p_exponent >= nP:
                    continue
                if x_degree > int_max - 1:
                    raise OverflowError("polynomial degree does not fit a C int")
                if coefficient:
                    coefficients[(p_exponent, x_degree)] = coefficient
                if allow_x and x_degree > degree:
                    degree = x_degree
            term_maps.append(coefficients)
    return term_maps, degree


def _ring_matrix_from_terms(terms, rows, dim, nP, modulus, p_name,
                            identity=False):
    coefficient_ring = Integers(modulus)
    polynomial_ring = PolynomialRing(coefficient_ring, p_name)
    entries = []
    for row in range(rows):
        for column in range(dim):
            coefficients = []
            for p in range(nP):
                if identity and p == 0 and row == column:
                    coefficient = 1
                elif terms is None:
                    coefficient = 0
                else:
                    coefficient = terms[row * dim + column].get((p, 0), 0)
                coefficients.append(coefficient_ring(coefficient))
            entries.append(polynomial_ring(coefficients))
    return Matrix(polynomial_ring, rows, dim, entries)


def _remainder_forest_ring(M, m, k, nP, kbase=0, indices=None, V=None,
                           kappa=None, return_state=False, fixed_p2=False,
                           initial_z=None):
    """Shared safe marshalling layer for the native P² and P^n forests."""
    cdef mpz_t *A1 = NULL
    cdef mpz_t *V1 = NULL
    cdef mpz_t *M1 = NULL
    cdef mpz_t *m1 = NULL
    cdef long *k1 = NULL
    cdef mpz_t z
    cdef Integer tmp
    cdef size_t A_initialized = 0
    cdef size_t V_initialized = 0
    cdef size_t M_initialized = 0
    cdef size_t m_initialized = 0
    cdef bint z_initialized = False
    cdef int rows_c, dim_c, deg_c, nP_c, kappa_c
    cdef long n_c, kbase_c
    cdef size_t A_count_c, V_count_c, M_count_c, m_count_c
    cdef size_t t, r, c, p, d

    precision = _ring_integer(nP, "nP")
    if precision < 1:
        raise ValueError("nP must be positive")
    nP = int(precision)
    int_max = Integer(2) ** (8 * sizeof(int) - 1) - 1
    long_min = -(Integer(2) ** (8 * sizeof(long) - 1))
    long_max = -long_min - 1
    if nP > int_max:
        raise OverflowError("nP does not fit a C int")

    if not hasattr(M, "is_square") or not M.is_square():
        raise ValueError("M must be a square Sage matrix")
    dim = M.dimensions()[0]
    if dim <= 0:
        raise ValueError("M must have positive dimension")
    native_dim = dim
    M_terms, deg = _ring_matrix_terms(M, nP, "M", True, True)

    if V is None:
        rows = dim
        V_terms = None
    else:
        if not hasattr(V, "dimensions"):
            raise TypeError("V must be a Sage matrix")
        rows = V.dimensions()[0]
        if rows <= 0 or V.dimensions()[1] != dim:
            raise ValueError("V must have positive rows and the same column count as M")
        V_terms, _ = _ring_matrix_terms(V, nP, "V", False)

    if indices is None:
        try:
            n = len(m)
            k_length = len(k)
        except TypeError:
            raise TypeError("m and k must be sized sequences when indices is omitted")
        if k_length != n:
            raise ValueError("m and k must have the same length")
        index_values = list(range(n))
    else:
        index_values = list(indices)
        n = len(index_values)

    kbase_value = _ring_integer(kbase, "kbase")
    if kbase_value < long_min or kbase_value > long_max:
        raise OverflowError("kbase does not fit a C long")
    kbase_c = int(kbase_value)
    m_values = []
    k_values = []
    modulus_product = Integer(1)
    previous_k = kbase_c
    for position in range(n):
        index = index_values[position]
        m_value = _ring_integer(_ring_batch_value(m, position, index), "each modulus")
        if m_value <= 0:
            raise ValueError("moduli must be positive")
        k_value = _ring_integer(_ring_batch_value(k, position, index), "each endpoint")
        if k_value < kbase_value or (position and k_value < previous_k):
            raise ValueError("k must be a nondecreasing sequence of endpoints not less than kbase")
        if k_value < long_min or k_value > long_max:
            raise OverflowError("endpoint does not fit a C long")
        if k_value - kbase_value > long_max or k_value - previous_k > long_max:
            raise OverflowError("endpoint span does not fit a C long")
        m_values.append(m_value)
        k_values.append(int(k_value))
        modulus_product *= m_value
        previous_k = int(k_value)

    if initial_z is None:
        initial_z_value = modulus_product
    else:
        initial_z_value = _ring_integer(initial_z, "z")
        if initial_z_value <= 0:
            raise ValueError("z must be positive")
        if initial_z_value % modulus_product:
            raise ValueError("z must be divisible by the product of the moduli")

    if kappa is None:
        kappa_value = 1 if n <= 1 else ceil(log(log(n, 2), 2)) + 1
    else:
        kappa_value = _ring_integer(kappa, "kappa")
    if kappa_value < 0 or kappa_value > int_max:
        raise ValueError("kappa must be a nonnegative C int")
    if n > long_max:
        raise OverflowError("batch size does not fit a C long")
    if max(rows, dim, native_dim, deg, nP, int(kappa_value)) > int_max:
        raise OverflowError("matrix dimensions, degree, precision, or kappa do not fit a C int")

    matrix_cells = nP * native_dim * native_dim
    vector_cells = nP * rows * native_dim
    if matrix_cells > int_max or vector_cells > int_max:
        raise OverflowError("ring matrix dimensions exceed the native C int limit")
    M_count = matrix_cells * (deg + 1)
    V_count = vector_cells
    A_count = n * vector_cells
    m_count = n
    for count, label in ((M_count, "M"), (V_count, "V"),
                         (A_count, "outputs"), (m_count, "moduli")):
        _check_mpz_array_size(count, label)
        if count > long_max:
            raise OverflowError("%s array exceeds the native long limit" % label)

    if n == 0:
        if not return_state:
            return {}
        p_name = M.base_ring().variable_names()[0]
        final_V = _ring_matrix_from_terms(V_terms, rows, dim, nP,
                                          initial_z_value, p_name,
                                          identity=(V_terms is None))
        return {}, {"z": initial_z_value, "final_V": final_V}

    rows_c = rows
    dim_c = native_dim
    deg_c = deg
    nP_c = nP
    n_c = n
    kappa_c = int(kappa_value)
    A_count_c = A_count
    V_count_c = V_count
    M_count_c = M_count
    m_count_c = m_count

    try:
        M1 = <mpz_t *>malloc(M_count_c * sizeof(mpz_t))
        V1 = <mpz_t *>malloc(V_count_c * sizeof(mpz_t))
        m1 = <mpz_t *>malloc(m_count_c * sizeof(mpz_t))
        k1 = <long *>malloc(m_count_c * sizeof(long))
        A1 = <mpz_t *>malloc(A_count_c * sizeof(mpz_t))
        if M1 == NULL or V1 == NULL or m1 == NULL or k1 == NULL or A1 == NULL:
            raise MemoryError("unable to allocate native ring forest buffers")

        for r in range(native_dim):
            for c in range(native_dim):
                if r < dim and c < dim:
                    entry_terms = M_terms[r * dim + c]
                else:
                    entry_terms = {}
                for p in range(nP):
                    for d in range(deg + 1):
                        if r == c and r >= dim:
                            coefficient = 1 if p == 0 and d == 0 else 0
                        else:
                            coefficient = entry_terms.get((p, d), 0)
                        tmp = Integer(coefficient)
                        mpz_init_set(M1[M_initialized], tmp.value)
                        M_initialized += 1

        for p in range(nP):
            for r in range(rows):
                for c in range(native_dim):
                    if V_terms is None:
                        coefficient = 1 if p == 0 and r == c else 0
                    elif c >= dim:
                        coefficient = 0
                    else:
                        coefficient = V_terms[r * dim + c].get((p, 0), 0)
                    tmp = Integer(coefficient)
                    mpz_init_set(V1[V_initialized], tmp.value)
                    V_initialized += 1

        for t in range(n):
            tmp = m_values[t]
            mpz_init_set(m1[t], tmp.value)
            m_initialized += 1
            k1[t] = k_values[t]

        for t in range(A_count_c):
            mpz_init(A1[t])
            A_initialized += 1
        mpz_init(z)
        z_initialized = True

        sig_on()
        try:
            if initial_z is None:
                tmp = modulus_product
                mpz_set(z, tmp.value)
            else:
                tmp = initial_z_value
                mpz_set(z, tmp.value)
            if fixed_p2:
                rforest_p2(A1, V1, rows_c, M1, deg_c, dim_c, m1,
                           kbase_c, k1, n_c, z, kappa_c)
            else:
                rforest_pn(A1, V1, rows_c, M1, deg_c, dim_c, nP_c, m1,
                           kbase_c, k1, n_c, z, kappa_c)
        finally:
            sig_off()

        results = {}
        p_name = M.base_ring().variable_names()[0]
        for t in range(n):
            output_ring = PolynomialRing(Integers(m_values[t]), p_name)
            output_entries = []
            for r in range(rows):
                for c in range(dim):
                    coefficients = []
                    for p in range(nP):
                        output_offset = ((t * nP + p) * rows + r) * native_dim + c
                        tmp.set_from_mpz(A1[output_offset])
                        coefficients.append(output_ring.base_ring()(tmp))
                    output_entries.append(output_ring(coefficients))
            results[index_values[t]] = Matrix(output_ring, rows, dim,
                                              output_entries)

        if not return_state:
            return results

        z_value = Integer(0)
        z_value.set_from_mpz(z)
        state_ring = PolynomialRing(Integers(z_value), p_name)
        state_entries = []
        for r in range(rows):
            for c in range(dim):
                coefficients = []
                for p in range(nP):
                    state_offset = (p * rows + r) * native_dim + c
                    tmp.set_from_mpz(V1[state_offset])
                    coefficients.append(state_ring.base_ring()(tmp))
                state_entries.append(state_ring(coefficients))
        final_V = Matrix(state_ring, rows, dim, state_entries)
        return results, {"z": z_value, "final_V": final_V}

    finally:
        if z_initialized:
            mpz_clear(z)
        for t in range(M_initialized):
            mpz_clear(M1[t])
        if M1 != NULL:
            free(M1)
        for t in range(V_initialized):
            mpz_clear(V1[t])
        if V1 != NULL:
            free(V1)
        for t in range(m_initialized):
            mpz_clear(m1[t])
        if m1 != NULL:
            free(m1)
        if k1 != NULL:
            free(k1)
        for t in range(A_initialized):
            mpz_clear(A1[t])
        if A1 != NULL:
            free(A1)


def remainder_forest_p2(M, m, k, kbase=0, indices=None, V=None,
                        kappa=None, return_state=False, *, z=None):
    r"""Compute a remainder forest over the truncated formal ring ``ZZ[P]/(P^2)``.

    ``M`` is a square Sage matrix over an exact two-generator integer
    polynomial ring. Its first generator is formal ``P`` and its second is
    transition variable ``x``, regardless of their Sage names. ``x`` is
    evaluated at each transition index; ``P`` remains formal and is truncated
    modulo ``P^2``. ``V`` may be a rectangular matrix over ``ZZ``, a
    one-generator integer polynomial ring (whose generator is formal ``P``),
    or the same two-generator ring, and must be independent of ``x``. Each output is
    a Sage matrix over ``(ZZ/m[i]ZZ)[P]`` with only the coefficients of 1 and
    P retained. The usual ``m``, ``k``, ``indices``, ``kbase`` and ``kappa``
    conventions of :func:`remainder_forest` apply, including exclusive
    endpoints and right multiplication.

    If ``return_state`` is true, return ``(results, state)`` where the state
    dict has keys ``'z'`` (the exact residual modulus) and ``'final_V'`` (the
    final ``V`` modulo that modulus).
    The keyword-only ``z`` may provide an initial positive multiple of all
    endpoint moduli; by default it is their product, so the residual modulus
    is one and ``final_V`` is zero in the modulus-one ring.
    """
    return _remainder_forest_ring(M, m, k, 2, kbase, indices, V, kappa,
                                  return_state, True, z)


def remainder_forest_pn(M, m, k, nP, kbase=0, indices=None, V=None,
                        kappa=None, return_state=False, *, z=None):
    r"""Compute a remainder forest over the truncated formal ring ``ZZ[P]/(P^nP)``.

    ``M`` and ``V`` use the exact coefficient and transition-variable contract
    documented by :func:`remainder_forest_p2`. ``nP`` is an integer at least
    one. Each output is a Sage matrix over ``(ZZ/m[i]ZZ)[P]`` truncated to
    exponents below ``nP``. With ``return_state=True``, the return value is
    ``(results, state)`` where the state dict has keys ``'z'`` and ``'final_V'``.
    The keyword-only ``z`` has the same initial-modulus meaning as in
    :func:`remainder_forest_p2`.
    """
    return _remainder_forest_ring(M, m, k, nP, kbase, indices, V, kappa,
                                  return_state, False, z)

def inflate_matrix(M, vars, e):
    """
    Compute an inflation of a matrix over a polynomial ring.

    The return value is a matrix over another polynomial ring with the specified variables omitted.
    It is a block matrix with blocks corresponding to monomials in the omitted variables.

    EXAMPLES::

        sage: from pyrforest.rforest import inflate_matrix
        sage: R.<x,y> = PolynomialRing(ZZ, 2)
        sage: M = Matrix([[x+y, x*y], [x-y, 1]])
        sage: inflate_matrix(M, [y], 2)
        [ x  1  0  x]
        [ 0  x  0  0]
        [ x -1  1  0]
        [ 0  x  0  1]
        sage: inflate_matrix(M, [x,y], 2)
        [ 0  0  1  0  0  0]
        [ 0  0  1  0  0  0]
        [ 0  0  0  0  0  0]
        [ 0  0  1  1  0  0]
        [ 0  0 -1  0  1  0]
        [ 0  0  0  0  0  1]
    """
    R = M.base_ring()
    gens = R.gens()
    if any(x not in gens for x in vars):
        raise ValueError("Invalid key in d")

    dim = M.dimensions()
    n = R.ngens()
    n1 = len(vars)
    vec = [tuple(v[:-1]) for v in IntegerVectors(n=e-1, k=n1+1)]
    nv = len(vec)

    vec_dict = {}
    for v1 in vec:
        for v2 in vec:
            v = tuple(v1[t] - v2[t] for t in range(n1))
            if v in vec:
                vec_dict[v1, v2] = v

    K = R.base_ring()
    vars2 = tuple(x for x in gens if x not in vars)
    if n == n1:
        R1 = K
    else:
        R1 = PolynomialRing(K, n-n1, names=vars2)
    R2 = PolynomialRing(R1, n1, names=vars)
    # Must explicitly define the coercion map from R to R2.
    var_target = [(R2.gens()[vars.index(x)] if x in vars else R2(R1.gens()[vars2.index(x)])) for x in gens]
    h = R.hom(var_target)
    mat_dict = {}
    for t1 in range(dim[0]):
        for t2 in range(dim[1]):
            tmp = h(M[t1, t2])
            for j1 in range(nv):
                v1 = vec[j1]
                for j2 in range(nv):
                    v2 = vec[j2]
                    if (v1, v2) in vec_dict:
                        mat_dict[t1*nv + j1, t2*nv + j2] = tmp[vec_dict[v1, v2]]
    return Matrix(mat_dict)

def deflate_matrix(M, vars, e, R1=None):
    """
    Undo the effect of ``inflate_matrix``.

    EXAMPLES::

        sage: from pyrforest.rforest import inflate_matrix, deflate_matrix
        sage: R.<x,y> = PolynomialRing(ZZ, 2)
        sage: M = Matrix([[x+y, x*y], [x-y, 1]])
        sage: deflate_matrix(inflate_matrix(M, [y], 2), [y], 2) == M
        True
        sage: deflate_matrix(inflate_matrix(M, [x,y], 2), [x,y], 2)
        [x + y     0]
        [x - y     1]
    """
    K = M.base_ring()
    dim = M.dimensions()
    n1 = len(vars)
    vec = [tuple(v[:-1]) for v in IntegerVectors(n=e-1, k=n1+1)]
    nv = len(vec)
    mat_dict = {}
    j2 = vec.index((0,)*n1)
    if not R1:
        R1 = PolynomialRing(K, n1, names=vars)
    for t1 in range(dim[0] // nv):
        for t2 in range(dim[1] // nv):
            d = {}
            for j1 in range(nv):
                d[vec[j1]] = M[t1*nv + j1, t2*nv + j2]
            mat_dict[t1, t2] = R1(d)
    return Matrix(mat_dict)

def remainder_forest_generic_prime(M, d, e, k, indices=None, m=None, kbase=0, V=None, ans=None, kappa=None):
    """
    Compute a remainder forest using one or more "generic primes".

    INPUT::

     - ``M``: a matrix of multivariate polynomials with integer coefficients.
     - ``d``: a dictionary keyed by variables in the base ring of ``M``. Exactly one
      variable must be omitted.
     - ``e``: a positive integer.
     - ``k``: a list or dict of integers, or a function (see below). This list must be monotone;
         if a dict or function is given, it must evaluate to a monotone list when applied to ``indices``.
     - ``indices``: a list or generator arbitrary values.
     - ``kbase``: an integer (defaults to 0).
     - ``m``: a lambda function on indices (see below).
     - ``V``: a matrix of polynomials. If omitted, use the identity matrix.
     - ``ans``: a dict of matrices (optional).
     - ``kappa``: a tuning parameter (optional). This controls the number of trees in the forest.

    OUTPUT::

     As for ``remainder_forest`` except that ``M`` is evaluated with the variables named in ``d``
     truncated modulo the ``e``-th power of the ideal they jointly generate, then specialized as per ``d``;
     the other variable is evaluated at successive integers as before. Also, the moduli may not be specified freely:
     if ``m`` is omitted, we impose that ``m[p] == p^e``.

     While ``d`` must include all but one variable, the values of ``d`` may omit keys. In these cases, the corresponding
     variable is left unevaluated.

    NOTE::

     In practice, it will usually be more efficient to call `remainder_forest` directly.

    EXAMPLES::

        sage: from pyrforest import remainder_forest_generic_prime
        sage: P.<x,y,z> = ZZ[]
        sage: M = Matrix([[x+y+1,z],[y+z,1]])
        sage: indices = prime_range(2, 50)
        sage: k = {p: p-1 for p in indices}
        sage: V = Matrix([[0,1+y]])
        sage: d = {x: {p: p for p in indices}, y: {p: -p for p in indices}}
        sage: ans, _ = remainder_forest_generic_prime(M, d, 2, k, indices=indices, V=V)
        sage: ans2 = {p: V.apply_map(lambda t: t.subs({x: p, y: -p})) for p in indices}
        sage: ans2 = {p: ans2[p] * prod(M.apply_map(lambda t: t.subs({x: p, y: -p, z: j})) for j in range(k[p])) for p in indices}
        sage: all((ans[p] - ans2[p]) % p^2 == 0 for p in indices)
        True

    """
    if not M.is_square():
        raise ValueError("Matrix must be square")
    dim = M.dimensions()[0]
    if V is None:
        rows = dim
    else:
        rows = V.dimensions()[0]
        if V.dimensions()[1] != dim:
            raise ValueError("Matrix dimension mismatch")

    R = M.base_ring()
    if any(x not in R.gens() for x in d):
        raise ValueError("Invalid key in d")
    x = [x for x in R.gens() if x not in d]
    if len(x) != 1:
       raise ValueError("d must omit exactly one variable of M")

    vars = list(d.keys())
    if V is None:
        V1 = None
    else:
        V1 = inflate_matrix(V.change_ring(R), vars, e)
    M1 = inflate_matrix(M, vars, e)
    if m is None:
        m = {p: p**e for p in indices}
    else:
        m = {p: m(p) for p in indices}
    tmp = remainder_forest(M1, m, k, kbase, indices, V1, kappa)
    if indices is None:
        indices = tmp.keys()
    ansdict = (ans is not None)
    if not ansdict:
        ans = {i: 1 for i in indices}

    R1 = PolynomialRing(R.base_ring(), len(vars), names=vars)
    for i in indices:
        M2 = deflate_matrix(tmp[i], vars, e, R1)
        d2 = {x: d[x][i] for x in R1.gens() if i in d[x]}
        ans[i] *= M2.apply_map(lambda t,d2=d2: t.subs(d2))

    if ansdict:
        return None
    return ans, R1
