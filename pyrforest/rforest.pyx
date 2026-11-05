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
from cpython.pycapsule cimport (
    PyCapsule_GetContext,
    PyCapsule_IsValid,
    PyCapsule_New,
    PyCapsule_SetContext,
)
from cpython.object cimport PyObject
from cpython.ref cimport Py_DECREF, Py_INCREF
from cysignals.signals cimport sig_on, sig_off
from cpython.bytes cimport PyBytes_FromStringAndSize
from libc.stddef cimport size_t
from libc.stdint cimport uint64_t
from libc.stdlib cimport malloc, free
from sage.libs.gmp.mpz cimport (
    mpz_clear,
    mpz_cmp,
    mpz_export,
    mpz_init,
    mpz_init_set,
    mpz_init_set_ui,
    mpz_mod,
    mpz_sgn,
    mpz_set,
)
from sage.libs.gmp.types cimport mpz_srcptr, mpz_t
from sage.combinat.integer_vector import IntegerVectors
from sage.functions.log import log
from sage.functions.other import ceil
from sage.matrix.constructor import Matrix
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.polynomial.laurent_polynomial_ring_base import LaurentPolynomialRing_generic
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from pyrforest.native_pn cimport (
    PYRFOREST_PN_GMP_RESULT_ABI_VERSION,
    PYRFOREST_PN_GMP_RESULT_CAPSULE_NAME,
    pyrforest_pn_gmp_result_v1,
    pyrforest_pn_gmp_result_v1_init,
)



PACKED_P2_UINT64_CAPABILITY = 1
PN_GMP_RESULT_ABI_VERSION = PYRFOREST_PN_GMP_RESULT_ABI_VERSION
PN_GMP_RESULT_CAPABILITY = "pn-gmp-result-v1"
PN_GMP_RESULT_CAPSULE_NAME = "pyrforest.PnGmpForestResult/v1"


class PackedP2UnsupportedDomain(ValueError):
    """The packed P² API cannot represent the requested prime domain."""


cdef object _PACKED_P2_RESULT_TOKEN = object()


cdef class PackedP2ForestResult:
    """Immutable owner of a validated row-major P⁰/P¹ uint64 result.

    Metadata sequences are copied to tuples. Direct construction validates
    their alignment, prime domain, canonical words, dimensions, and buffer size.
    """

    cdef object _buffer
    cdef readonly object primes
    cdef readonly object indices
    cdef readonly object endpoints
    cdef readonly object moduli
    cdef readonly int rows
    cdef readonly int columns

    def __cinit__(self, object buffer, object primes, object indices,
                  object endpoints, object moduli, object rows, object columns,
                  object construction_token=None):
        cdef const uint64_t[::1] coefficient_words
        cdef uint64_t modulus_word
        cdef Py_ssize_t plane_words_c
        cdef Py_ssize_t first_word
        cdef Py_ssize_t last_word
        cdef Py_ssize_t word_index
        cdef Py_ssize_t batch_position

        if not isinstance(buffer, bytes):
            raise TypeError("packed result storage must be immutable bytes")

        try:
            prime_values = tuple(primes)
            index_values = tuple(indices)
            endpoint_values = tuple(endpoints)
            modulus_values = tuple(moduli)
        except TypeError:
            raise TypeError("packed result metadata must be iterable sequences")

        batch_count = len(prime_values)
        if (len(index_values) != batch_count
                or len(endpoint_values) != batch_count
                or len(modulus_values) != batch_count):
            raise ValueError("packed result metadata sequences must have equal lengths")

        rows_value = _ring_integer(rows, "rows")
        columns_value = _ring_integer(columns, "columns")
        plane_count, packed_count, expected_nbytes = _packed_p2_output_size(
            batch_count, rows_value, columns_value, sizeof(int),
            sizeof(long), sizeof(size_t), sizeof(Py_ssize_t)
        )
        expected_nbytes = int(expected_nbytes)
        if len(buffer) != expected_nbytes:
            raise ValueError("packed result storage does not match its shape")

        normalized_primes = []
        normalized_endpoints = []
        normalized_moduli = []
        for position in range(batch_count):
            prime_value = _ring_integer(prime_values[position], "each prime")
            if not _packed_p2_endpoint_supported(prime_value, sizeof(long)):
                raise PackedP2UnsupportedDomain(
                    "packed result prime endpoint is outside the native C long/uint32 domain"
                )
            endpoint_value = _ring_integer(
                endpoint_values[position], "each endpoint"
            )
            modulus_value = _ring_integer(
                modulus_values[position], "each modulus"
            )
            if endpoint_value != prime_value:
                raise ValueError("packed result endpoints must equal their primes")
            if modulus_value != prime_value * prime_value:
                raise ValueError("packed result moduli must be the prime squares")
            if (construction_token is not _PACKED_P2_RESULT_TOKEN
                    and not prime_value.is_prime(proof=True)):
                raise PackedP2UnsupportedDomain(
                    "packed result metadata requires prime values"
                )
            normalized_primes.append(prime_value)
            normalized_endpoints.append(endpoint_value)
            normalized_moduli.append(modulus_value)

        if construction_token is not _PACKED_P2_RESULT_TOKEN:
            coefficient_words = memoryview(buffer).cast("Q")
            plane_words_c = <Py_ssize_t>plane_count
            for batch_position in range(<Py_ssize_t>batch_count):
                modulus_word = <uint64_t>int(
                    normalized_moduli[batch_position]
                )
                first_word = batch_position * plane_words_c
                last_word = first_word + plane_words_c
                for word_index in range(first_word, last_word):
                    if coefficient_words[word_index] >= modulus_word:
                        raise ValueError(
                            "packed result coefficient is not canonical modulo p²"
                        )

        self._buffer = buffer
        self.primes = tuple(normalized_primes)
        self.indices = index_values
        self.endpoints = tuple(normalized_endpoints)
        self.moduli = tuple(normalized_moduli)
        self.rows = int(rows_value)
        self.columns = int(columns_value)

    @property
    def buffer(self):
        """Immutable native-endian bytes backing the machine-word view."""
        return self._buffer

    @property
    def coefficients(self):
        """Read-only flat memoryview with format ``Q`` (uint64)."""
        return memoryview(self._buffer).cast("Q")

    @property
    def shape(self):
        """Logical shape ``(batch, 2, rows, columns)``."""
        return (len(self.primes), 2, self.rows, self.columns)

    @property
    def plane_names(self):
        return ("P^0", "P^1")

    def matrix(self, batch_position, plane):
        """Return one read-only row-major ``(rows, columns)`` plane view."""
        batch_position = int(batch_position)
        plane = int(plane)
        if batch_position < 0 or batch_position >= len(self.primes):
            raise IndexError("packed result batch position is out of range")
        if plane not in (0, 1):
            raise IndexError("packed P² plane must be 0 or 1")
        cells = self.rows * self.columns
        first = (batch_position * 2 + plane) * cells * 8
        last = first + cells * 8
        return memoryview(self._buffer)[first:last].cast(
            "Q", shape=(self.rows, self.columns)
        )

    def __repr__(self):
        return "PackedP2ForestResult(primes=%d, shape=%r, readonly=True)" % (
            len(self.primes), self.shape
        )


def _construct_packed_p2_result(*args):
    """Private construction seam used to exercise post-native cleanup."""
    return PackedP2ForestResult(*args, _PACKED_P2_RESULT_TOKEN)


cdef object _NATIVE_PN_RESULT_TOKEN = object()


cdef void _native_pn_capsule_destructor(object capsule) noexcept:
    cdef void *owner_pointer

    if not PyCapsule_IsValid(capsule, PYRFOREST_PN_GMP_RESULT_CAPSULE_NAME):
        return
    owner_pointer = PyCapsule_GetContext(capsule)
    if owner_pointer != NULL:
        Py_DECREF(<object><PyObject *>owner_pointer)


cdef class NativePnForestResult:
    """Immutable owner for an arbitrary-precision PN GMP forest result.

    Instances can only be created by the checked native producer. GMP storage
    is not exposed as Python integers or a byte buffer. ``as_capsule()`` makes
    a fresh, read-only ABI view whose lifetime pins this owner.
    """

    cdef mpz_t *_coefficients
    cdef mpz_t *_moduli
    cdef long *_endpoints
    cdef size_t _coefficient_count
    cdef size_t _batch_count
    cdef bint _ready
    cdef pyrforest_pn_gmp_result_v1 _payload
    cdef readonly object indices
    cdef readonly object endpoints
    cdef readonly object moduli
    cdef readonly int nP
    cdef readonly int rows
    cdef readonly int columns

    def __cinit__(self, object construction_token=None):
        self._coefficients = NULL
        self._moduli = NULL
        self._endpoints = NULL
        self._coefficient_count = 0
        self._batch_count = 0
        self._ready = False
        self.indices = ()
        self.endpoints = ()
        self.moduli = ()
        self.nP = 0
        self.rows = 0
        self.columns = 0
        pyrforest_pn_gmp_result_v1_init(
            &self._payload, 0, 0, 0, 0, NULL, NULL, NULL
        )
        if construction_token is not _NATIVE_PN_RESULT_TOKEN:
            raise TypeError(
                "NativePnForestResult instances are created by the PN forest"
            )

    cdef void _adopt(self, mpz_t *coefficients, mpz_t *moduli,
                     long *endpoints, size_t coefficient_count,
                     size_t batch_count, int nP, int rows, int columns):
        self._coefficients = coefficients
        self._moduli = moduli
        self._endpoints = endpoints
        self._coefficient_count = coefficient_count
        self._batch_count = batch_count
        self.nP = nP
        self.rows = rows
        self.columns = columns
        pyrforest_pn_gmp_result_v1_init(
            &self._payload,
            batch_count,
            <size_t>nP,
            <size_t>rows,
            <size_t>columns,
            <const long *>endpoints,
            <mpz_srcptr>moduli,
            <mpz_srcptr>coefficients,
        )
        self._ready = True

    def __dealloc__(self):
        cdef size_t i
        for i in range(self._coefficient_count):
            mpz_clear(self._coefficients[i])
        if self._coefficients != NULL:
            free(self._coefficients)
        for i in range(self._batch_count):
            mpz_clear(self._moduli[i])
        if self._moduli != NULL:
            free(self._moduli)
        if self._endpoints != NULL:
            free(self._endpoints)

    @property
    def shape(self):
        """Logical coefficient layout ``(batch, nP, rows, columns)``."""
        return (self._batch_count, self.nP, self.rows, self.columns)

    def as_capsule(self):
        """Return a fresh v1 capsule that holds this owner alive."""
        cdef object capsule
        if not self._ready:
            raise RuntimeError("native PN result storage has not been adopted")
        capsule = PyCapsule_New(
            <void *>&self._payload,
            PYRFOREST_PN_GMP_RESULT_CAPSULE_NAME,
            _native_pn_capsule_destructor,
        )
        Py_INCREF(self)
        try:
            PyCapsule_SetContext(capsule, <void *>self)
        except:
            Py_DECREF(self)
            raise
        return capsule

    def __repr__(self):
        return "NativePnForestResult(shape=%r, endpoints=%r)" % (
            self.shape, self.endpoints
        )


cdef NativePnForestResult _new_native_pn_result(
        size_t batch_count, int nP, int rows, int columns,
        object indices, object endpoint_values, object modulus_values):
    cdef NativePnForestResult result
    cdef object index_tuple = tuple(indices)
    cdef object endpoint_tuple = tuple(endpoint_values)
    cdef object modulus_tuple = tuple(modulus_values)

    if (len(index_tuple) != batch_count
            or len(endpoint_tuple) != batch_count
            or len(modulus_tuple) != batch_count):
        raise ValueError("native PN result metadata lengths are inconsistent")

    result = NativePnForestResult(_NATIVE_PN_RESULT_TOKEN)
    result.indices = index_tuple
    result.endpoints = endpoint_tuple
    result.moduli = modulus_tuple
    return result


def _construct_native_pn_result(batch_count, nP, rows, columns, indices,
                                endpoint_values, modulus_values):
    """Private metadata construction seam used by cleanup tests."""
    return _new_native_pn_result(
        batch_count, nP, rows, columns, indices, endpoint_values, modulus_values
    )


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
    if isinstance(parent, LaurentPolynomialRing_generic):
        raise TypeError("%s must use an ordinary polynomial ring, not a Laurent polynomial ring" % label)
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


def _packed_p2_endpoint_limit(long_size):
    """Return the packed prime cap for a signed C long of ``long_size`` bytes."""
    uint32_max = Integer(2) ** 32 - 1
    long_max = Integer(2) ** (8 * long_size - 1) - 1
    return min(uint32_max, long_max)


def _packed_p2_endpoint_supported(prime, long_size):
    """Whether ``prime`` can be used as a packed endpoint on this C ABI."""
    return 2 <= prime <= _packed_p2_endpoint_limit(long_size)


def _packed_p2_output_size(batch_count, rows, columns, int_size, long_size,
                           size_t_size, py_ssize_t_size):
    """Validate packed P² dimensions and return plane words, total words, bytes.

    The plane and total counts use the same C-int, native-long, size_t, and
    Python-buffer limits for direct result construction and factory preflight.
    Width arguments are explicit so allocation-free boundary tests can model
    other ABIs without claiming to run them.
    """
    batch_count = Integer(batch_count)
    rows = Integer(rows)
    columns = Integer(columns)
    if batch_count < 0:
        raise ValueError("packed result batch count cannot be negative")
    if rows <= 0 or columns <= 0:
        raise ValueError("packed result rows and columns must be positive")
    if min(int_size, long_size, size_t_size, py_ssize_t_size) <= 0:
        raise ValueError("packed result ABI widths must be positive")

    int_max = Integer(2) ** (8 * int_size - 1) - 1
    if rows > int_max or columns > int_max:
        raise OverflowError("packed result dimensions do not fit a C int")
    plane_count = 2 * rows * columns
    if plane_count > int_max:
        raise OverflowError("packed result plane exceeds the native C int limit")

    packed_count = batch_count * plane_count
    long_max = Integer(2) ** (8 * long_size - 1) - 1
    if packed_count > long_max:
        raise OverflowError("packed result exceeds the native C long limit")

    size_t_max = Integer(2) ** (8 * size_t_size) - 1
    if packed_count > size_t_max // sizeof(uint64_t):
        raise OverflowError("packed result is too large for native size_t")
    py_ssize_max = Integer(2) ** (8 * py_ssize_t_size - 1) - 1
    if packed_count > py_ssize_max // sizeof(uint64_t):
        raise OverflowError("packed result is too large for a Python buffer")
    return plane_count, packed_count, packed_count * sizeof(uint64_t)


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
                           initial_z=None, packed_p2=False, native_pn=False):
    """Shared safe marshalling layer for the native P² and P^n forests."""
    cdef mpz_t *A1 = NULL
    cdef mpz_t *V1 = NULL
    cdef mpz_t *M1 = NULL
    cdef mpz_t *m1 = NULL
    cdef long *k1 = NULL
    cdef uint64_t *packed_words = NULL
    cdef mpz_t z
    cdef Integer tmp
    cdef Integer prime_value
    cdef NativePnForestResult native_result
    cdef size_t A_initialized = 0
    cdef size_t V_initialized = 0
    cdef size_t M_initialized = 0
    cdef size_t m_initialized = 0
    cdef size_t packed_count_c = 0
    cdef size_t exported_words = 0
    cdef bint z_initialized = False
    cdef int rows_c, dim_c, deg_c, nP_c, kappa_c
    cdef long n_c, kbase_c
    cdef size_t A_count_c, V_count_c, M_count_c, m_count_c
    cdef Py_ssize_t packed_nbytes_c = 0
    cdef size_t t, r, c, p, d, packed_offset

    precision = _ring_integer(nP, "nP")
    if precision < 1:
        raise ValueError("nP must be positive")
    nP = int(precision)
    int_max = Integer(2) ** (8 * sizeof(int) - 1) - 1
    long_min = -(Integer(2) ** (8 * sizeof(long) - 1))
    long_max = -long_min - 1
    if nP > int_max:
        raise OverflowError("nP does not fit a C int")
    if packed_p2 and (not fixed_p2 or nP != 2):
        raise ValueError("packed output is available only for the fixed P² forest")
    if native_pn and (fixed_p2 or packed_p2):
        raise ValueError("native GMP output is available only for generic P^n")
    if native_pn and return_state:
        raise ValueError("native GMP output does not expose final_V state")
    packed_endpoint_limit = (
        _packed_p2_endpoint_limit(sizeof(long)) if packed_p2 else None
    )

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
            if not packed_p2:
                k_length = len(k)
        except TypeError:
            if packed_p2:
                raise TypeError("primes must be a sized sequence when indices is omitted")
            raise TypeError("m and k must be sized sequences when indices is omitted")
        if not packed_p2 and k_length != n:
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
    prime_values = []
    modulus_product = Integer(1)
    previous_k = kbase_c
    for position in range(n):
        index = index_values[position]
        if packed_p2:
            prime_value = _ring_integer(
                _ring_batch_value(m, position, index), "each prime"
            )
            if (not _packed_p2_endpoint_supported(prime_value, sizeof(long))
                    or not prime_value.is_prime(proof=True)):
                raise PackedP2UnsupportedDomain(
                    "packed P² output requires primes in the native C long/uint32 endpoint domain (maximum %s)"
                    % packed_endpoint_limit
                )
            prime_values.append(prime_value)
            m_value = prime_value * prime_value
            k_value = prime_value
        else:
            m_value = _ring_integer(
                _ring_batch_value(m, position, index), "each modulus"
            )
            if m_value <= 0:
                raise ValueError("moduli must be positive")
            k_value = _ring_integer(
                _ring_batch_value(k, position, index), "each endpoint"
            )
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

    if packed_p2:
        _, packed_count, packed_nbytes = (
            _packed_p2_output_size(
                n, rows, dim, sizeof(int), sizeof(long), sizeof(size_t),
                sizeof(Py_ssize_t)
            )
        )
        packed_count_c = <size_t>packed_count
        packed_nbytes_c = <Py_ssize_t>packed_nbytes

    if n == 0:
        if packed_p2:
            return PackedP2ForestResult(
                b"", tuple(), tuple(index_values), tuple(), tuple(), rows, dim
            )
        if native_pn:
            native_result_factory = globals()["_construct_native_pn_result"]
            native_result = native_result_factory(
                0, nP, rows, dim, index_values, (), ()
            )
            native_result._adopt(NULL, NULL, NULL, 0, 0,
                                 nP, rows, dim)
            return native_result
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
        if packed_p2 and packed_count_c:
            packed_words = <uint64_t *>malloc(
                packed_count_c * sizeof(uint64_t)
            )
            if packed_words == NULL:
                raise MemoryError("unable to allocate packed uint64 output buffer")

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

        if packed_p2:
            packed_offset = 0
            for t in range(n):
                for p in range(2):
                    for r in range(rows):
                        for c in range(dim):
                            output_offset = ((t * 2 + p) * rows + r) * native_dim + c
                            if (mpz_sgn(A1[output_offset]) < 0
                                    or mpz_cmp(A1[output_offset], m1[t]) >= 0):
                                raise ArithmeticError(
                                    "native P² output is not canonical modulo p²"
                                )
                            packed_words[packed_offset] = 0
                            exported_words = 0
                            mpz_export(
                                &packed_words[packed_offset], &exported_words,
                                -1, sizeof(uint64_t), 0, 0, A1[output_offset]
                            )
                            if exported_words > 1:
                                raise OverflowError(
                                    "canonical P² coefficient does not fit uint64"
                                )
                            packed_offset += 1
            packed_raw = PyBytes_FromStringAndSize(
                <char *>packed_words, packed_nbytes_c
            )
            packed_result_factory = globals()["_construct_packed_p2_result"]
            return packed_result_factory(
                packed_raw,
                tuple(int(value) for value in prime_values),
                tuple(index_values),
                tuple(k_values),
                tuple(int(value) for value in m_values),
                rows,
                dim,
            )

        if native_pn:
            # The native algorithm returns modular residues. Normalize in GMP
            # before publishing so all branches preserve Sage's canonical
            # residue representatives, including any future signed native
            # intermediate path.
            for t in range(n):
                for p in range(nP):
                    for r in range(rows):
                        for c in range(dim):
                            output_offset = (
                                ((t * nP + p) * rows + r) * native_dim + c
                            )
                            mpz_mod(A1[output_offset], A1[output_offset], m1[t])
                            if (mpz_sgn(A1[output_offset]) < 0
                                    or mpz_cmp(A1[output_offset], m1[t]) >= 0):
                                raise ArithmeticError(
                                    "native PN output is not canonical modulo its modulus"
                                )

            native_result_factory = globals()["_construct_native_pn_result"]
            native_result = native_result_factory(
                n, nP, rows, dim, index_values, k_values, m_values
            )
            native_result._adopt(A1, m1, k1, A_count_c, m_count_c,
                                 nP_c, rows_c, dim_c)
            # The owner is now solely responsible for clearing these arrays.
            A1 = NULL
            A_initialized = 0
            m1 = NULL
            m_initialized = 0
            k1 = NULL
            return native_result

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
        if packed_words != NULL:
            free(packed_words)


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


def remainder_forest_p2_packed(M, primes, kbase=0, indices=None, V=None,
                               kappa=None):
    r"""Return a packed P² forest for primes supported by the native ABI.

    For each prime ``p``, compute ``V*M(kbase)*...*M(p-1)`` modulo ``p^2``
    with formal ``P^2 = 0``. The input matrix and optional rectangular ``V``
    follow :func:`remainder_forest_p2`; when ``V`` is omitted the identity is
    used. ``primes`` may be a sized sequence, dict, or callable, following the
    existing batch conventions. With ``indices``, those keys are visited and
    retained in exactly that order. Since the native forest requires
    nondecreasing endpoints, the corresponding prime values must also be
    nondecreasing. The packed endpoint domain is ``2 <= p <=
    min(2^32-1, LONG_MAX)``, where ``LONG_MAX`` is the largest signed value of
    the host's native C ``long``. Thus hosts with a 64-bit ``long`` support the
    full uint32 prime range, while narrower hosts classify larger endpoints
    as :class:`PackedP2UnsupportedDomain` before native marshalling.

    Each result is a :class:`PackedP2ForestResult`. Its immutable native-endian
    ``uint64`` buffer has logical shape ``(batch, 2, rows, columns)``. The
    planes are ordered ``P^0`` then ``P^1`` and each plane is row-major. The
    flat read-only ``coefficients`` memoryview uses format ``Q``; its linear
    index is ``(((batch_position*2 + plane)*rows + row)*columns + column)``.
    The ``primes``, ``endpoints``, and ``moduli`` tuples retain input order,
    with endpoints ``p`` and moduli ``p^2``. Each coefficient is checked to be
    canonical in ``[0,p^2)`` before checked GMP export. The bytes own a copy,
    so they remain valid after native temporaries are freed and after later
    forest calls.

    Raises :class:`PackedP2UnsupportedDomain` if a value is composite or
    outside the host's supported endpoint domain. A consumer can catch that
    exception and use :func:`remainder_forest_p2` for its arbitrary-precision
    fallback.
    """
    return _remainder_forest_ring(
        M,
        primes,
        None,
        2,
        kbase=kbase,
        indices=indices,
        V=V,
        kappa=kappa,
        return_state=False,
        fixed_p2=True,
        initial_z=None,
        packed_p2=True,
    )


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


def remainder_forest_pn_native(M, m, k, nP, kbase=0, indices=None, V=None,
                               kappa=None, return_state=False, *, z=None):
    r"""Return an owned arbitrary-precision PN result for compiled consumers.

    This opt-in API has the same input and endpoint semantics as
    :func:`remainder_forest_pn`, but returns :class:`NativePnForestResult`
    rather than constructing Sage polynomial matrices. Its capsule exposes a
    versioned, read-only GMP ABI; the result owner remains alive for as long as
    any capsule made by :meth:`NativePnForestResult.as_capsule` is alive.

    ``nP`` is the formal truncation bound and each value in ``m`` is its
    separate arbitrary-precision arithmetic modulus. ``endpoints`` in the
    result preserve the requested native endpoints; they need not equal a
    prime or modulus. ``return_state=True`` is unsupported because the native
    result contract covers forest outputs only.
    """
    if return_state:
        raise ValueError("native GMP output does not expose final_V state")
    return _remainder_forest_ring(
        M,
        m,
        k,
        nP,
        kbase=kbase,
        indices=indices,
        V=V,
        kappa=kappa,
        return_state=False,
        fixed_p2=False,
        initial_z=z,
        packed_p2=False,
        native_pn=True,
    )

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
