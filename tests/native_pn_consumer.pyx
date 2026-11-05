# cython: language_level=3
"""Separate Cython consumer used to validate the public PN GMP ABI."""

from cpython.bytes cimport PyBytes_FromString
from cpython.pycapsule cimport PyCapsule_GetPointer
from libc.stddef cimport size_t
from libc.stdlib cimport free, malloc
from sage.libs.gmp.types cimport mpz_srcptr
from pyrforest.native_pn cimport (
    PYRFOREST_PN_GMP_RESULT_ABI_VERSION,
    PYRFOREST_PN_GMP_RESULT_CAPSULE_NAME,
    pyrforest_pn_gmp_result_v1,
    pyrforest_pn_gmp_result_v1_build_major,
    pyrforest_pn_gmp_result_v1_build_minor,
    pyrforest_pn_gmp_result_v1_build_patch,
    pyrforest_pn_gmp_result_v1_build_limb_bits,
    pyrforest_pn_gmp_result_v1_build_limb_bytes,
)

cdef extern from "gmp.h":
    size_t mpz_sizeinbase(mpz_srcptr value, int base)
    char *mpz_get_str(char *text, int base, mpz_srcptr value)
    int mpz_cmp(mpz_srcptr left, mpz_srcptr right)


cdef bytes _decimal(mpz_srcptr value):
    cdef size_t capacity = mpz_sizeinbase(value, 10) + 3
    cdef char *text = <char *>malloc(capacity)
    cdef bytes result
    if text == NULL:
        raise MemoryError("unable to allocate GMP decimal scratch buffer")
    try:
        if mpz_get_str(text, 10, value) == NULL:
            raise MemoryError("GMP could not format a result value")
        result = PyBytes_FromString(text)
    finally:
        free(text)
    return result


def read_capsule(object capsule):
    """Validate a capsule and copy its metadata/values for tests only."""
    cdef const pyrforest_pn_gmp_result_v1 *payload
    cdef size_t max_size = <size_t>-1
    cdef size_t total
    cdef size_t batch_position
    cdef size_t row
    cdef size_t column
    cdef size_t degree
    cdef size_t offset
    cdef mpz_srcptr modulus_pointer
    cdef mpz_srcptr coefficient_pointer
    cdef object moduli = []
    cdef object coefficients = []
    cdef object endpoints = []
    cdef size_t bits

    payload = <const pyrforest_pn_gmp_result_v1 *>PyCapsule_GetPointer(
        capsule, PYRFOREST_PN_GMP_RESULT_CAPSULE_NAME
    )
    if payload == NULL:
        raise ValueError("capsule has no PN GMP result payload")
    if payload.abi_version != PYRFOREST_PN_GMP_RESULT_ABI_VERSION:
        raise ValueError("unsupported PN GMP result ABI version")
    if payload.struct_size < sizeof(pyrforest_pn_gmp_result_v1):
        raise ValueError("PN GMP result payload is shorter than the v1 prefix")
    if (payload.gmp_version_major != pyrforest_pn_gmp_result_v1_build_major()
            or payload.gmp_version_minor != pyrforest_pn_gmp_result_v1_build_minor()
            or payload.gmp_version_patch != pyrforest_pn_gmp_result_v1_build_patch()
            or payload.gmp_limb_bits != pyrforest_pn_gmp_result_v1_build_limb_bits()
            or payload.gmp_limb_bytes != pyrforest_pn_gmp_result_v1_build_limb_bytes()):
        raise ValueError("PN GMP result payload does not match consumer GMP ABI")

    total = payload.batch_count
    for bits in (payload.nP, payload.rows, payload.columns):
        if bits and total > max_size // bits:
            raise OverflowError("PN GMP result dimensions overflow size_t")
        total *= bits
    if payload.batch_count and (
            payload.endpoints == NULL or payload.moduli == NULL):
        raise ValueError("PN GMP result metadata arrays are missing")
    if total and payload.coefficients == NULL:
        raise ValueError("PN GMP result coefficient array is missing")

    for batch_position in range(payload.batch_count):
        endpoints.append(payload.endpoints[batch_position])
        modulus_pointer = <mpz_srcptr>(payload.moduli + batch_position)
        if modulus_pointer[0]._mp_size <= 0:
            raise ValueError("PN GMP result modulus must be positive")
        moduli.append(_decimal(modulus_pointer))

    for batch_position in range(payload.batch_count):
        modulus_pointer = <mpz_srcptr>(payload.moduli + batch_position)
        for degree in range(payload.nP):
            for row in range(payload.rows):
                for column in range(payload.columns):
                    offset = (((batch_position * payload.nP + degree)
                               * payload.rows + row) * payload.columns + column)
                    coefficient_pointer = <mpz_srcptr>(
                        payload.coefficients + offset
                    )
                    if (coefficient_pointer[0]._mp_size < 0
                            or mpz_cmp(coefficient_pointer, modulus_pointer) >= 0):
                        raise ValueError("PN GMP coefficient is not canonical")
                    coefficients.append(_decimal(coefficient_pointer))

    return {
        "abi_version": payload.abi_version,
        "struct_size": payload.struct_size,
        "gmp": (
            payload.gmp_version_major,
            payload.gmp_version_minor,
            payload.gmp_version_patch,
            payload.gmp_limb_bits,
            payload.gmp_limb_bytes,
        ),
        "shape": (
            payload.batch_count,
            payload.nP,
            payload.rows,
            payload.columns,
        ),
        "endpoints": endpoints,
        "moduli": moduli,
        "coefficients": coefficients,
    }


def read_result(object result):
    """Hold the owner and capsule through a full compiled read."""
    cdef object capsule = result.as_capsule()
    return read_capsule(capsule)
