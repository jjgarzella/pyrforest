from libc.stddef cimport size_t
from libc.stdint cimport uint32_t
from sage.libs.gmp.types cimport mpz_srcptr

cdef extern from "native_pn.h":
    enum:
        PYRFOREST_PN_GMP_RESULT_ABI_VERSION
    const char *PYRFOREST_PN_GMP_RESULT_CAPSULE_NAME

    ctypedef struct pyrforest_pn_gmp_result_v1:
        uint32_t abi_version
        uint32_t struct_size
        uint32_t gmp_version_major
        uint32_t gmp_version_minor
        uint32_t gmp_version_patch
        uint32_t gmp_limb_bits
        uint32_t gmp_limb_bytes
        size_t batch_count
        size_t nP
        size_t rows
        size_t columns
        const long *endpoints
        mpz_srcptr moduli
        mpz_srcptr coefficients

    void pyrforest_pn_gmp_result_v1_init(
        pyrforest_pn_gmp_result_v1 *result,
        size_t batch_count,
        size_t nP,
        size_t rows,
        size_t columns,
        const long *endpoints,
        mpz_srcptr moduli,
        mpz_srcptr coefficients,
    )

    uint32_t pyrforest_pn_gmp_result_v1_build_major()
    uint32_t pyrforest_pn_gmp_result_v1_build_minor()
    uint32_t pyrforest_pn_gmp_result_v1_build_patch()
    uint32_t pyrforest_pn_gmp_result_v1_build_limb_bits()
    uint32_t pyrforest_pn_gmp_result_v1_build_limb_bytes()
