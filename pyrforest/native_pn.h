#ifndef PYRFOREST_NATIVE_PN_H
#define PYRFOREST_NATIVE_PN_H

#include <gmp.h>
#include <stddef.h>
#include <stdint.h>

#ifndef __GNU_MP_VERSION
#error "pyrforest PN GMP result ABI requires GMP headers"
#endif

#define PYRFOREST_PN_GMP_RESULT_ABI_VERSION 1
#define PYRFOREST_PN_GMP_RESULT_CAPSULE_NAME \
    "pyrforest.PnGmpForestResult/v1"

/*
 * Read-only, versioned view of a NativePnForestResult's owned GMP arrays.
 * Keep this prefix stable for ABI version 1. Advance the ABI version if a
 * field or its meaning changes; trailing fields may be added with a new
 * struct_size while preserving this prefix.
 */
typedef struct pyrforest_pn_gmp_result_v1 {
    uint32_t abi_version;
    uint32_t struct_size;
    uint32_t gmp_version_major;
    uint32_t gmp_version_minor;
    uint32_t gmp_version_patch;
    uint32_t gmp_limb_bits;
    uint32_t gmp_limb_bytes;
    size_t batch_count;
    size_t nP;
    size_t rows;
    size_t columns;
    const long *endpoints;
    mpz_srcptr moduli;
    mpz_srcptr coefficients;
} pyrforest_pn_gmp_result_v1;

static inline uint32_t pyrforest_pn_gmp_result_v1_build_major(void)
{
    return __GNU_MP_VERSION;
}

static inline uint32_t pyrforest_pn_gmp_result_v1_build_minor(void)
{
    return __GNU_MP_VERSION_MINOR;
}

static inline uint32_t pyrforest_pn_gmp_result_v1_build_patch(void)
{
    return __GNU_MP_VERSION_PATCHLEVEL;
}

static inline uint32_t pyrforest_pn_gmp_result_v1_build_limb_bits(void)
{
    return GMP_NUMB_BITS;
}

static inline uint32_t pyrforest_pn_gmp_result_v1_build_limb_bytes(void)
{
    return (uint32_t)sizeof(mp_limb_t);
}

static inline void pyrforest_pn_gmp_result_v1_init(
    pyrforest_pn_gmp_result_v1 *result,
    size_t batch_count,
    size_t nP,
    size_t rows,
    size_t columns,
    const long *endpoints,
    mpz_srcptr moduli,
    mpz_srcptr coefficients)
{
    result->abi_version = PYRFOREST_PN_GMP_RESULT_ABI_VERSION;
    result->struct_size = (uint32_t)sizeof(*result);
    result->gmp_version_major = __GNU_MP_VERSION;
    result->gmp_version_minor = __GNU_MP_VERSION_MINOR;
    result->gmp_version_patch = __GNU_MP_VERSION_PATCHLEVEL;
    result->gmp_limb_bits = GMP_NUMB_BITS;
    result->gmp_limb_bytes = (uint32_t)sizeof(mp_limb_t);
    result->batch_count = batch_count;
    result->nP = nP;
    result->rows = rows;
    result->columns = columns;
    result->endpoints = endpoints;
    result->moduli = moduli;
    result->coefficients = coefficients;
}

#endif /* PYRFOREST_NATIVE_PN_H */
