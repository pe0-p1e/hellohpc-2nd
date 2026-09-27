#pragma once

#include "main.hpp"

#if defined(__AVX512F__)
#include <immintrin.h>

struct miniclash_q1_result {
    uint32 q1;
    uint32 m0;
    uint32 m1;
    uint32 q17;
    uint32 q18;
    uint32 q19;
    uint32 q20;
    uint32 q21;
    uint32 block5;
};

struct miniclash_block0_q17_result {
    uint32 q17;
    uint32 q18;
    uint32 q19;
    uint32 q20;
};

static inline __attribute__((always_inline)) __m512i mc_ff(
    __m512i b, __m512i c, __m512i d)
{
    return _mm512_xor_si512(d, _mm512_and_si512(b, _mm512_xor_si512(c, d)));
}

static inline __attribute__((always_inline)) __m512i mc_gg(
    __m512i b, __m512i c, __m512i d)
{
    return _mm512_xor_si512(c, _mm512_and_si512(d, _mm512_xor_si512(b, c)));
}

template<bool NeedQ21>
static inline __attribute__((always_inline)) bool miniclash_find_q1_avx512(
    const uint32 q1base,
    const uint32 random_mask,
    const uint32 Q[68],
    const uint32 tt0,
    const uint32 tt1,
    const uint32 tt17,
    const uint32 tt18,
    const uint32 tt19,
    const uint32 q17_mask,
    const uint32 q17_expected,
    const uint32 q18_expected,
    const uint32 q19_expected,
    miniclash_q1_result& out)
{
    const uint32 q0  = Q[Qoff + 0];
    const uint32 qm1 = Q[Qoff - 1];
    const uint32 q2  = Q[Qoff + 2];
    const uint32 q15 = Q[Qoff + 15];
    const uint32 q16 = Q[Qoff + 16];

    uint32 block5 = 0;
    if constexpr (NeedQ21) {
        block5 = Q[Qoff + 6] - Q[Qoff + 5];
        block5 = RR(block5, 12)
            - FF(Q[Qoff + 5], Q[Qoff + 4], Q[Qoff + 3])
            - Q[Qoff + 2] - 0x4787c62au;
    }

    const __m512i vq1base = _mm512_set1_epi32((int)q1base);
    const __m512i vrmask = _mm512_set1_epi32((int)random_mask);
    const __m512i vq0 = _mm512_set1_epi32((int)q0);
    const __m512i vqm1 = _mm512_set1_epi32((int)qm1);
    const __m512i vq2 = _mm512_set1_epi32((int)q2);
    const __m512i vq15 = _mm512_set1_epi32((int)q15);
    const __m512i vq16 = _mm512_set1_epi32((int)q16);
    const __m512i vtt0 = _mm512_set1_epi32((int)tt0);
    const __m512i vtt1 = _mm512_set1_epi32((int)tt1);
    const __m512i vtt17 = _mm512_set1_epi32((int)tt17);
    const __m512i vtt18 = _mm512_set1_epi32((int)tt18);
    const __m512i vtt19 = _mm512_set1_epi32((int)tt19);

    const __m512i vq17mask = _mm512_set1_epi32((int)q17_mask);
    const __m512i vq17expected = _mm512_set1_epi32((int)q17_expected);
    const __m512i vq18mask = _mm512_set1_epi32((int)0xa0020000u);
    const __m512i vq18expected = _mm512_set1_epi32((int)q18_expected);
    const __m512i vq19mask = _mm512_set1_epi32((int)0x80020000u);
    const __m512i vq19expected = _mm512_set1_epi32((int)q19_expected);
    const __m512i vq20mask = _mm512_set1_epi32((int)0x80040000u);
    const __m512i vq20expected = _mm512_set1_epi32((int)0x00040000u);
    const __m512i bit17 = _mm512_set1_epi32((int)0x00020000u);

    alignas(64) uint32 rnd[16];
    alignas(64) uint32 lq1[16], lm0[16], lm1[16];
    alignas(64) uint32 lq17[16], lq18[16], lq19[16], lq20[16], lq21[16];

    for (unsigned batch = 0; batch < 256; ++batch) {
        if ((batch & 15u) == 0 && miniclash_cancelled())
            return false;

        for (unsigned lane = 0; lane < 16; ++lane)
            rnd[lane] = xrng64();

        const __m512i vrnd = _mm512_load_si512((const void*)rnd);
        const __m512i vq1 = _mm512_or_si512(
            vq1base, _mm512_and_si512(vrnd, vrmask));

        __m512i vm1 = _mm512_ror_epi32(_mm512_sub_epi32(vq2, vq1), 12);
        vm1 = _mm512_sub_epi32(vm1, mc_ff(vq1, vq0, vqm1));
        vm1 = _mm512_sub_epi32(vm1, vtt1);

        __m512i vq17 = _mm512_rol_epi32(_mm512_add_epi32(vtt17, vm1), 5);
        vq17 = _mm512_add_epi32(vq17, vq16);

        __mmask16 active = _mm512_cmpeq_epi32_mask(
            _mm512_and_si512(_mm512_xor_si512(vq17, vq16), vq17mask),
            vq17expected);
        active &= _mm512_testn_epi32_mask(vq17, bit17);
        if (!active)
            continue;

        __m512i vq18 = _mm512_add_epi32(mc_gg(vq17, vq16, vq15), vtt18);
        vq18 = _mm512_add_epi32(_mm512_rol_epi32(vq18, 9), vq17);

        active &= _mm512_cmpeq_epi32_mask(
            _mm512_and_si512(_mm512_xor_si512(vq18, vq17), vq18mask),
            vq18expected);
        if (!active)
            continue;

        __m512i vq19 = _mm512_add_epi32(mc_gg(vq18, vq17, vq16), vtt19);
        vq19 = _mm512_add_epi32(_mm512_rol_epi32(vq19, 14), vq18);

        active &= _mm512_cmpeq_epi32_mask(
            _mm512_and_si512(vq19, vq19mask), vq19expected);
        if (!active)
            continue;

        __m512i vm0 = _mm512_ror_epi32(_mm512_sub_epi32(vq1, vq0), 7);
        vm0 = _mm512_sub_epi32(vm0, vtt0);

        __m512i vq20 = _mm512_add_epi32(
            mc_gg(vq19, vq18, vq17),
            _mm512_add_epi32(vq16,
                _mm512_add_epi32(
                    _mm512_set1_epi32((int)0xe9b6c7aau), vm0)));
        vq20 = _mm512_add_epi32(_mm512_rol_epi32(vq20, 20), vq19);

        active &= _mm512_cmpeq_epi32_mask(
            _mm512_and_si512(_mm512_xor_si512(vq20, vq19), vq20mask),
            vq20expected);
        if (!active)
            continue;

        __m512i vq21 = _mm512_setzero_si512();
        if constexpr (NeedQ21) {
            vq21 = _mm512_add_epi32(
                mc_gg(vq20, vq19, vq18),
                _mm512_add_epi32(vq17,
                    _mm512_set1_epi32((int)(0xd62f105du + block5))));
            vq21 = _mm512_add_epi32(_mm512_rol_epi32(vq21, 5), vq20);

            active &= _mm512_testn_epi32_mask(
                _mm512_xor_si512(vq21, vq20),
                _mm512_set1_epi32((int)0x80020000u));
            if (!active)
                continue;
        }

        const unsigned lane = (unsigned)__builtin_ctz((unsigned)active);
        _mm512_store_si512((void*)lq1, vq1);
        _mm512_store_si512((void*)lm0, vm0);
        _mm512_store_si512((void*)lm1, vm1);
        _mm512_store_si512((void*)lq17, vq17);
        _mm512_store_si512((void*)lq18, vq18);
        _mm512_store_si512((void*)lq19, vq19);
        _mm512_store_si512((void*)lq20, vq20);
        if constexpr (NeedQ21)
            _mm512_store_si512((void*)lq21, vq21);

        out.q1 = lq1[lane];
        out.m0 = lm0[lane];
        out.m1 = lm1[lane];
        out.q17 = lq17[lane];
        out.q18 = lq18[lane];
        out.q19 = lq19[lane];
        out.q20 = lq20[lane];
        out.q21 = NeedQ21 ? lq21[lane] : 0;
        out.block5 = block5;
        return true;
    }

    return false;
}

static inline __attribute__((always_inline)) bool miniclash_find_block0_q17_avx512(
    const uint32 q16,
    const uint32 q15,
    const uint32 tt18,
    const uint32 tt19,
    const uint32 tt20,
    miniclash_block0_q17_result& out)
{
    const __m512i vq16 = _mm512_set1_epi32((int)q16);
    const __m512i vq15 = _mm512_set1_epi32((int)q15);
    const __m512i vtt18 = _mm512_set1_epi32((int)tt18);
    const __m512i vtt19 = _mm512_set1_epi32((int)tt19);
    const __m512i vtt20 = _mm512_set1_epi32((int)tt20);

    const __m512i vrmask = _mm512_set1_epi32((int)0x3ffd7ff7u);
    const __m512i vkeep = _mm512_set1_epi32((int)(q16 & 0xc0008008u));
    const __m512i vxor = _mm512_set1_epi32((int)0x40000000u);

    alignas(64) uint32 rnd[16];
    alignas(64) uint32 lq17[16], lq18[16], lq19[16], lq20[16];

    for (unsigned batch = 0; batch < 8; ++batch) {
        for (unsigned lane = 0; lane < 16; ++lane)
            rnd[lane] = xrng64();

        __m512i vq17 = _mm512_or_si512(
            _mm512_and_si512(_mm512_load_si512((const void*)rnd), vrmask),
            vkeep);
        vq17 = _mm512_xor_si512(vq17, vxor);

        __m512i vq18 = _mm512_add_epi32(mc_gg(vq17, vq16, vq15), vtt18);
        vq18 = _mm512_add_epi32(_mm512_rol_epi32(vq18, 9), vq17);
        __mmask16 active = _mm512_cmpeq_epi32_mask(
            _mm512_and_si512(_mm512_xor_si512(vq18, vq17),
                             _mm512_set1_epi32((int)0xa0020000u)),
            _mm512_set1_epi32((int)0x00020000u));
        if (!active)
            continue;

        __m512i vq19 = _mm512_add_epi32(mc_gg(vq18, vq17, vq16), vtt19);
        vq19 = _mm512_add_epi32(_mm512_rol_epi32(vq19, 14), vq18);
        active &= _mm512_cmpeq_epi32_mask(
            _mm512_and_si512(vq19, _mm512_set1_epi32((int)0x80020000u)),
            _mm512_set1_epi32((int)0x80000000u));
        if (!active)
            continue;

        __m512i vq20 = _mm512_add_epi32(mc_gg(vq19, vq18, vq17), vtt20);
        vq20 = _mm512_add_epi32(_mm512_rol_epi32(vq20, 20), vq19);
        active &= _mm512_cmpeq_epi32_mask(
            _mm512_and_si512(_mm512_xor_si512(vq20, vq19),
                             _mm512_set1_epi32((int)0x80040000u)),
            _mm512_set1_epi32((int)0x00040000u));
        if (!active)
            continue;

        const unsigned lane = (unsigned)__builtin_ctz((unsigned)active);
        _mm512_store_si512((void*)lq17, vq17);
        _mm512_store_si512((void*)lq18, vq18);
        _mm512_store_si512((void*)lq19, vq19);
        _mm512_store_si512((void*)lq20, vq20);

        out.q17 = lq17[lane];
        out.q18 = lq18[lane];
        out.q19 = lq19[lane];
        out.q20 = lq20[lane];
        return true;
    }

    return false;
}

#endif
