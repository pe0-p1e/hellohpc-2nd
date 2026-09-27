#include <iostream>
#include <vector>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif
#include "main.hpp"

thread_local uint32 seed32_1, seed32_2;
thread_local const std::atomic<bool>* miniclash_cancel_flag = nullptr;


static inline __attribute__((always_inline)) bool finish_block0_candidate(
    uint32 block[], const uint32 IV[], uint32 a, uint32 b, uint32 c, uint32 d)
{
    MD5_STEP(II, a, b, c, d, block[12], 0x655b59c3, 6);
    if (0 != ((a^c) >> 31)) return false;
    MD5_STEP(II, d, a, b, c, block[3], 0x8f0ccc92, 10);
    if (0 != ((b^d) >> 31)) return false;
    MD5_STEP(II, c, d, a, b, block[10], 0xffeff47d, 15);
    if (0 != ((a^c) >> 31)) return false;
    MD5_STEP(II, b, c, d, a, block[1], 0x85845dd1, 21);
    if (0 != ((b^d) >> 31)) return false;
    MD5_STEP(II, a, b, c, d, block[8], 0x6fa87e4f, 6);
    if (0 != ((a^c) >> 31)) return false;
    MD5_STEP(II, d, a, b, c, block[15], 0xfe2ce6e0, 10);
    if (0 != ((b^d) >> 31)) return false;
    MD5_STEP(II, c, d, a, b, block[6], 0xa3014314, 15);
    if (0 != ((a^c) >> 31)) return false;
    MD5_STEP(II, b, c, d, a, block[13], 0x4e0811a1, 21);
    if (0 == ((b^d) >> 31)) return false;
    MD5_STEP(II, a, b, c, d, block[4], 0xf7537e82, 6);
    if (0 != ((a^c) >> 31)) return false;
    MD5_STEP(II, d, a, b, c, block[11], 0xbd3af235, 10);
    if (0 != ((b^d) >> 31)) return false;
    MD5_STEP(II, c, d, a, b, block[2], 0x2ad7d2bb, 15);
    if (0 != ((a^c) >> 31)) return false;
    MD5_STEP(II, b, c, d, a, block[9], 0xeb86d391, 21);

    uint32 IHV1 = b + IV[1];
    uint32 IHV2 = c + IV[2];
    uint32 IHV3 = d + IV[3];

    bool wang = true;
    if (0x02000000 != ((IHV2^IHV1) & 0x86000000)) wang = false;
    if (0 != ((IHV1^IHV3) & 0x82000000)) wang = false;
    if (0 != (IHV1 & 0x06000020)) wang = false;

    bool stevens = true;
    if (((IHV1^IHV2)>>31)!=0 || ((IHV1^IHV3)>>31)!=0) stevens = false;
    if ((IHV3&(1<<25))!=0 || (IHV2&(1<<25))!=0 || (IHV1&(1<<25))!=0
        || ((IHV2^IHV1)&1)!=0) stevens = false;
    if (!(wang || stevens)) return false;

    uint32 V1[4], V2[4];
    for (int t = 0; t < 4; ++t) V2[t] = V1[t] = IV[t];
    uint32 block2[16];
    for (int t = 0; t < 16; ++t) block2[t] = block[t];
    block2[4] += 1u<<31;
    block2[11] += 1u<<15;
    block2[14] += 1u<<31;
    md5_compress(V1, block);
    md5_compress(V2, block2);
    return (V2[0] == V1[0] + (1u<<31))
        && (V2[1] == V1[1] + (1u<<31) + (1u<<25))
        && (V2[2] == V1[2] + (1u<<31) + (1u<<25))
        && (V2[3] == V1[3] + (1u<<31) + (1u<<25));
}

namespace {
const std::vector<uint32> q4mask = [] {
	std::vector<uint32> v(1<<4);
	for (unsigned k = 0; k < v.size(); ++k)
		v[k] = ((k<<2) ^ (k<<26)) & 0x38000004;
	return v;
}();

const std::vector<uint32> q9q10mask = [] {
	std::vector<uint32> v(1<<3);
	for (unsigned k = 0; k < v.size(); ++k)
		v[k] = ((k<<13) ^ (k<<4)) & 0x2060;
	return v;
}();

const std::vector<uint32> q9mask = [] {
	std::vector<uint32> v(1<<16);
	for (unsigned k = 0; k < v.size(); ++k)
		v[k] = ((k<<1) ^ (k<<2) ^ (k<<5) ^ (k<<7) ^ (k<<8) ^ (k<<10) ^ (k<<11) ^ (k<<13)) & 0x0eb94f16;
	return v;
}();
}

void find_block0(uint32 block[], const uint32 IV[])
{
	uint32 Q[68] = { IV[0], IV[3], IV[2], IV[1] };

	while (true)
	{
		Q[Qoff + 1] = xrng64();
		Q[Qoff + 3] = (xrng64() & 0xfe87bc3f) | 0x017841c0;
		Q[Qoff + 4] = (xrng64() & 0x44000033) | 0x000002c0 | (Q[Qoff + 3] & 0x0287bc00);
		Q[Qoff + 5] = 0x41ffffc8 | (Q[Qoff + 4] & 0x04000033);
		Q[Qoff + 6] = 0xb84b82d6;
		Q[Qoff + 7] = (xrng64() & 0x68000084) | 0x02401b43;
		Q[Qoff + 8] = (xrng64() & 0x2b8f6e04) | 0x005090d3 | (~Q[Qoff + 7] & 0x40000000);
		Q[Qoff + 9] = 0x20040068 | (Q[Qoff + 8] & 0x00020000) | (~Q[Qoff + 8] & 0x40000000);
		Q[Qoff + 10] = (xrng64() & 0x40000000) | 0x1040b089;
		Q[Qoff + 11] = (xrng64() & 0x10408008) | 0x0fbb7f16 | (~Q[Qoff + 10] & 0x40000000);
		Q[Qoff + 12] = (xrng64() & 0x1ed9df7f) | 0x00022080 | (~Q[Qoff + 11] & 0x40200000);
		Q[Qoff + 13] = (xrng64() & 0x5efb4f77) | 0x20049008;
		Q[Qoff + 14] = (xrng64() & 0x1fff5f77) | 0x0000a088 | (~Q[Qoff + 13] & 0x40000000);
		Q[Qoff + 15] = (xrng64() & 0x5efe7ff7) | 0x80008000 | (~Q[Qoff + 14] & 0x00010000);
		Q[Qoff + 16] = (xrng64() & 0x1ffdffff) | 0xa0000000 | (~Q[Qoff + 15] & 0x40020000);

		MD5_REVERSE_STEP(0, 0xd76aa478, 7);
		MD5_REVERSE_STEP(6, 0xa8304613, 17);
		MD5_REVERSE_STEP(7, 0xfd469501, 22);
		MD5_REVERSE_STEP(11, 0x895cd7be, 22);
		MD5_REVERSE_STEP(14, 0xa679438e, 17);
		MD5_REVERSE_STEP(15, 0x49b40821, 22);

		const uint32 tt1 = FF(Q[Qoff + 1], Q[Qoff + 0], Q[Qoff - 1]) + Q[Qoff - 2] + 0xe8c7b756;		
		const uint32 tt17 = GG(Q[Qoff + 16], Q[Qoff + 15], Q[Qoff + 14]) + Q[Qoff + 13] + 0xf61e2562;
		const uint32 tt18 = Q[Qoff + 14] + 0xc040b340 + block[6];
		const uint32 tt19 = Q[Qoff + 15] + 0x265e5a51 + block[11];
		const uint32 tt20 = Q[Qoff + 16] + 0xe9b6c7aa + block[0];
		const uint32 tt5 = RR(Q[Qoff + 6] - Q[Qoff + 5], 12) - FF(Q[Qoff + 5], Q[Qoff + 4], Q[Qoff + 3]) - 0x4787c62a;

		// change q17 until conditions are met on q18, q19 and q20
		unsigned counter = 0;
		while (counter < (1 << 7))
		{
			const uint32 q16 = Q[Qoff + 16];
			uint32 q17 = ((xrng64() & 0x3ffd7ff7) | (q16&0xc0008008)) ^ 0x40000000;
			++counter;

			uint32 q18 = GG(q17, q16, Q[Qoff + 15]) + tt18;
			q18 = RL(q18, 9); q18 += q17;
			if (0x00020000 != ((q18^q17)&0xa0020000))
				continue;

			uint32 q19 = GG(q18, q17, q16) + tt19;
			q19 = RL(q19, 14); q19 += q18;
			if (0x80000000 != (q19 & 0x80020000))
				continue;
			
			uint32 q20 = GG(q19, q18, q17) + tt20;
			q20 = RL(q20, 20); q20 += q19;
			if (0x00040000 != ((q20^q19) & 0x80040000))
				continue;

			block[1] = q17-q16; block[1] = RR(block[1], 5); block[1] -= tt17;
			uint32 q2 = block[1] + tt1; q2 = RL(q2, 12); q2 += Q[Qoff + 1];
			block[5] = tt5 - q2;

			Q[Qoff + 2] = q2;
			Q[Qoff + 17] = q17;
			Q[Qoff + 18] = q18;
			Q[Qoff + 19] = q19;
			Q[Qoff + 20] = q20;
			MD5_REVERSE_STEP(2, 0x242070db, 17);

			counter = 0;
			break;
		}
		if (counter != 0)
			continue;

		const uint32 q4 = Q[Qoff + 4];
		const uint32 q9backup = Q[Qoff + 9];
		const uint32 tt21 = GG(Q[Qoff+20], Q[Qoff+19], Q[Qoff+18]) + Q[Qoff+17] + 0xd62f105d;

		// iterate over possible changes of q4 
		// while keeping all conditions on q1-q20 intact
		// this changes m3, m4, m5 and m7
		unsigned counter2 = 0;
		while (counter2 < (1<<4))
		{
			Q[Qoff+4] = q4 ^ q4mask[counter2];
			++counter2;
			MD5_REVERSE_STEP(5, 0x4787c62a, 12);
			uint32 q21 = tt21 + block[5];
			q21 = RL(q21,5); q21 += Q[Qoff+20];
			if (0 != ((q21^Q[Qoff+20]) & 0x80020000))
				continue;

			Q[Qoff + 21] = q21;
			MD5_REVERSE_STEP(3, 0xc1bdceee, 22);
			MD5_REVERSE_STEP(4, 0xf57c0faf, 7);
			MD5_REVERSE_STEP(7, 0xfd469501, 22);
			
			const uint32 tt22 = GG(Q[Qoff + 21], Q[Qoff + 20], Q[Qoff + 19]) + Q[Qoff + 18] + 0x02441453;
			const uint32 tt23 = Q[Qoff + 19] + 0xd8a1e681 + block[15];
			const uint32 tt24 = Q[Qoff + 20] + 0xe7d3fbc8 + block[4];

			const uint32 tt9 = Q[Qoff + 6] + 0x8b44f7af;
			const uint32 tt10 = Q[Qoff + 7] + 0xffff5bb1;
			const uint32 tt8 = FF(Q[Qoff + 8], Q[Qoff + 7], Q[Qoff + 6]) + Q[Qoff + 5] + 0x698098d8;		
			const uint32 tt12 = RR(Q[Qoff+13]-Q[Qoff+12],7) - 0x6b901122;
			const uint32 tt13 = RR(Q[Qoff+14]-Q[Qoff+13],12) - FF(Q[Qoff+13],Q[Qoff+12],Q[Qoff+11]) - 0xfd987193;

			// iterate over possible changes of q9 and q10
			// while keeping conditions on q1-q21 intact
			// this changes m8, m9, m10, m12 and m13 (and not m11!)
			// the possible changes of q9 that also do not change m10 are used below
			for (unsigned counter3 = 0; counter3 < (1<<3);)
			{
				uint32 q10 = Q[Qoff+10] ^ (q9q10mask[counter3] & 0x60);
				Q[Qoff + 9] = q9backup ^ (q9q10mask[counter3] & 0x2000);
				++counter3;
				uint32 m10 = RR(Q[Qoff+11]-q10,17);
				m10 -= FF(q10, Q[Qoff+9], Q[Qoff+8]) + tt10;

				uint32 aa = Q[Qoff + 21];
				uint32 dd = tt22+m10; dd = RL(dd, 9) + aa;
				if (0x80000000 != (dd & 0x80000000)) continue;

				uint32 bb = Q[Qoff + 20];
				uint32 cc = tt23 + GG(dd, aa, bb); 
				if (0 != (cc & 0x20000)) continue;
				cc = RL(cc, 14) + dd;
				if (0 != (cc & 0x80000000)) continue;

				bb = tt24 + GG(cc, dd, aa); bb = RL(bb, 20) + cc;
				if (0 == (bb & 0x80000000)) continue;

				block[10] = m10;
				block[13] = tt13 - q10;

				// iterate over possible changes of q9
				// while keeping intact conditions on q1-q24
				// this changes m8, m9 and m12 (but not m10!)
                // SIMD-filter q9 candidates through the first ten MD5 steps.
#if defined(__AVX512F__)
                const __m512i vaa0 = _mm512_set1_epi32((int)aa);
                const __m512i vbb0 = _mm512_set1_epi32((int)bb);
                const __m512i vcc0 = _mm512_set1_epi32((int)cc);
                const __m512i vdd0 = _mm512_set1_epi32((int)dd);
                const __m512i vq9base = _mm512_set1_epi32((int)Q[Qoff + 9]);
                const __m512i vq8 = _mm512_set1_epi32((int)Q[Qoff + 8]);
                const __m512i vq7 = _mm512_set1_epi32((int)Q[Qoff + 7]);
                const __m512i vq12 = _mm512_set1_epi32((int)Q[Qoff + 12]);
                const __m512i vq11 = _mm512_set1_epi32((int)Q[Qoff + 11]);
                const __m512i vq10 = _mm512_set1_epi32((int)q10);
                const __m512i vtt8 = _mm512_set1_epi32((int)tt8);
                const __m512i vtt9 = _mm512_set1_epi32((int)tt9);
                const __m512i vtt12 = _mm512_set1_epi32((int)tt12);
                const __m512i vc15 = _mm512_set1_epi32(1 << 15);

                alignas(64) uint32 lane_q9[16], lane_b8[16], lane_b9[16], lane_b12[16];
                alignas(64) uint32 lane_a[16], lane_b[16], lane_c[16], lane_d[16];

#define VGG(x,y,z) _mm512_xor_si512((y), _mm512_and_si512((z), _mm512_xor_si512((x),(y))))
#define VHH(x,y,z) _mm512_xor_si512(_mm512_xor_si512((x),(y)),(z))
#define VII(x,y,z) _mm512_xor_si512((y), _mm512_or_si512((x), _mm512_xor_si512((z), _mm512_set1_epi32(-1))))
#define VSTEP_GG(A,B,C,D,M,K,R) do { \
    (A) = _mm512_add_epi32((A), _mm512_add_epi32(VGG((B),(C),(D)), _mm512_add_epi32((M), _mm512_set1_epi32((int)(K))))); \
    (A) = _mm512_add_epi32(_mm512_rol_epi32((A),(R)), (B)); \
} while(0)
#define VSTEP_HH(A,B,C,D,M,K,R) do { \
    (A) = _mm512_add_epi32((A), _mm512_add_epi32(VHH((B),(C),(D)), _mm512_add_epi32((M), _mm512_set1_epi32((int)(K))))); \
    (A) = _mm512_add_epi32(_mm512_rol_epi32((A),(R)), (B)); \
} while(0)
#define VSTEP_II(A,B,C,D,M,K,R) do { \
    (A) = _mm512_add_epi32((A), _mm512_add_epi32(VII((B),(C),(D)), _mm512_add_epi32((M), _mm512_set1_epi32((int)(K))))); \
    (A) = _mm512_add_epi32(_mm512_rol_epi32((A),(R)), (B)); \
} while(0)

                for (unsigned counter4 = 0; counter4 < (1u<<16); counter4 += 16)
                {
                    __m512i q9v = _mm512_xor_si512(vq9base,
                        _mm512_loadu_si512((const void*)(&q9mask[counter4])));
                    __m512i b12v = _mm512_sub_epi32(
                        _mm512_sub_epi32(vtt12, _mm512_xor_si512(vq10,
                            _mm512_and_si512(vq12, _mm512_xor_si512(vq11, vq10)))), q9v);
                    __m512i m8v = _mm512_sub_epi32(q9v, vq8);
                    __m512i b8v = _mm512_sub_epi32(_mm512_ror_epi32(m8v, 7), vtt8);
                    __m512i m9v = _mm512_sub_epi32(vq10, q9v);
                    __m512i ff9 = _mm512_xor_si512(vq7,
                        _mm512_and_si512(q9v, _mm512_xor_si512(vq8, vq7)));
                    __m512i b9v = _mm512_sub_epi32(
                        _mm512_sub_epi32(_mm512_ror_epi32(m9v, 12), ff9), vtt9);

                    __m512i a = vaa0, b = vbb0, c = vcc0, d = vdd0;
                    VSTEP_GG(a,b,c,d,b9v,0x21e1cde6,5);
                    VSTEP_GG(d,a,b,c,_mm512_set1_epi32((int)block[14]),0xc33707d6,9);
                    VSTEP_GG(c,d,a,b,_mm512_set1_epi32((int)block[3]),0xf4d50d87,14);
                    VSTEP_GG(b,c,d,a,b8v,0x455a14ed,20);
                    VSTEP_GG(a,b,c,d,_mm512_set1_epi32((int)block[13]),0xa9e3e905,5);
                    VSTEP_GG(d,a,b,c,_mm512_set1_epi32((int)block[2]),0xfcefa3f8,9);
                    VSTEP_GG(c,d,a,b,_mm512_set1_epi32((int)block[7]),0x676f02d9,14);
                    VSTEP_GG(b,c,d,a,b12v,0x8d2a4c8a,20);
                    VSTEP_HH(a,b,c,d,_mm512_set1_epi32((int)block[5]),0xfffa3942,4);
                    VSTEP_HH(d,a,b,c,b8v,0x8771f681,11);
                    c = _mm512_add_epi32(c, _mm512_add_epi32(VHH(d,a,b),
                        _mm512_set1_epi32((int)(block[11] + 0x6d9d6122u))));
                    const __mmask16 survivors = _mm512_testn_epi32_mask(c, vc15);
                    if (!survivors) continue;
                    c = _mm512_add_epi32(_mm512_rol_epi32(c,16), d);
                    VSTEP_HH(b,c,d,a,_mm512_set1_epi32((int)block[14]),0xfde5380c,23);
                    VSTEP_HH(a,b,c,d,_mm512_set1_epi32((int)block[1]),0xa4beea44,4);
                    VSTEP_HH(d,a,b,c,_mm512_set1_epi32((int)block[4]),0x4bdecfa9,11);
                    VSTEP_HH(c,d,a,b,_mm512_set1_epi32((int)block[7]),0xf6bb4b60,16);
                    VSTEP_HH(b,c,d,a,_mm512_set1_epi32((int)block[10]),0xbebfbc70,23);
                    VSTEP_HH(a,b,c,d,_mm512_set1_epi32((int)block[13]),0x289b7ec6,4);
                    VSTEP_HH(d,a,b,c,_mm512_set1_epi32((int)block[0]),0xeaa127fa,11);
                    VSTEP_HH(c,d,a,b,_mm512_set1_epi32((int)block[3]),0xd4ef3085,16);
                    VSTEP_HH(b,c,d,a,_mm512_set1_epi32((int)block[6]),0x04881d05,23);
                    VSTEP_HH(a,b,c,d,b9v,0xd9d4d039,4);
                    VSTEP_HH(d,a,b,c,b12v,0xe6db99e5,11);
                    VSTEP_HH(c,d,a,b,_mm512_set1_epi32((int)block[15]),0x1fa27cf8,16);
                    VSTEP_HH(b,c,d,a,_mm512_set1_epi32((int)block[2]),0xc4ac5665,23);
                    __mmask16 active2 = survivors & _mm512_testn_epi32_mask(
                        _mm512_xor_si512(b,d), _mm512_set1_epi32((int)0x80000000u));
                    if (!active2) continue;
                    const __m512i vsign = _mm512_set1_epi32((int)0x80000000u);
                    VSTEP_II(a,b,c,d,_mm512_set1_epi32((int)block[0]),0xf4292244,6);
                    active2 &= _mm512_testn_epi32_mask(_mm512_xor_si512(a,c), vsign);
                    if (!active2) continue;
                    VSTEP_II(d,a,b,c,_mm512_set1_epi32((int)block[7]),0x432aff97,10);
                    active2 &= _mm512_test_epi32_mask(_mm512_xor_si512(b,d), vsign);
                    if (!active2) continue;
                    VSTEP_II(c,d,a,b,_mm512_set1_epi32((int)block[14]),0xab9423a7,15);
                    active2 &= _mm512_testn_epi32_mask(_mm512_xor_si512(a,c), vsign);
                    if (!active2) continue;
                    VSTEP_II(b,c,d,a,_mm512_set1_epi32((int)block[5]),0xfc93a039,21);
                    active2 &= _mm512_testn_epi32_mask(_mm512_xor_si512(b,d), vsign);
                    if (!active2) continue;

                    _mm512_store_si512((void*)lane_q9,q9v);
                    _mm512_store_si512((void*)lane_b8,b8v);
                    _mm512_store_si512((void*)lane_b9,b9v);
                    _mm512_store_si512((void*)lane_b12,b12v);
                    _mm512_store_si512((void*)lane_a,a);
                    _mm512_store_si512((void*)lane_b,b);
                    _mm512_store_si512((void*)lane_c,c);
                    _mm512_store_si512((void*)lane_d,d);

                    unsigned sm = (unsigned)active2;
                    while (sm) {
                        const unsigned lane = __builtin_ctz(sm);
                        sm &= sm - 1;
                        Q[Qoff + 9] = lane_q9[lane];
                        block[8] = lane_b8[lane];
                        block[9] = lane_b9[lane];
                        block[12] = lane_b12[lane];
                        if (finish_block0_candidate(block, IV,
                            lane_a[lane], lane_b[lane], lane_c[lane], lane_d[lane]))
                            return;
                    }
                }
#undef VSTEP_II
#undef VSTEP_HH
#undef VSTEP_GG
#undef VII
#undef VHH
#undef VGG
#else
                for (unsigned counter4 = 0; counter4 < (1u<<16); ++counter4)
                {
                    uint32 q9 = Q[Qoff + 9] ^ q9mask[counter4];
                    block[12] = tt12 - FF(Q[Qoff + 12], Q[Qoff + 11], q10) - q9;
                    uint32 m8 = q9 - Q[Qoff + 8];
                    block[8] = RR(m8, 7) - tt8;
                    uint32 m9 = q10 - q9;
                    block[9] = RR(m9, 12) - FF(q9, Q[Qoff + 8], Q[Qoff + 7]) - tt9;

                    uint32 a = aa, b = bb, c = cc, d = dd;
                    MD5_STEP(GG, a, b, c, d, block[9], 0x21e1cde6, 5);
                    MD5_STEP(GG, d, a, b, c, block[14], 0xc33707d6, 9);
                    MD5_STEP(GG, c, d, a, b, block[3], 0xf4d50d87, 14);
                    MD5_STEP(GG, b, c, d, a, block[8], 0x455a14ed, 20);
                    MD5_STEP(GG, a, b, c, d, block[13], 0xa9e3e905, 5);
                    MD5_STEP(GG, d, a, b, c, block[2], 0xfcefa3f8, 9);
                    MD5_STEP(GG, c, d, a, b, block[7], 0x676f02d9, 14);
                    MD5_STEP(GG, b, c, d, a, block[12], 0x8d2a4c8a, 20);
                    MD5_STEP(HH, a, b, c, d, block[5], 0xfffa3942, 4);
                    MD5_STEP(HH, d, a, b, c, block[8], 0x8771f681, 11);
                    c += HH(d, a, b) + block[11] + 0x6d9d6122;
                    if (0 != (c & (1 << 15))) continue;
                    c = (c<<16 | c>>16) + d;
                    MD5_STEP(HH, b, c, d, a, block[14], 0xfde5380c, 23);
                    MD5_STEP(HH, a, b, c, d, block[1], 0xa4beea44, 4);
                    MD5_STEP(HH, d, a, b, c, block[4], 0x4bdecfa9, 11);
                    MD5_STEP(HH, c, d, a, b, block[7], 0xf6bb4b60, 16);
                    MD5_STEP(HH, b, c, d, a, block[10], 0xbebfbc70, 23);
                    MD5_STEP(HH, a, b, c, d, block[13], 0x289b7ec6, 4);
                    MD5_STEP(HH, d, a, b, c, block[0], 0xeaa127fa, 11);
                    MD5_STEP(HH, c, d, a, b, block[3], 0xd4ef3085, 16);
                    MD5_STEP(HH, b, c, d, a, block[6], 0x04881d05, 23);
                    MD5_STEP(HH, a, b, c, d, block[9], 0xd9d4d039, 4);
                    MD5_STEP(HH, d, a, b, c, block[12], 0xe6db99e5, 11);
                    MD5_STEP(HH, c, d, a, b, block[15], 0x1fa27cf8, 16);
                    MD5_STEP(HH, b, c, d, a, block[2], 0xc4ac5665, 23);
                    if (0 != ((b^d) & 0x80000000)) continue;
                    MD5_STEP(II, a, b, c, d, block[0], 0xf4292244, 6);
                    if (0 != ((a^c) >> 31)) continue;
                    MD5_STEP(II, d, a, b, c, block[7], 0x432aff97, 10);
                    if (0 == ((b^d) >> 31)) continue;
                    MD5_STEP(II, c, d, a, b, block[14], 0xab9423a7, 15);
                    if (0 != ((a^c) >> 31)) continue;
                    MD5_STEP(II, b, c, d, a, block[5], 0xfc93a039, 21);
                    if (0 != ((b^d) >> 31)) continue;
                    if (finish_block0_candidate(block, IV, a,b,c,d)) return;
                }
#endif
			}
		}
	}
}
