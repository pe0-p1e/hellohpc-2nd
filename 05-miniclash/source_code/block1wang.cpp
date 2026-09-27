#include <iostream>
#include <vector>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif
#include "main.hpp"


static inline __attribute__((always_inline)) bool finish_wang_candidate(
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

    uint32 block2[16];
    uint32 V1[4], V2[4];
    for (int t=0;t<4;++t) { V1[t]=IV[t]; V2[t]=IV[t]+(1u<<31); }
    V2[1] += (1u<<25); V2[2] += (1u<<25); V2[3] += (1u<<25);
    for (int t=0;t<16;++t) block2[t]=block[t];
    block2[4] += 1u<<31; block2[11] -= 1u<<15; block2[14] += 1u<<31;
    md5_compress(V1,block); md5_compress(V2,block2);
    return V2[0]==V1[0] && V2[1]==V1[1] && V2[2]==V1[2] && V2[3]==V1[3];
}

namespace {
const std::vector<uint32> q4mask = [] {
	std::vector<uint32> v(1<<6);
	for (unsigned k = 0; k < v.size(); ++k)
		v[k] = ((k<<13) ^ (k<<19)) & 0x01c0e000;
	return v;
}();

struct WangQMasks {
	std::vector<uint32> q9, q10;
	WangQMasks(): q9(1<<5), q10(1<<5) {
		for (unsigned k = 0; k < q9.size(); ++k) {
			uint32 msk = (k<<5) ^ (k<<13) ^ (k<<17) ^ (k<<24);
			q9[k] = msk & 0x00084000;
			q10[k] = msk & 0x18000020;
		}
	}
};
const WangQMasks wang_q_masks;

const std::vector<uint32> q9mask2 = [] {
	std::vector<uint32> v(1<<10);
	for (unsigned k = 0; k < v.size(); ++k)
		v[k] = ((k<<1) ^ (k<<7) ^ (k<<14) ^ (k<<15) ^ (k<<22)) & 0x6074041c;
	return v;
}();
}

void find_block1_wang(uint32 block[], const uint32 IV[])
{
	uint32 Q[68] = { IV[0], IV[3], IV[2], IV[1] };

	while (true) 
	{
		if (miniclash_cancelled()) return;
		uint32 aa = Q[Qoff] & 0x80000000;
		uint32 bb = 0x80000000 ^ aa;

		Q[Qoff + 2] = (xrng64() & 0x71de7799) | 0x0c008840 | bb;
		Q[Qoff + 3] = (xrng64() & 0x01c06601) | 0x3e1f0966 | (Q[Qoff + 2] & 0x80000018);
		Q[Qoff + 4] = 0x3a040010 | (Q[Qoff + 3] & 0x80000601);
		Q[Qoff + 5] = (xrng64() & 0x03c0e000) | 0x482f0e50 | aa;
		Q[Qoff + 6] = (xrng64() & 0x600c0000) | 0x05e2ec56 | aa;
		Q[Qoff + 7] = (xrng64() & 0x604c203e) | 0x16819e01 | bb | (Q[Qoff + 6] & 0x01000000);
		Q[Qoff + 8] = (xrng64() & 0x604c7c1c) | 0x043283e0 | (Q[Qoff + 7] & 0x80000002);
		Q[Qoff + 9] =  (xrng64() & 0x00002800) | 0x1c0101c1 | (Q[Qoff + 8] & 0x80001000);
		Q[Qoff + 10] = 0x078bcbc0 | bb;
		Q[Qoff + 11] = (xrng64() & 0x07800000) | 0x607dc7df | bb;
		Q[Qoff + 12] = (xrng64() & 0x00f00f7f) | 0x00081080 | (Q[Qoff + 11] & 0xe7000000);
		Q[Qoff + 13] = (xrng64() & 0x00701f77) | 0x3f0fe008 | aa;
		Q[Qoff + 14] = (xrng64() & 0x00701f77) | 0x408be088 | aa;
		Q[Qoff + 15] = (xrng64() & 0x00ff3ff7) | 0x7d000000;
		Q[Qoff + 16] = (xrng64() & 0x4ffdffff) | 0x20000000 | (~Q[Qoff + 15] & 0x00020000);
	    
		MD5_REVERSE_STEP(5, 0x4787c62a, 12);
		MD5_REVERSE_STEP(6, 0xa8304613, 17);
		MD5_REVERSE_STEP(7, 0xfd469501, 22);
		MD5_REVERSE_STEP(11, 0x895cd7be, 22);
		MD5_REVERSE_STEP(14, 0xa679438e, 17);
		MD5_REVERSE_STEP(15, 0x49b40821, 22);

		const uint32 tt17 = GG(Q[Qoff + 16], Q[Qoff + 15], Q[Qoff + 14]) + Q[Qoff + 13] + 0xf61e2562;
		const uint32 tt18 = Q[Qoff + 14] + 0xc040b340 + block[6];
		const uint32 tt19 = Q[Qoff + 15] + 0x265e5a51 + block[11];

		const uint32 tt0 = FF(Q[Qoff + 0], Q[Qoff - 1], Q[Qoff - 2]) + Q[Qoff - 3] + 0xd76aa478;
		const uint32 tt1 = Q[Qoff - 2] + 0xe8c7b756;		

		const uint32 q1a = 0x04200040 | (Q[Qoff + 2] & 0xf01e1080);
		
		unsigned counter = 0;
		while (counter < (1 << 12))
		{
			++counter;

			uint32 q1 = q1a | (xrng64() & 0x01c0e71f);
			uint32 m1 = Q[Qoff+2] - q1;
			m1 = RR(m1, 12) - FF(q1, Q[Qoff+0], Q[Qoff-1]) - tt1;

			const uint32 q16 = Q[Qoff+16];
			uint32 q17 = tt17 + m1;
			q17 = RL(q17, 5) + q16;
			if (0x40000000 != ((q17^q16) & 0xc0008008)) continue;
			if (0 != (q17 & 0x00020000)) continue;

			uint32 q18 = GG(q17, q16, Q[Qoff+15]) + tt18;
			q18 = RL(q18, 9); q18 += q17;
			if (0x00020000 != ((q18^q17) & 0xa0020000)) continue;

			uint32 q19 = GG(q18, q17, q16) + tt19;
			q19 = RL(q19, 14); q19 += q18;
			if (0 != (q19 & 0x80020000)) continue;

			uint32 m0 = q1 - Q[Qoff + 0];
			m0 = RR(m0, 7) - tt0;

			uint32 q20 = GG(q19, q18, q17) + q16 + 0xe9b6c7aa + m0;
			q20 = RL(q20, 20); q20 += q19;
			if (0x00040000 != ((q20^q19) & 0x80040000))	continue;
			
			Q[Qoff + 1] = q1;
			Q[Qoff + 17] = q17;
			Q[Qoff + 18] = q18;
			Q[Qoff + 19] = q19;
			Q[Qoff + 20] = q20;

			block[0] = m0;
			block[1] = m1;
			MD5_REVERSE_STEP(2, 0x242070db, 17);

			counter = 0;
			break;
		}
		if (counter != 0)
			continue;

		const uint32 q4b = Q[Qoff + 4];
		const uint32 q9b = Q[Qoff + 9];
		const uint32 q10b = Q[Qoff + 10];
		const uint32 tt21 = GG(Q[Qoff+20], Q[Qoff+19], Q[Qoff+18]) + Q[Qoff+17] + 0xd62f105d;

		counter = 0;
		while (counter < (1<<6))
		{
			Q[Qoff + 4] = q4b ^ q4mask[counter];
			++counter;
			MD5_REVERSE_STEP(5, 0x4787c62a, 12);
			uint32 q21 = tt21 + block[5];
			q21 = RL(q21, 5); q21 += Q[Qoff+20];
			if (0 != ((q21^Q[Qoff+20]) & 0x80020000)) continue;

			Q[Qoff+21] = q21;
			MD5_REVERSE_STEP(3, 0xc1bdceee, 22);
			MD5_REVERSE_STEP(4, 0xf57c0faf, 7);
			MD5_REVERSE_STEP(7, 0xfd469501, 22);

			const uint32 tt10 = Q[Qoff + 7] + 0xffff5bb1;
			const uint32 tt22 = GG(Q[Qoff + 21], Q[Qoff + 20], Q[Qoff + 19]) + Q[Qoff + 18] + 0x02441453;
			const uint32 tt23 = Q[Qoff + 19] + 0xd8a1e681 + block[15];
			const uint32 tt24 = Q[Qoff + 20] + 0xe7d3fbc8 + block[4];
	 
			unsigned counter2 = 0;
			while (counter2 < (1<<5))
			{
				uint32 q10 = q10b ^ wang_q_masks.q10[counter2];
				uint32 m10 = RR(Q[Qoff+11]-q10,17);
				uint32 q9 = q9b ^ wang_q_masks.q9[counter2];
				++counter2;

				m10 -= FF(q10, q9, Q[Qoff+8]) + tt10;

				uint32 aa = Q[Qoff + 21];
				uint32 dd = tt22+m10; dd = RL(dd, 9) + aa;
				if (0 != (dd & 0x80000000)) continue;			

				uint32 bb = Q[Qoff + 20];
				uint32 cc = tt23 + GG(dd, aa, bb); 
				if (0 != (cc & 0x20000)) continue;
				cc = RL(cc, 14) + dd;
				if (0 != (cc & 0x80000000)) continue;

				bb = tt24 + GG(cc, dd, aa); bb = RL(bb, 20) + cc;
				if (0 == (bb & 0x80000000)) continue;

				block[10] = m10;
				Q[Qoff + 9] = q9;
				Q[Qoff + 10] = q10;
				MD5_REVERSE_STEP(13, 0xfd987193, 12);


#if defined(__AVX512F__)
                const uint32 tt8v_s = FF(Q[Qoff+8], Q[Qoff+7], Q[Qoff+6]) + Q[Qoff+5] + 0x698098d8u;
                const uint32 tt9v_s = Q[Qoff+6] + 0x8b44f7afu;
                const uint32 tt12v_s = RR(Q[Qoff+13]-Q[Qoff+12],7)
                    - FF(Q[Qoff+12],Q[Qoff+11],q10) - 0x6b901122u;
                const __m512i vq9base=_mm512_set1_epi32((int)q9), vq8=_mm512_set1_epi32((int)Q[Qoff+8]), vq7=_mm512_set1_epi32((int)Q[Qoff+7]), vq10=_mm512_set1_epi32((int)q10);
                const __m512i vtt8=_mm512_set1_epi32((int)tt8v_s), vtt9=_mm512_set1_epi32((int)tt9v_s), vtt12=_mm512_set1_epi32((int)tt12v_s);
                const __m512i vaa0=_mm512_set1_epi32((int)aa), vbb0=_mm512_set1_epi32((int)bb), vcc0=_mm512_set1_epi32((int)cc), vdd0=_mm512_set1_epi32((int)dd);
                const __m512i signbit=_mm512_set1_epi32((int)0x80000000u), bit15=_mm512_set1_epi32(1<<15);
                alignas(64) uint32 la[16],lb[16],lc[16],ld[16],l8[16],l9[16],l12[16];
#define WG(x,y,z) _mm512_xor_si512((y),_mm512_and_si512((z),_mm512_xor_si512((x),(y))))
#define WH(x,y,z) _mm512_xor_si512(_mm512_xor_si512((x),(y)),(z))
#define WI(x,y,z) _mm512_xor_si512((y),_mm512_or_si512((x),_mm512_xor_si512((z),_mm512_set1_epi32(-1))))
#define WSG(A,B,C,D,M,K,R) do{(A)=_mm512_add_epi32((A),_mm512_add_epi32(WG((B),(C),(D)),_mm512_add_epi32((M),_mm512_set1_epi32((int)(K)))));(A)=_mm512_add_epi32(_mm512_rol_epi32((A),(R)),(B));}while(0)
#define WSH(A,B,C,D,M,K,R) do{(A)=_mm512_add_epi32((A),_mm512_add_epi32(WH((B),(C),(D)),_mm512_add_epi32((M),_mm512_set1_epi32((int)(K)))));(A)=_mm512_add_epi32(_mm512_rol_epi32((A),(R)),(B));}while(0)
#define WSI(A,B,C,D,M,K,R) do{(A)=_mm512_add_epi32((A),_mm512_add_epi32(WI((B),(C),(D)),_mm512_add_epi32((M),_mm512_set1_epi32((int)(K)))));(A)=_mm512_add_epi32(_mm512_rol_epi32((A),(R)),(B));}while(0)
                for (unsigned k9 = 0; k9 < q9mask2.size(); k9 += 16) {
                    if ((k9 & 255u) == 0 && miniclash_cancelled()) return;
                    __m512i q9v=_mm512_xor_si512(vq9base,_mm512_loadu_si512((const void*)(&q9mask2[k9])));
                    __m512i b12v=_mm512_sub_epi32(vtt12,q9v);
                    __m512i b8v=_mm512_sub_epi32(_mm512_ror_epi32(_mm512_sub_epi32(q9v,vq8),7),vtt8);
                    __m512i ff9=_mm512_xor_si512(vq7,_mm512_and_si512(q9v,_mm512_xor_si512(vq8,vq7)));
                    __m512i b9v=_mm512_sub_epi32(_mm512_sub_epi32(_mm512_ror_epi32(_mm512_sub_epi32(vq10,q9v),12),ff9),vtt9);
                    __m512i a=vaa0,b=vbb0,c=vcc0,d=vdd0;
                    WSG(a,b,c,d,b9v,0x21e1cde6,5); WSG(d,a,b,c,_mm512_set1_epi32((int)block[14]),0xc33707d6,9);
                    WSG(c,d,a,b,_mm512_set1_epi32((int)block[3]),0xf4d50d87,14); WSG(b,c,d,a,b8v,0x455a14ed,20);
                    WSG(a,b,c,d,_mm512_set1_epi32((int)block[13]),0xa9e3e905,5); WSG(d,a,b,c,_mm512_set1_epi32((int)block[2]),0xfcefa3f8,9);
                    WSG(c,d,a,b,_mm512_set1_epi32((int)block[7]),0x676f02d9,14); WSG(b,c,d,a,b12v,0x8d2a4c8a,20);
                    WSH(a,b,c,d,_mm512_set1_epi32((int)block[5]),0xfffa3942,4); WSH(d,a,b,c,b8v,0x8771f681,11);
                    c=_mm512_add_epi32(c,_mm512_add_epi32(WH(d,a,b),_mm512_set1_epi32((int)(block[11]+0x6d9d6122u))));
                    __mmask16 active=_mm512_test_epi32_mask(c,bit15); if(!active) continue;
                    c=_mm512_add_epi32(_mm512_rol_epi32(c,16),d);
                    WSH(b,c,d,a,_mm512_set1_epi32((int)block[14]),0xfde5380c,23); WSH(a,b,c,d,_mm512_set1_epi32((int)block[1]),0xa4beea44,4);
                    WSH(d,a,b,c,_mm512_set1_epi32((int)block[4]),0x4bdecfa9,11); WSH(c,d,a,b,_mm512_set1_epi32((int)block[7]),0xf6bb4b60,16);
                    WSH(b,c,d,a,_mm512_set1_epi32((int)block[10]),0xbebfbc70,23); WSH(a,b,c,d,_mm512_set1_epi32((int)block[13]),0x289b7ec6,4);
                    WSH(d,a,b,c,_mm512_set1_epi32((int)block[0]),0xeaa127fa,11); WSH(c,d,a,b,_mm512_set1_epi32((int)block[3]),0xd4ef3085,16);
                    WSH(b,c,d,a,_mm512_set1_epi32((int)block[6]),0x04881d05,23); WSH(a,b,c,d,b9v,0xd9d4d039,4);
                    WSH(d,a,b,c,b12v,0xe6db99e5,11); WSH(c,d,a,b,_mm512_set1_epi32((int)block[15]),0x1fa27cf8,16);
                    WSH(b,c,d,a,_mm512_set1_epi32((int)block[2]),0xc4ac5665,23); active&=_mm512_testn_epi32_mask(_mm512_xor_si512(b,d),signbit); if(!active) continue;
                    WSI(a,b,c,d,_mm512_set1_epi32((int)block[0]),0xf4292244,6); active&=_mm512_testn_epi32_mask(_mm512_xor_si512(a,c),signbit); if(!active) continue;
                    WSI(d,a,b,c,_mm512_set1_epi32((int)block[7]),0x432aff97,10); active&=_mm512_test_epi32_mask(_mm512_xor_si512(b,d),signbit); if(!active) continue;
                    WSI(c,d,a,b,_mm512_set1_epi32((int)block[14]),0xab9423a7,15); active&=_mm512_testn_epi32_mask(_mm512_xor_si512(a,c),signbit); if(!active) continue;
                    WSI(b,c,d,a,_mm512_set1_epi32((int)block[5]),0xfc93a039,21); active&=_mm512_testn_epi32_mask(_mm512_xor_si512(b,d),signbit); if(!active) continue;
                    _mm512_store_si512((void*)la,a);_mm512_store_si512((void*)lb,b);_mm512_store_si512((void*)lc,c);_mm512_store_si512((void*)ld,d);
                    _mm512_store_si512((void*)l8,b8v);_mm512_store_si512((void*)l9,b9v);_mm512_store_si512((void*)l12,b12v);
                    unsigned sm=(unsigned)active; while(sm){unsigned lane=__builtin_ctz(sm);sm&=sm-1;block[8]=l8[lane];block[9]=l9[lane];block[12]=l12[lane];if(finish_wang_candidate(block,IV,la[lane],lb[lane],lc[lane],ld[lane]))return;}
                }
#undef WSI
#undef WSH
#undef WSG
#undef WI
#undef WH
#undef WG
#else
				for (unsigned k9 = 0; k9 < (1<<10);)
				{
					uint32 a = aa, b = bb, c = cc, d = dd;
					Q[Qoff + 9] = q9 ^ q9mask2[k9]; ++k9;
					MD5_REVERSE_STEP(8, 0x698098d8, 7);
					MD5_REVERSE_STEP(9, 0x8b44f7af, 12);
					MD5_REVERSE_STEP(12, 0x6b901122, 7);

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
					if (0 == (c & (1 << 15))) 
						continue;
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
					if (0 != ((b^d) & 0x80000000))
						continue;

					MD5_STEP(II, a, b, c, d, block[0], 0xf4292244, 6);
					if (0 != ((a^c) >> 31)) continue;
					MD5_STEP(II, d, a, b, c, block[7], 0x432aff97, 10);
					if (0 == ((b^d) >> 31)) continue;
					MD5_STEP(II, c, d, a, b, block[14], 0xab9423a7, 15);
					if (0 != ((a^c) >> 31)) continue;
					MD5_STEP(II, b, c, d, a, block[5], 0xfc93a039, 21);
					if (0 != ((b^d) >> 31)) continue;
					MD5_STEP(II, a, b, c, d, block[12], 0x655b59c3, 6);
					if (0 != ((a^c) >> 31)) continue;
					MD5_STEP(II, d, a, b, c, block[3], 0x8f0ccc92, 10);
					if (0 != ((b^d) >> 31)) continue;
					MD5_STEP(II, c, d, a, b, block[10], 0xffeff47d, 15);
					if (0 != ((a^c) >> 31)) continue;
					MD5_STEP(II, b, c, d, a, block[1], 0x85845dd1, 21);
					if (0 != ((b^d) >> 31)) continue;
					MD5_STEP(II, a, b, c, d, block[8], 0x6fa87e4f, 6);
					if (0 != ((a^c) >> 31)) continue;
					MD5_STEP(II, d, a, b, c, block[15], 0xfe2ce6e0, 10);
					if (0 != ((b^d) >> 31)) continue;
					MD5_STEP(II, c, d, a, b, block[6], 0xa3014314, 15);
					if (0 != ((a^c) >> 31)) continue;
					MD5_STEP(II, b, c, d, a, block[13], 0x4e0811a1, 21);
					if (0 == ((b^d) >> 31)) continue;
					MD5_STEP(II, a, b, c, d, block[4], 0xf7537e82, 6);
					if (0 != ((a^c) >> 31)) continue;
					MD5_STEP(II, d, a, b, c, block[11], 0xbd3af235, 10);
					if (0 != ((b^d) >> 31)) continue;
					MD5_STEP(II, c, d, a, b, block[2], 0x2ad7d2bb, 15);
					if (0 != ((a^c) >> 31)) continue;
					MD5_STEP(II, b, c, d, a, block[9], 0xeb86d391, 21);


					uint32 block2[16];
					uint32 IV1[4], IV2[4];
					for (int t = 0; t < 4; ++t)
					{
						IV1[t] = IV[t];
						IV2[t] = IV[t] + (1 << 31);
					}
					IV2[1] += (1 << 25);
					IV2[2] += (1 << 25);
					IV2[3] += (1 << 25);

					for (int t = 0; t < 16; ++t)
						block2[t] = block[t];
					block2[4] += 1<<31;
					block2[11] -= 1<<15;
					block2[14] += 1<<31;

					md5_compress(IV1, block);
					md5_compress(IV2, block2);
					if (IV2[0]==IV1[0] && IV2[1]==IV1[1] && IV2[2]==IV1[2] && IV2[3]==IV1[3])
						return;

					if (IV2[0] != IV1[0]) { }
				}
#endif

			}
		}
	}
}
