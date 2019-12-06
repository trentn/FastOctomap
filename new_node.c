static inline __m256i new_node0(__m256d txmv, __m256d tymv, __m256d tzmv)
{
    __m256d tmp1v = _mm256_sub_pd(txmv, tymv);
    __m256d tmp2v = _mm256_sub_pd(txmv, tzmv);
    __m256d tmp3v = _mm256_sub_pd(tymv, tzmv);

    __m256d signbitv = _mm256_broadcast_sd((double*)&signbit);

    //X
    //currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & 0x8000000000000000)>>;
    __m256d tmpv = _mm256_and_pd(tmp1v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    __m256i currentNodev = _mm256_srli_epi64(_mm256_castpd_si256(tmpv),59);

    
    //Y
    //currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>;
    tmpv = _mm256_andnot_pd(tmp1v,tmp3v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),61));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    
    //Z
    //currentNode |= (unsigned long long)( ~*(unsigned long long*)(&tmp3) & ~*(unsigned long long*)(&tmp2) & 0x8000000000000000 )>>;
    __m256i onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&all_ones));
    tmp2v = _mm256_castsi256_pd(_mm256_xor_si256(_mm256_castpd_si256(tmp2v),onesv));
    tmpv = _mm256_andnot_pd(tmp3v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),2));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    return currentNodev;
}

static inline __m256i new_node1(__m256d txmv, __m256d tymv, __m256d tzmv)
{
    /*
    long x_shift = (63-x);
    long y_shift = (63-y);
    long z_shift = (63-z);
    */

    __m256d tmp1v = _mm256_sub_pd(txmv, tymv);
    __m256d tmp2v = _mm256_sub_pd(txmv, tzmv);
    __m256d tmp3v = _mm256_sub_pd(tymv, tzmv);

    __m256d signbitv = _mm256_broadcast_sd((double*)&signbit);

    //X
    //currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & 0x8000000000000000)>>;
    __m256d tmpv = _mm256_and_pd(tmp1v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    __m256i currentNodev = _mm256_srli_epi64(_mm256_castpd_si256(tmpv),58);

    
    //Y
    //currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>;
    tmpv = _mm256_andnot_pd(tmp1v,tmp3v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),60));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    
    //Z
    //currentNode |= (unsigned long long)( ~*(unsigned long long*)(&tmp3) & ~*(unsigned long long*)(&tmp2) & 0x8000000000000000 )>>;
    __m256i onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&all_ones));
    tmp2v = _mm256_castsi256_pd(_mm256_xor_si256(_mm256_castpd_si256(tmp2v),onesv));
    tmpv = _mm256_andnot_pd(tmp3v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),55));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    return currentNodev;
}

static inline __m256i new_node2(__m256d txmv, __m256d tymv, __m256d tzmv)
{
    /*
    long x_shift = (63-x);
    long y_shift = (63-y);
    long z_shift = (63-z);
    */

    __m256d tmp1v = _mm256_sub_pd(txmv, tymv);
    __m256d tmp2v = _mm256_sub_pd(txmv, tzmv);
    __m256d tmp3v = _mm256_sub_pd(tymv, tzmv);

    __m256d signbitv = _mm256_broadcast_sd((double*)&signbit);

    //X
    //currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & 0x8000000000000000)>>;
    __m256d tmpv = _mm256_and_pd(tmp1v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    __m256i currentNodev = _mm256_srli_epi64(_mm256_castpd_si256(tmpv),57);

    
    //Y
    //currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>;
    tmpv = _mm256_andnot_pd(tmp1v,tmp3v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),55));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    
    //Z
    //currentNode |= (unsigned long long)( ~*(unsigned long long*)(&tmp3) & ~*(unsigned long long*)(&tmp2) & 0x8000000000000000 )>>;
    __m256i onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&all_ones));
    tmp2v = _mm256_castsi256_pd(_mm256_xor_si256(_mm256_castpd_si256(tmp2v),onesv));
    tmpv = _mm256_andnot_pd(tmp3v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),60));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    return currentNodev;
}

static inline __m256i new_node3(__m256d txmv, __m256d tymv, __m256d tzmv)
{
    /*
    long x_shift = (63-x);
    long y_shift = (63-y);
    long z_shift = (63-z);
    */

    __m256d tmp1v = _mm256_sub_pd(txmv, tymv);
    __m256d tmp2v = _mm256_sub_pd(txmv, tzmv);
    __m256d tmp3v = _mm256_sub_pd(tymv, tzmv);

    __m256d signbitv = _mm256_broadcast_sd((double*)&signbit);

    //X
    //currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & 0x8000000000000000)>>;
    __m256d tmpv = _mm256_and_pd(tmp1v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    __m256i currentNodev = _mm256_srli_epi64(_mm256_castpd_si256(tmpv),56);

    
    //Y
    //currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>;
    tmpv = _mm256_andnot_pd(tmp1v,tmp3v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),55));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    
    //Z
    //currentNode |= (unsigned long long)( ~*(unsigned long long*)(&tmp3) & ~*(unsigned long long*)(&tmp2) & 0x8000000000000000 )>>;
    __m256i onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&all_ones));
    tmp2v = _mm256_castsi256_pd(_mm256_xor_si256(_mm256_castpd_si256(tmp2v),onesv));
    tmpv = _mm256_andnot_pd(tmp3v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),55));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    return currentNodev;
}

static inline __m256i new_node4(__m256d txmv, __m256d tymv, __m256d tzmv)
{
    /*
    long x_shift = (63-x);
    long y_shift = (63-y);
    long z_shift = (63-z);
    */

    __m256d tmp1v = _mm256_sub_pd(txmv, tymv);
    __m256d tmp2v = _mm256_sub_pd(txmv, tzmv);
    __m256d tmp3v = _mm256_sub_pd(tymv, tzmv);

    __m256d signbitv = _mm256_broadcast_sd((double*)&signbit);

    //X
    //currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & 0x8000000000000000)>>;
    __m256d tmpv = _mm256_and_pd(tmp1v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    __m256i currentNodev = _mm256_srli_epi64(_mm256_castpd_si256(tmpv),55);

    
    //Y
    //currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>;
    tmpv = _mm256_andnot_pd(tmp1v,tmp3v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),57));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    
    //Z
    //currentNode |= (unsigned long long)( ~*(unsigned long long*)(&tmp3) & ~*(unsigned long long*)(&tmp2) & 0x8000000000000000 )>>;
    __m256i onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&all_ones));
    tmp2v = _mm256_castsi256_pd(_mm256_xor_si256(_mm256_castpd_si256(tmp2v),onesv));
    tmpv = _mm256_andnot_pd(tmp3v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),58));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    return currentNodev;
}

static inline __m256i new_node5(__m256d txmv, __m256d tymv, __m256d tzmv)
{
    /*
    long x_shift = (63-x);
    long y_shift = (63-y);
    long z_shift = (63-z);
    */

    __m256d tmp1v = _mm256_sub_pd(txmv, tymv);
    __m256d tmp2v = _mm256_sub_pd(txmv, tzmv);
    __m256d tmp3v = _mm256_sub_pd(tymv, tzmv);

    __m256d signbitv = _mm256_broadcast_sd((double*)&signbit);

    //X
    //currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & 0x8000000000000000)>>;
    __m256d tmpv = _mm256_and_pd(tmp1v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    __m256i currentNodev = _mm256_srli_epi64(_mm256_castpd_si256(tmpv),55);

    
    //Y
    //currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>;
    tmpv = _mm256_andnot_pd(tmp1v,tmp3v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),56));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    
    //Z
    //currentNode |= (unsigned long long)( ~*(unsigned long long*)(&tmp3) & ~*(unsigned long long*)(&tmp2) & 0x8000000000000000 )>>;
    __m256i onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&all_ones));
    tmp2v = _mm256_castsi256_pd(_mm256_xor_si256(_mm256_castpd_si256(tmp2v),onesv));
    tmpv = _mm256_andnot_pd(tmp3v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),55));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    return currentNodev;
}

static inline __m256i new_node6(__m256d txmv, __m256d tymv, __m256d tzmv)
{
    /*
    long x_shift = (63-x);
    long y_shift = (63-y);
    long z_shift = (63-z);
    */

    __m256d tmp1v = _mm256_sub_pd(txmv, tymv);
    __m256d tmp2v = _mm256_sub_pd(txmv, tzmv);
    __m256d tmp3v = _mm256_sub_pd(tymv, tzmv);

    __m256d signbitv = _mm256_broadcast_sd((double*)&signbit);

    //X
    //currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & 0x8000000000000000)>>;
    __m256d tmpv = _mm256_and_pd(tmp1v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    __m256i currentNodev = _mm256_srli_epi64(_mm256_castpd_si256(tmpv),55);

    
    //Y
    //currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>;
    tmpv = _mm256_andnot_pd(tmp1v,tmp3v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),55));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    
    //Z
    //currentNode |= (unsigned long long)( ~*(unsigned long long*)(&tmp3) & ~*(unsigned long long*)(&tmp2) & 0x8000000000000000 )>>;
    __m256i onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&all_ones));
    tmp2v = _mm256_castsi256_pd(_mm256_xor_si256(_mm256_castpd_si256(tmp2v),onesv));
    tmpv = _mm256_andnot_pd(tmp3v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmpv),56));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    return currentNodev;
}