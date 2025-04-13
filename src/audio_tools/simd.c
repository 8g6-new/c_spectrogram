__m512 atan2_ps(__m512 y, __m512 x) {
    __m512 zero     = _mm512_set1_ps(0.0f);
    __m512 one      = _mm512_set1_ps(1.0f);
    __m512 pi       = _mm512_set1_ps(3.14159265358979323846f);
    __m512 pi_2     = _mm512_set1_ps(1.57079632679489661923f);  // pi/2
    __m512 pi_4     = _mm512_set1_ps(0.785398163397448309616f); // pi/4
    __m512 c        = _mm512_set1_ps(0.273f);                   // correction factor

    __mmask16 x_pos = _mm512_cmp_ps_mask(x, zero, _CMP_GT_OS);
    __mmask16 x_neg = _mm512_cmp_ps_mask(x, zero, _CMP_LT_OS);
    __mmask16 y_pos = _mm512_cmp_ps_mask(y, zero, _CMP_GT_OS);
    __mmask16 y_neg = _mm512_cmp_ps_mask(y, zero, _CMP_LT_OS);
    __mmask16 x_zero = _mm512_cmp_ps_mask(x, zero, _CMP_EQ_OS);

    __m512 z = _mm512_div_ps(y, x);
    __m512 abs_z = _mm512_abs_ps(z);

    // atan approximation
    __m512 atan = _mm512_fmadd_ps(
        c, _mm512_mul_ps(_mm512_sub_ps(abs_z, one), z),
        _mm512_mul_ps(z, pi_4)
    );

    // quadrant correction
    __m512 result = atan;

    // x < 0 and y >= 0: atan + pi
    result = _mm512_mask_add_ps(result, x_neg & ~y_neg, atan, pi);

    // x < 0 and y < 0: atan - pi
    result = _mm512_mask_sub_ps(result, x_neg & y_neg, atan, pi);

    // x == 0 and y > 0: pi/2
    result = _mm512_mask_mov_ps(result, x_zero & y_pos, pi_2);

    // x == 0 and y < 0: -pi/2
    result = _mm512_mask_mov_ps(result, x_zero & y_neg, _mm512_sub_ps(zero, pi_2));

    return result;
}

