#include <Rasterer2d.h>
#include "gtest/gtest.h"

void calc_coeffs00(float *coeffs, float* v)
{
    float norm[2] = {(v[3] - v[1])*.125f, (v[0] - v[2])*.125f};
    float lin[2] = {v[0] + v[2], v[1] + v[3]};

    coeffs[1 - 1] += (lin[0])*norm[0]; // 1,0
    coeffs[3 - 1] += (lin[0])*norm[0]; // 1,1
    coeffs[2 - 1] += (lin[1])*norm[1]; // 0,1
}

void calc_coeffs01(float *coeffs, float* v)
{
    float norm[2] = {(v[3] - v[1])*.125f, (v[0] - v[2])*.125f};
    float lin[2] = {v[0] + v[2], v[1] + v[3]};

    coeffs[1 - 1] += (lin[0])*norm[0]; // 1,0
    coeffs[3 - 1] -= (lin[0])*norm[0]; // 1,1
    coeffs[2 - 1] += (2 - lin[1])*norm[1]; // 0,1
}

void calc_coeffs10(float *coeffs, float* v)
{
        float norm[2] = {(v[3] - v[1])*.125f, (v[0] - v[2])*.125f};
        float lin[2] = {v[0] + v[2], v[1] + v[3]};

        coeffs[1 - 1] += (2 - lin[0])*norm[0]; // 1,0
        coeffs[3 - 1] += (2 - lin[0])*norm[0]; // 1,1
        coeffs[2 - 1] += (lin[1])*norm[1]; // 0,1
}

void calc_coeffs11(float *coeffs, float* v)
{
        float norm[2] = {(v[3] - v[1])*.125f, (v[0] - v[2])*.125f};
        float lin[2] = {v[0] + v[2], v[1] + v[3]};

        coeffs[1 - 1] += (2 - lin[0])*norm[0]; // 1,0
        coeffs[3 - 1] -= (2 - lin[0])*norm[0]; // 1,1
        coeffs[2 - 1] += (2 - lin[1])*norm[1]; // 0,1
}

void calc_coeffs_line(float *coeffs, float* v, int i, int j)
{
        if ( i == 0 )
        {
	        if ( j == 0 )
               calc_coeffs00(coeffs,v);
            else
               calc_coeffs01(coeffs,v);
        }
        else
        {
            if ( j == 0 )
                calc_coeffs10(coeffs,v);
            else
               calc_coeffs11(coeffs,v);
        }
}

int main(int argc, char **argv ) {
::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
