#include <Rasterer2d.h>
#include <Tree2PS.h>

#include "gtest/gtest.h"
#include <iostream>
#include <math.h>

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

bool
equivalent(double x, double y, double tol = 1.0e-05 )
{
        return fabs(x-y) < tol;
}

TEST(HaarLinearTest, Quad0) {

	float cval[3];
	double px, py;
	double qx, qy;
	float v[4];

	cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
    v[0] = 0.4; v[1] = 0.3;
    v[2] = 0.2; v[3] = 0.4 ;
    calc_coeffs_line(cval, v, 0, 0 );

    px = 0.2; py = 0.15;
    qx = 0.1; qy = 0.2;
	Rasterer2d rast(1);
	rast.rasterize(px, py, qx, qy);

	EXPECT_FLOAT_EQ(rast.m_root->m_haar[0], cval[0]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[1], cval[1]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[2], cval[2]);
#if 0
    std::cout << "Haar0 (T1:Q0):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
    std::cout << "Haar1 (T1:Q0):" << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
    std::cout << "Haar2 (T1:Q0):" << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
}

TEST(HaarLinearTest, Quad1) {

	float cval[3];
	double px, py;
	double qx, qy;
	float v[4];

	cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
    v[0] = 0.4; v[1] = 0.3;
    v[2] = 0.2; v[3] = 0.4 ;
    calc_coeffs_line(cval, v, 1, 0 );

    px = 0.7; py = 0.15;
    qx = 0.6; qy = 0.2;
    Rasterer2d rast(1);
    rast.rasterize(px, py, qx, qy);
    assert(rast.m_root);

	EXPECT_FLOAT_EQ(rast.m_root->m_haar[0], cval[0]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[1], cval[1]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[2], cval[2]);
#if 1
    std::cout << "Haar0 (T1:Q1):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
    std::cout << "Haar1 (T1:Q1):" << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
    std::cout << "Haar2 (T1:Q1):" << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
}

TEST(HaarLinearTest, Quad2) {

	float cval[3];
	double px, py;
	double qx, qy;
	float v[4];

    cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
    v[0] = 0.4; v[1] = 0.3;
    v[2] = 0.2; v[3] = 0.4 ;
    calc_coeffs_line(cval, v, 0, 1 );

    px = 0.2; py = 0.65;
    qx = 0.1; qy = 0.7;
    Rasterer2d rast(1);
    rast.rasterize(px, py, qx, qy);
    assert(rast.m_root);

	EXPECT_FLOAT_EQ(rast.m_root->m_haar[0], cval[0]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[1], cval[1]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[2], cval[2]);

#if 1
    std::cout << "Haar0 (T1:Q1):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
    std::cout << "Haar1 (T1:Q1):" << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
    std::cout << "Haar2 (T1:Q1):" << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
}

TEST(HaarLinearTest, Quad3) {

	float cval[3];
	double px, py;
	double qx, qy;
	float v[4];

	cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
    v[0] = 0.4; v[1] = 0.3;
    v[2] = 0.2; v[3] = 0.4 ;
    calc_coeffs_line(cval, v, 1, 1 );

    px = 0.7; py = 0.65;
    qx = 0.6; qy = 0.7;
    Rasterer2d rast(1);
    rast.rasterize(px, py, qx, qy);
    assert(rast.m_root);

	EXPECT_FLOAT_EQ(rast.m_root->m_haar[0], cval[0]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[1], cval[1]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[2], cval[2]);

#if 1
    std::cout << "Haar0 (T1:Q1):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
    std::cout << "Haar1 (T1:Q1):" << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
    std::cout << "Haar2 (T1:Q1):" << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
}


TEST(HaarLinearTest, Quad01) {

	float cval[3];
	double px, py;
	double qx, qy;
	float v[4];

	// Q0 - Q1
    cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
    v[0] = 0.0; v[1] = 0.0;
    v[2] = 1.0; v[3] = 0.5 ;
    calc_coeffs_line(cval, v, 0, 0 );
    v[0] = 0.0; v[1] = 0.5;
    v[2] = 0.5; v[3] = 0.75 ;
    calc_coeffs_line(cval, v, 1, 0 );

    px = 0.0; py = 0.0;
    qx = 0.75; qy = 0.375;
    Rasterer2d rast(1);
    rast.rasterize(px, py, qx, qy);

	EXPECT_FLOAT_EQ(rast.m_root->m_haar[0], cval[0]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[1], cval[1]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[2], cval[2]);
#if 0
    cout << "Haar0 (T1:Q1):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
    cout << "Haar1 (T1:Q1):" << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
    cout << "Haar2 (T1:Q1):" << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
}

TEST(HaarLinearTest, Quad13) {

	float cval[3];
	double px, py;
	double qx, qy;
	float v[4];

	// Q1 - Q3
    cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
    v[0] = 0.0; v[1] = 0.0;
    v[2] = 0.5; v[3] = 1.0 ;
    calc_coeffs_line(cval, v, 1, 0 );
    v[0] = 0.5; v[1] = 0.0;
    v[2] = 0.8; v[3] = 0.6 ;
    calc_coeffs_line(cval, v, 1, 1 );

    px = 0.5; py = 0.0;
    qx = 0.9; qy = 0.8;
    Rasterer2d rast(1);
    rast.rasterize(px, py, qx, qy);
    assert(rast.m_root);

	EXPECT_FLOAT_EQ(rast.m_root->m_haar[0], cval[0]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[1], cval[1]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[2], cval[2]);
}

TEST(HaarLinearTest, Quad23) {

	float cval[3];
	double px, py;
	double qx, qy;
	float v[4];

    cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
    v[0] = 0.0; v[1] = 0.0;
    v[2] = 1.0; v[3] = 0.5 ;
    calc_coeffs_line(cval, v, 0, 1 );
    v[0] = 0.0; v[1] = 0.5;
    v[2] = 0.5; v[3] = 0.75 ;
    calc_coeffs_line(cval, v, 1, 1 );

    px = 0.0; py = 0.5;
    qx = 0.75; qy = 0.875;
    Rasterer2d rast(1);
    rast.rasterize(px, py, qx, qy);
    assert(rast.m_root);

	EXPECT_FLOAT_EQ(rast.m_root->m_haar[0], cval[0]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[1], cval[1]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[2], cval[2]);
}

TEST(HaarLinearTest, Quad02) {

	float cval[3];
	double px, py;
	double qx, qy;
	float v[4];

    cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
    v[0] = 0.0; v[1] = 0.0;
    v[2] = 0.5; v[3] = 1.0 ;
    calc_coeffs_line(cval, v, 0, 0 );
    v[0] = 0.5; v[1] = 0.0;
    v[2] = 0.8; v[3] = 0.6 ;
    calc_coeffs_line(cval, v, 0, 1 );

    px = 0.0; py = 0.0;
    qx = 0.4; qy = 0.8;
    Rasterer2d rast(1);
    rast.rasterize(px, py, qx, qy);
    assert(rast.m_root);

	EXPECT_FLOAT_EQ(rast.m_root->m_haar[0], cval[0]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[1], cval[1]);
	EXPECT_FLOAT_EQ(rast.m_root->m_haar[2], cval[2]);
}

TEST(RasterizeTest, square) {

	float x[4];
	float y[4];
	
	x[0] = y[0] = 0.;
	x[1] = 1.0; y[1] = 0.;
	x[2] = 1.0; y[2] = 1.;
	x[3] = 0.0; y[3] = 1.;

	int depth = 2;
	Rasterer2d rast(depth);
	for( int i = 0; i < 4; i++ ) {
   		rast.rasterize(x[i], y[i], x[(i+1)%4], y[(i+1)%4]);
	}

	std::vector<std::vector<double>> img;
	rast.render(img);
	EXPECT_FLOAT_EQ(img.size(),exp2(depth));
	for( int i = 0; i < img.size(); i++ ) {
		EXPECT_FLOAT_EQ(img[i].size(), exp2(depth));
		for( int j = 0; j <  img[i].size(); j++ )
			EXPECT_FLOAT_EQ(img[i][j], 1.0);
	}

	rast.toPostscript(img, "square");

	{
		std::string name("squareQ.ps");
		Tree2PS t2ps;
		t2ps.write(rast,name);
	}
}

#if 0
TEST(RasterizeTest, DISABLE_Letter_L) {

	float x[6];
	float y[6];
	
	x[0] = y[0] = 0.;
	x[1] = 1.0; y[1] = 0.;
	x[2] = 1.0; y[2] = 0.35;
	x[3] = 0.45; y[3] = 0.25;
	x[4] = 0.25; y[4] = 1.0;
	x[5] = 0.0; y[5] = 1.0;

	int depth = 5;
	Rasterer2d rast(depth);
	for( int i = 0; i < 6; i++ ) {
   		rast.rasterize(x[i], y[i], x[(i+1)%6], y[(i+1)%6]);
	}

	std::vector<std::vector<double>> img;
	rast.render(img);
	EXPECT_FLOAT_EQ(img.size(),exp2(depth));
	for( int i = 0; i < img.size(); i++ ) {
		EXPECT_FLOAT_EQ(img[i].size(), exp2(depth));
		for( int j = 0; j < img[i].size(); j++ ) {
			if( i < 2 ) EXPECT_FLOAT_EQ(img[i][j], 1.0);
		}
	}
	rast.toPostscript(img, "letterL");
	rast.treeToPostscript("letterLQ");
}
#endif

TEST(RasterizeTest, Triangle) {

	float x[3];
	float y[3];
	
	x[0] = 0.75; y[0] = 1.;
	x[1] = 0.25; y[1] = 1.;
	x[2] = 0.5; y[2] = 0.0;

	int depth = 4;
	Rasterer2d rast(depth);
	for( int i = 0; i < 3; i++ ) {
   		rast.rasterize(x[i], y[i], x[(i+1)%3], y[(i+1)%3]);
	}

	std::vector<std::vector<double>> img;
	rast.render(img);
	EXPECT_FLOAT_EQ(img.size(),exp2(depth));
	for( int i = 0; i < img.size(); i++ ) {
		EXPECT_FLOAT_EQ(img[i].size(), exp2(depth));
	}
	rast.toPostscript(img, "triangle.ps");
	rast.treeToPostscript("triangleQ.ps");
}

#if 0
TEST(RasterizeTest, DISABLE_Hexagon) {

	float x[6];
	float y[6];
	
	x[0] = 0.25; y[0] = 0.;
	x[1] = 0.75; y[1] = 0.;
	x[2] = 1.0; y[2] = 0.5;
	x[3] = 0.75; y[3] = 1.0;
	x[4] = 0.25; y[4] = 1.0;
	x[5] = 0.0; y[5] = 0.5;

	int depth = 4;
	Rasterer2d rast(depth);
	for( int i = 0; i < 6; i++ ) {
   		rast.rasterize(x[i], y[i], x[(i+1)%6], y[(i+1)%6]);
	}

	std::vector<std::vector<double>> img;
	rast.render(img);
	EXPECT_FLOAT_EQ(img.size(),exp2(depth));
	for( int i = 0; i < img.size(); i++ ) {
		EXPECT_FLOAT_EQ(img[i].size(), exp2(depth));
	}
	rast.toPostscript(img, "Hexagon");
}
#endif
int main(int argc, char **argv ) {
::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
