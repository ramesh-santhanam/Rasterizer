#include <Rasterer2d.h>
#include <Tree2PS.h>
#include <WriteTGA.h>

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

void makeStar( float cx, float cy,
	std::vector<double>& px,
	std::vector<double>& py,
	double r1, double r2,
	int n)
{
	px.resize(n);
	py.resize(n);

	float r;
	double ang = 0.0;
	double inc = 2*M_PI/n;

   	for( int i = 0; i < n; i++ ) {
       	if( i % 2 == 0 )
       		r = r1;
       	else
       		r = r2;
       	px[i] = cx + r * cos(ang);
        py[i] = cy + r * sin(ang);
        ang = ang + inc;
	}
	
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

	int depth = 9;
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
	writeTGA(img, "square.tga");
	{
		std::string name("squareQ.ps");
		Tree2PS t2ps;
		t2ps.write(rast,name);
	}
}

TEST(CheckTest, checkL )
{

	float x[6], y[6];

	x[0] = 0.0; y[0] = 0.0;
	x[1] = 1.0; y[1] = 0.0;
	x[2] = 1.0; y[2] = 0.25;
	x[3] = 0.45; y[3] = 0.25;
	x[4] = 0.45; y[4] = 1.0;
	x[5] = 0.00; y[5] = 1.0;

	int depth = 8;
	Rasterer2d rast(depth);
   	rast.rasterize(x[0], y[0], x[1], y[1]);
   	rast.rasterize(x[1], y[1], x[2], y[2]);
   	rast.rasterize(x[2], y[2], x[3], y[3]);
   	rast.rasterize(x[3], y[3], x[4], y[4]);
   	rast.rasterize(x[4], y[4], x[5], y[5]);
   	rast.rasterize(x[5], y[5], x[0], y[0]);

    std::vector<std::vector<double>> img;
    rast.render(img);

	writeTGA(img, "LL.tga");
	{
		std::string name("LLQ.ps");
		Tree2PS t2ps;
		t2ps.write(rast,name);
	}
}

TEST(RasterizeTest, Letter_L) {

	float x[6];
	float y[6];
	
	x[0] = y[0] = 0.;
	x[1] = 1.0; y[1] = 0.;
	x[2] = 1.0; y[2] = 0.35;
	x[3] = 0.45; y[3] = 0.25;
	x[4] = 0.25; y[4] = 1.0;
	x[5] = 0.0; y[5] = 1.0;

	int depth = 7;
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
	writeTGA(img, "letterL.tga");
	{
		std::string name("letterLQ.ps");
		Tree2PS t2ps;
		t2ps.write(rast,name);
	}
}

TEST(RasterizeTest, Triangle) {

	float x[3];
	float y[3];
	
	x[0] = 0.75; y[0] = 1.;
	x[1] = 0.25; y[1] = 1.;
	x[2] = 0.5; y[2] = 0.0;

	int depth = 5;
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
	{
		std::string name("triangleQ.ps");
		Tree2PS t2ps;
		t2ps.write(rast,name);
	}
}

TEST(RasterizeTest, TriangleRT) {

	float x[3];
	float y[3];
	
	x[0] = 0.0; y[0] = 0.;
	x[1] = 1.0; y[1] = 0.;
	x[2] = 0.0; y[2] = 1.0;

	int depth = 7;
	Rasterer2d rast(depth);
	for( int i = 0; i < 3; i++ ) {
   		rast.rasterize(x[i], y[i], x[(i+1)%3], y[(i+1)%3]);
	}
	std::vector<std::vector<double>> img;
	rast.render(img);
	rast.toPostscript(img, "triangleRT.ps");
	{
		std::string name("triangleRTQ.ps");
		Tree2PS t2ps;
		t2ps.write(rast,name);
	}
}

TEST(RasterizeTest, Shape) {

	float x[16];
	float y[16];

	x[0] = 0.3; y[0] = 0.0;	
	x[1] = 0.5; y[1] = 0.2;
	x[2] = 0.7; y[2] = 0.0;
	x[3] = 0.7; y[3] = 0.3;
	x[4] = 1.0; y[4] = 0.3;

	x[5] = 0.8; y[5] = 0.5;
	x[6] = 1.0; y[6] = 0.7;
	x[7] = 0.7; y[7] = 0.7;
	x[8] = 0.7; y[8] = 1.0;
	x[9] = 0.5; y[9] = 0.8;
	x[10] = 0.3; y[10] = 1.0;
	x[11] = 0.3; y[11] = 0.7;
	x[12] = 0.0; y[12] = 0.7;
	x[13] = 0.2; y[13] = 0.5;
	x[14] = 0.0; y[14] = 0.3;
	x[15] = 0.3; y[15] = 0.3;

	int depth = 6;
	Rasterer2d rast(depth);
	for( int i = 0; i < 16; i++ ) {
   		rast.rasterize(x[i], y[i], x[(i+1)%16], y[(i+1)%16]);
	}
	std::vector<std::vector<double>> img;
	rast.render(img);

	rast.toPostscript(img, "Shape.ps");
	{
		std::string name("ShapeQ.ps");
		Tree2PS t2ps;
		t2ps.write(rast,name);
	}
}

TEST(RasterizeTest, Star) {
    
    int n = 50;
    std::vector <double> x, y;

	float r1 = 0.47;
	float r2 = 0.1;
   
	makeStar( 0.5, 0.5, x, y, r1, r2, n );  
    int depth = 10;
    Rasterer2d rast(depth);
    for( int i = 0; i < n; i++ ) {
        rast.rasterize(x[i], y[i], x[(i+1)%n], y[(i+1)%n]);
    }
    std::vector<std::vector<double>> img;
    rast.render(img);
	writeTGA(img, "star.tga");
    {
        std::string name("starQ.ps");
        Tree2PS t2ps;
        t2ps.write(rast,name);
    }
}

TEST(RasterizeTest, Hexagon) {

	float x[6];
	float y[6];
	
	x[0] = 0.25; y[0] = 0.;
	x[1] = 0.75; y[1] = 0.;
	x[2] = 1.0; y[2] = 0.5;
	x[3] = 0.75; y[3] = 1.0;
	x[4] = 0.25; y[4] = 1.0;
	x[5] = 0.0; y[5] = 0.5;

	int depth = 8;
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
	writeTGA(img, "hexagon.tga");
	{
		std::string name("hexagonQ.ps");
		Tree2PS t2ps;
		t2ps.write(rast,name);
	}
}

TEST(RasterizeTest, LetterT) {

    float x[8], y[8];

	x[0] = 0.0; y[0] = 1.0;
	x[1] = 0.0; y[1] = 0.8;
	x[2] = 0.4; y[2] = 0.8;
	x[3] = 0.4; y[3] = 0.0;
	x[4] = 0.6; y[4] = 0.0;
	x[5] = 0.6; y[5] = 0.8;
	x[6] = 1.; y[6] = 0.8;
	x[7] = 1.0; y[7] = 1.0;

	int depth = 7;
    Rasterer2d rast(depth);
    for( int i = 0; i < 8; i++ )
        rast.rasterize(x[i],y[i],x[(i+1)%8], y[(i+1)%8]);

    std::vector<std::vector<double>> img;
    rast.render(img);

    rast.toPostscript(img, "letterT.ps");
    writeTGA(img, "letterT.tga");
    {
        std::string name("letterTQ.ps");
        Tree2PS t2ps;
        t2ps.write(rast,name);
	}
}

TEST(RasterizeTest, Shape1) {
    
    float x[7], y[7];
    
    x[0] = 0.; y[0] = 0.;
    x[1] = 1.0; y[1] = 0.0;
    x[2] = 0.5; y[2] = 0.5;
    
    x[3] = 0.5; y[3] = 0.75;
    x[4] = 0.25; y[4] = 0.75;
    x[5] = 0.25; y[5] = 1.0;
    x[6] = 0.0; y[6] = 1.0;
    
    int depth = 7;
    Rasterer2d rast(depth);
    for( int i = 0; i < 7; i++ )
        rast.rasterize(x[i],y[i],x[(i+1)%7], y[(i+1)%7]);
    
    std::vector<std::vector<double>> img;
    rast.render(img);
    writeTGA(img, "shape1.tga");

    {
        std::string name("shape1Q.ps");
        Tree2PS t2ps;
        t2ps.write(rast,name);
    }
}

TEST(RasterizeTest, Shape2) {

   float x[4], y[4];
   x[0] = 0.; y[0] = 0.;
   x[1] = 1.0; y[1] = 0.0;
   x[2] = 1.0; y[2] = 1.;
   x[3] = 0.0; y[3] = 1.;

    int depth = 12;
    Rasterer2d rast(depth);
    for( int i = 0; i < 4; i++ )
        rast.rasterize(x[i],y[i],x[(i+1)%4], y[(i+1)%4]);

   x[0] = 0.25; y[0] = 0.25;
   x[1] = 0.25; y[1] = 0.75;
   x[2] = 0.75; y[2] = 0.75;
   x[3] = 0.75; y[3] = 0.25;

   for( int i = 0; i < 4; i++ )
        rast.rasterize(x[i],y[i],x[(i+1)%4], y[(i+1)%4]);

    std::vector<std::vector<double>> img;
    rast.render(img);
    writeTGA(img, "shape2.tga");
    {
        std::string name("shape2Q.ps");
        Tree2PS t2ps;
        t2ps.write(rast,name);
    }
}

TEST(RasterizeTest, Shape3) {

   float x[4], y[4];
   x[0] = 0.; y[0] = 0.;
   x[1] = 1.0; y[1] = 0.0;
   x[2] = 1.0; y[2] = 1.;
   x[3] = 0.0; y[3] = 1.;

    int depth = 12;
    Rasterer2d rast(depth);
    for( int i = 0; i < 4; i++ )
        rast.rasterize(x[i],y[i],x[(i+1)%4], y[(i+1)%4]);

	// inner
	{
		std::vector<double> px, py;
		double r1 = 0.2;
		double r2 = 0.05;	
		int n = 50;
		makeStar( 0.25, 0.25,px, py, r1, r2, 50);
		for(int i = 0; i < n; i++ ) 
        	rast.rasterize( px[(i+1)%n], py[(i+1)%n],px[i], py[i]);
	}

	// inner
	{
		std::vector<double> px, py;
		double r1 = 0.22;
		double r2 = 0.15;	
		int n = 30;
		makeStar( 0.75, 0.75,px, py, r1, r2, n);
		for(int i = 0; i < n; i++ ) 
        	rast.rasterize( px[(i+1)%n], py[(i+1)%n],px[i], py[i]);
	}
	{
		int n = 3;
   		x[0] = 0.6; y[0] = 0.1;
   		x[1] = 0.75; y[1] = 0.4;
   		x[2] = 0.9; y[2] = 0.1;
		for(int i = 0; i < n; i++ ) 
        	rast.rasterize( x[i], y[i], x[(i+1)%n], y[(i+1)%n] );

	}
	{
		float xx[] = {
			0.234578,
			0.242965, 0.246212, 0.25, 0.248377, 0.246212, 0.242154, 0.234578, 0.22592, 0.220238, 0.214015, 0.209145, 0.208063,
			0.208063, 0.209957, 0.217262, 0.224567, 0.22592, 0.227814, 0.231331, 0.233225, 0.237284, 0.241342, 0.239448, 0.231061,
			0.223755, 0.215909, 0.209686, 0.195346, 0.186147, 0.182359, 0.176677, 0.170455, 0.166937, 0.161797, 0.158279,
			0.155574, 0.155032, 0.158279, 0.16369, 0.162608, 0.148268, 0.146104, 0.146104, 0.145833, 0.147727, 0.148268,
			0.150703, 0.152868, 0.148539, 0.135823, 0.12013, 0.104167, 0.103084, 0.0968615, 0.090368, 0.0873918, 0.0838745,
			0.08171, 0.0844156, 0.0890152, 0.094697, 0.0979437, 0.092803, 0.084145, 0.0776515, 0.0733225, 0.0676407, 0.0676407,
			0.0716991, 0.0727814, 0.0662879, 0.0579004, 0.0508658, 0.0446429, 0.0359848, 0.0286797, 0.0286797, 0.0308442, 0.0316558,
			0.0251623, 0.0205628, 0.0205628, 0.0205628, 0.030303, 0.0373377, 0.0403139, 0.0462662, 0.0573593, 0.0611472, 0.0641234,
			0.0587121, 0.0443723, 0.021645, 0.00568182, 0.000541126, 0.0105519, 0.0221861, 0.0211039, 0.0183983, 0.0162338, 0.0151515,
			0.0116342, 0.00730519, 0.00811688, 0.0183983, 0.0419372, 0.0714286, 0.0873918, 0.102543, 0.117695, 0.143669, 0.166937,
			0.17289, 0.175595, 0.179383, 0.181818, 0.185336, 0.189664, 0.196158, 0.215097, 0.226732, 0.228084, 0.228626, 0.231061,
			0.233225, 0.234578, 0.234578};


		float yy[] = { 0.566306,
			0.576621, 0.591847, 0.607564, 0.613458, 0.617878, 0.62279, 0.619352, 0.611493, 0.612475, 0.619843, 0.623281,
			0.630157, 0.648821, 0.666503, 0.694499, 0.716601, 0.722004, 0.729371, 0.733792, 0.73723, 0.73723, 0.742633,
			0.749509, 0.749509, 0.746071, 0.730354, 0.715619, 0.685658, 0.665521, 0.655697, 0.658153, 0.662574, 0.674361,
			0.691552, 0.715128, 0.721513, 0.732809, 0.734774, 0.738703, 0.749018, 0.749018, 0.740668, 
			0.727898, 0.722004, 0.700884, 0.69057, 0.668959, 0.654715, 0.652259, 0.653242, 0.646365,
			0.634578, 0.639489, 0.660609, 0.679764, 0.688605, 0.697937, 0.709725, 0.731336, 0.737721,
			0.736248, 0.742141, 0.749018, 0.748035, 0.739195, 0.721513, 0.701866, 0.695481, 0.68222, 0.657662,
			0.664047, 0.67387, 0.680747, 0.686641, 0.702849, 0.728389, 0.738703, 0.744106, 0.749018, 0.748035,
			0.742141, 0.734774, 0.727898, 0.687623, 0.666503, 0.663556, 0.65668, 0.627701,
			0.602652, 0.591356, 0.59332, 0.599705, 0.592829, 0.564833, 0.527505, 0.5, 0.513752,
			0.528978, 0.534872, 0.529961, 0.518664, 0.515226, 0.527014, 0.546169, 0.570236, 0.582024,
			0.562868, 0.557466, 0.556974, 0.564342, 0.56778, 0.573183, 0.571218, 0.566798, 0.566306, 0.564342,
			0.564833, 0.568271, 0.568762, 0.561395, 0.55943, 0.557957, 0.559921, 0.55501, 0.555992, 0.566306, 0.566306};

		int n = 128;
		for( int i = 0; i < 128; i++ ) {
        		rast.rasterize( xx[(i+1)%n], yy[(i+1)%n],xx[i], yy[i]);
		}
	

	}

    std::vector<std::vector<double>> img;
    rast.render(img);
    writeTGA(img, "shape3.tga");
    {
        std::string name("shape3Q.ps");
        Tree2PS t2ps;
        t2ps.write(rast,name);
	}
}

int main(int argc, char **argv ) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
