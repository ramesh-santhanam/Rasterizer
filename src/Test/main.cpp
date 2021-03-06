#include <iostream>
#include <Rasterize2d.h>
#include <Rasterer2d.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>


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
                {
                        calc_coeffs00(coeffs,v);
                }
                else
                {
                        calc_coeffs01(coeffs,v);
                }
        }
        else
        {
                if ( j == 0 )
                {
                        calc_coeffs10(coeffs,v);
                }
                else
                {
                        calc_coeffs11(coeffs,v);
                }
        }
}

bool 
equivalent(double x, double y, double tol = 1.0e-05 )
{
	return fabs(x-y) < tol;
}

void 
runTest1()
{
	float cval[3];
	double px, py;
    double qx, qy;

    // Quad 0
    float v[4];
    {
		cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
        v[0] = 0.4; v[1] = 0.3;
        v[2] = 0.2; v[3] = 0.4 ;
        calc_coeffs_line(cval, v, 0, 0 );

		px = 0.2; py = 0.15;
        qx = 0.1; qy = 0.2;
        Rasterer2d rast(1);
        rast.rasterize(px, py, qx, qy);
        assert(rast.m_root);

		assert(equivalent(rast.m_root->m_haar[0], cval[0] ));
		assert(equivalent(rast.m_root->m_haar[1], cval[1] ));
		assert(equivalent(rast.m_root->m_haar[2], cval[2] ));
#if 0
		cout << "Haar0 (T1:Q0):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
		cout << "Haar1 (T1:Q0):" << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
		cout << "Haar2 (T1:Q0):" << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
   }

	// Quad 1
	{
		cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
        v[0] = 0.4; v[1] = 0.3;
        v[2] = 0.2; v[3] = 0.4 ;
        calc_coeffs_line(cval, v, 1, 0 );

		px = 0.7; py = 0.15;
        qx = 0.6; qy = 0.2;
        Rasterer2d rast(1);
        rast.rasterize(px, py, qx, qy);
        assert(rast.m_root);

		assert(equivalent(rast.m_root->m_haar[0], cval[0] ));
		assert(equivalent(rast.m_root->m_haar[1], cval[1] ));
		assert(equivalent(rast.m_root->m_haar[2], cval[2] ));
#if 0
		cout << "Haar0 (T1:Q1):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
		cout << "Haar1 (T1:Q1):" << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
		cout << "Haar2 (T1:Q1):" << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif

	}

	// Quad 2
	{
		cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
        v[0] = 0.4; v[1] = 0.3;
        v[2] = 0.2; v[3] = 0.4 ;
        calc_coeffs_line(cval, v, 0, 1 );

		px = 0.2; py = 0.65;
        qx = 0.1; qy = 0.7;
        Rasterer2d rast(1);
        rast.rasterize(px, py, qx, qy);
        assert(rast.m_root);

		assert(equivalent(rast.m_root->m_haar[0], cval[0] ));
		assert(equivalent(rast.m_root->m_haar[1], cval[1] ));
		assert(equivalent(rast.m_root->m_haar[2], cval[2] ));
#if 0
		cout << "Haar0 (T1:Q2) : " << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
		cout << "Haar1 (T1:Q2) : " << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
		cout << "Haar2 (T1:Q2) : " << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
	}

	// Quad 3
	{
		cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
        v[0] = 0.4; v[1] = 0.3;
        v[2] = 0.2; v[3] = 0.4 ;
        calc_coeffs_line(cval, v, 1, 1 );

		px = 0.7; py = 0.65;
        qx = 0.6; qy = 0.7;
        Rasterer2d rast(1);
        rast.rasterize(px, py, qx, qy);
        assert(rast.m_root);

		assert(equivalent(rast.m_root->m_haar[0], cval[0] ));
		assert(equivalent(rast.m_root->m_haar[1], cval[1] ));
		assert(equivalent(rast.m_root->m_haar[2], cval[2] ));
#if 0
		cout << "Haar0 (T1:Q3) : " << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
		cout << "Haar1 (T1:Q3) : " << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
		cout << "Haar2 (T1:Q3) : " << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
	}
	
}

void
runTest2() 
{
	float cval[3];
	double px, py;
    double qx, qy;

    // Quad 0
    float v[4];
    {
		// Q0 - Q2
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
#if 0
		cout << "Haar0 (T2:Q0-Q2):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
		cout << "Haar1 (T2:Q0-Q2): " << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
		cout << "Haar2 (T2:Q0-Q2): " << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
		assert(equivalent(rast.m_root->m_haar[0], cval[0] ));
		assert(equivalent(rast.m_root->m_haar[1], cval[1] ));
		assert(equivalent(rast.m_root->m_haar[2], cval[2] ));
	}

    {
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
#if 0
		cout << "Haar0 (T2:Q1-Q3):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
		cout << "Haar1 (T2:Q1-Q3): " << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
		cout << "Haar2 (T2:Q1-Q3): " << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
		assert(equivalent(rast.m_root->m_haar[0], cval[0] ));
		assert(equivalent(rast.m_root->m_haar[1], cval[1] ));
		assert(equivalent(rast.m_root->m_haar[2], cval[2] ));
	}

    {
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
        assert(rast.m_root);

		cout << "Haar0 (T2:Q0-Q1):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
		cout << "Haar1 (T2:Q0-Q1): " << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
		cout << "Haar2 (T2:Q0-Q1): " << rast.m_root->m_haar[2] << "," << cval[2] << "\n";

		assert(equivalent(rast.m_root->m_haar[0], cval[0] ));
		assert(equivalent(rast.m_root->m_haar[1], cval[1] ));
		assert(equivalent(rast.m_root->m_haar[2], cval[2] ));
	}

    {
		// Q2 - Q3
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
#if 0
		cout << "Haar0 (T2:Q2-Q3):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
		cout << "Haar1 (T2:Q2-Q3): " << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
		cout << "Haar2 (T2:Q2-Q3): " << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
		assert(equivalent(rast.m_root->m_haar[0], cval[0] ));
		assert(equivalent(rast.m_root->m_haar[1], cval[1] ));
		assert(equivalent(rast.m_root->m_haar[2], cval[2] ));
	}
}

void
runTest3()
{
	float cval[3];
	double px, py;
    double qx, qy;
    float v[4];
	{
		cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
        v[0] = 0.0; v[1] = 0.5;
        v[2] = 0.5; v[3] = 1.0;
        calc_coeffs_line(cval, v, 0, 0 );

        v[0] = 0.5; v[1] = 0.0;
        v[2] = 1.0; v[3] = 0.5 ;
        calc_coeffs_line(cval, v, 0, 1 );
		
        v[0] = 0.0; v[1] = 0.5;
        v[2] = 0.3; v[3] = 0.8 ;
        calc_coeffs_line(cval, v, 1, 1 );

		px = 0.0; py = 0.25;
        qx = 0.65; qy = 0.9;
        Rasterer2d rast(1);
        rast.rasterize(px, py, qx, qy);
        assert(rast.m_root);
#if 0
		cout << "Haar0 (T3:Q0-Q2-Q3):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
		cout << "Haar1 (T3:Q0-Q2-Q3): " << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
		cout << "Haar2 (T3:Q0-Q2-Q3): " << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif
		assert(equivalent(rast.m_root->m_haar[0], cval[0] ));
		assert(equivalent(rast.m_root->m_haar[1], cval[1] ));
		assert(equivalent(rast.m_root->m_haar[2], cval[2] ));
	}

	{
		cval[0] = 0.; cval[1] = 0.; cval[2] = 0.;
        v[0] = 0.50; v[1] = 0.50;
        v[2] = 1.0; v[3] = 1.0;
        calc_coeffs_line(cval, v, 0, 0 );

        v[0] = 0.0; v[1] = 0.0;
        v[2] = 0.5; v[3] = 0.5 ;
        calc_coeffs_line(cval, v, 1, 1 );
		
		px = 0.25; py = 0.25;
        qx = 0.75; qy = 0.75;
        Rasterer2d rast(1);
        rast.rasterize(px, py, qx, qy);
        assert(rast.m_root);
#if 0
		cout << "Haar0 (T3:Q0-Q3):" << rast.m_root->m_haar[0] << "," << cval[0] << "\n";
		cout << "Haar1 (T3:Q0-Q3): " << rast.m_root->m_haar[1] << "," << cval[1] << "\n";
		cout << "Haar2 (T3:Q0-Q3): " << rast.m_root->m_haar[2] << "," << cval[2] << "\n";
#endif

		assert(equivalent(rast.m_root->m_haar[0], cval[0] ));
		assert(equivalent(rast.m_root->m_haar[1], cval[1] ));
		assert(equivalent(rast.m_root->m_haar[2], cval[2] ));
	}
}

int
main()
{
	cout << "Test Rasterizer....\n";

	// lines confined to a single quadrant
	runTest1();

	// lines confined to two quadrant
	runTest2();

	// lines confined to three quadrant
	runTest3();

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->InsertNextPoint(0.0, 0.0, 0.0);
	points->InsertNextPoint(1.0, 0.0, 0.0);
	points->InsertNextPoint(1.0, 1.0, 0.0);
	points->InsertNextPoint(0.0, 1.0, 0.0);
	points->InsertNextPoint(0.5, 0.75, 0.0);
	points->InsertNextPoint(0.25, 0.25, 0.0);
	points->InsertNextPoint(0.25, 0.75, 0.0);

	// outer polygon
	vtkSmartPointer<vtkPolygon> polyO = vtkSmartPointer<vtkPolygon>::New();
	polyO->GetPointIds()->SetNumberOfIds(4);
	polyO->GetPointIds()->SetId(0,0);
	polyO->GetPointIds()->SetId(1,1);
	polyO->GetPointIds()->SetId(2,2);
	polyO->GetPointIds()->SetId(3,3);
	vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();
	polygons->InsertNextCell(polyO);

	// inner polygon
#if defined INNER_POLYGON
	vtkSmartPointer<vtkPolygon> polyI = vtkSmartPointer<vtkPolygon>::New();
	polyI->GetPointIds()->SetNumberOfIds(3);
	polyI->GetPointIds()->SetId(0,4);
	polyI->GetPointIds()->SetId(1,5);
	polyI->GetPointIds()->SetId(2,6);
	polygons->InsertNextCell(polyI);
#endif
	vtkSmartPointer<vtkPolyData> pdata = vtkSmartPointer<vtkPolyData>::New();
	pdata->SetPoints(points);
	pdata->SetPolys(polygons);

	vtkSmartPointer<Rasterize2d> raster = vtkSmartPointer<Rasterize2d>::New();
	raster->SetInputData(pdata);
	raster->SetDepth(1);
	raster->Update();
	return 0;
}
