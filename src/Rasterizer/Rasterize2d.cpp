#include "Rasterize2d.h"
#include "Rasterer2d.h"

#include "QuadTree.h"
#include "vtkPolyData.h"
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkPolygon.h"
#include "vtkHyperOctree.h"
#include "vtkHyperOctreeCursor.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <cstddef>

vtkStandardNewMacro(Rasterize2d);

using namespace std;


///////
Rasterize2d::Rasterize2d()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
}

Rasterize2d::~Rasterize2d()
{
}

void
Rasterize2d::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

int
Rasterize2d::ProcessRequest(vtkInformation* request,
			    vtkInformationVector** inputVector,
			    vtkInformationVector* outputVector)
{

	cout << "Process Request\n";

	// generate data.
	if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
	{
		return this->RequestData(request, inputVector, outputVector);
	}
	return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

vtkHyperOctree*
Rasterize2d::GetOutput()
{
	return this->GetOutput(0);
}

vtkHyperOctree* 
Rasterize2d::GetOutput(int port)
{
	return vtkHyperOctree::SafeDownCast(this->GetOutputDataObject(port));
}

int
Rasterize2d::RequestData( vtkInformation*,
			 vtkInformationVector** inputVector,
			 vtkInformationVector* outputVector )
{
	cout << "request data...\n";
#if 0
	float cvalues[3];	
	float v[4];

	v[0] = 0.4;
	v[1] = 0.4;
	
	v[2] = 0.8;
	v[3] = 1.0;
	cvalues[0] = 0.0; cvalues[1] = 0.0; cvalues[2] = 0.0;
	calc_coeffs_line(cvalues, v, 0, 0);
	cout << "Quad step1: " << cvalues[0] << "," << cvalues[1] << "," << cvalues[2] << "\n";	

	v[0] = 0.8;
	v[1] = 0.0;
	
	v[2] = 1.0;
	v[3] = 0.3;
	float vv[3];
	vv[0] = 0.; vv[1] = 0.; vv[2] = 0.;
	calc_coeffs_line(cvalues, v, 0, 1);
	calc_coeffs_line(vv, v, 0, 1);
	cout << "Quad step2: " << vv[0] << "," << vv[1] << "," << vv[2] << "\n";	

	v[0] = 0.0;
	v[1] = 0.3;
	
	v[2] = 0.2;
	v[3] = 0.6;
	calc_coeffs_line(cvalues, v, 1, 1);
	vv[0] = 0.; vv[1] = 0.; vv[2] = 0.;
	calc_coeffs_line(vv, v, 1, 1);
	cout << "Quad step3: " << vv[0] << "," << vv[1] << "," << vv[2] << "\n";	
	cout << "Quad 0+3: " << cvalues[0] << "," << cvalues[1] << "," << cvalues[2] << "\n";	
	cout <<"+++++++++++++++++\n";
#endif
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// code is the quadrant containing the point
	// 00 => x < 0.5, y < 0.5
	// 01 => x >= 0.5, y < 0.5
	// 10 => x < 0.5, y >= 0.5
	// 11 => x >= 0.5, y >= 0.5

	// insert edge into quadtree.

	// extract n-gons to rasterize
	vtkPolyData *input = vtkPolyData::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()));
	if( input  ) {

		Rasterer2d rastOp( Depth );	

		input->GetPolys()->InitTraversal();
  		vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
#if 0
		{
			Rasterer2d rastOp1( 1 );	
			double p1[2], p2[2];
			p1[0] = 0.2;
			p1[1] = 0.2;
			p2[0] = 0.4;
			p2[1] = 0.4;

			rastOp1.rasterize(p1[0], p1[1], p2[0], p2[1]);
			rastOp1.debug();
		}

		{

			Rasterer2d rastOp1( 1 );	
			double p1[2], p2[2];
			p1[0] = 0.7;
			p1[1] = 0.2;
			p2[0] = 0.9;
			p2[1] = 0.4;

			rastOp1.rasterize(p1[0], p1[1], p2[0], p2[1]);
			rastOp1.debug();
		}

		{

			Rasterer2d rastOp1( 1 );	
			double p1[2], p2[2];
			p1[0] = 0.2;
			p1[1] = 0.7;
			p2[0] = 0.4;
			p2[1] = 0.9;

			rastOp1.rasterize(p1[0], p1[1], p2[0], p2[1]);
			rastOp1.debug();
		}

		{

			Rasterer2d rastOp1( 1 );	
			double p1[2], p2[2];
			p1[0] = 0.7;
			p1[1] = 0.7;
			p2[0] = 0.9;
			p2[1] = 0.9;

			rastOp1.rasterize(p1[0], p1[1], p2[0], p2[1]);
			rastOp1.debug();
		}

		{
			Rasterer2d rastOp1( 1 );	
			double p1[2], p2[2];
			p1[0] = 0.2;
			p1[1] = 0.2;
			p2[0] = 0.8;
			p2[1] = 0.8;

			rastOp1.rasterize(p1[0], p1[1], p2[0], p2[1]);
			rastOp1.debug();
		}
#endif
		{
			Rasterer2d rastOp1( 1 );	
			double p1[2], p2[2];
			p1[0] = 0.2;
			p1[1] = 0.2;
			p2[0] = 0.6;
			p2[1] = 0.8;
			rastOp1.rasterize(p1[0], p1[1], p2[0], p2[1]);
			rastOp1.debug();
		}
#if 0
		// traverse polygon
		while(input->GetPolys()->GetNextCell(idList)) {

			// Haar each n-gon edge.
			vtkIdType nids = idList->GetNumberOfIds();
			for( vtkIdType vid = 0; vid < nids; vid++ )
			{
				vtkIdType id1 = idList->GetId(vid);
				double p1[3];
				input->GetPoint(id1, p1);
				double p2[3];
				vtkIdType id2 = (id1+1)%nids;
				input->GetPoint(id2, p2);
				rastOp.rasterize(p1[0], p1[1], p2[0], p2[1]);
			}	
			cout << "Haar'ed the n-gon edges:: " << idList->GetNumberOfIds() << "\n"; 	
		}
		rastOp.debug();
#endif
	}
	return 1;
}

int
Rasterize2d::FillOutputPortInformation(
	int vtkNotUsed(port),
	vtkInformation* info )
{
	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkHyperOctree");
	return 1;
}

int
Rasterize2d::FillInputPortInformation(
	int vtkNotUsed(port),
	vtkInformation* info )
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}

void
Rasterize2d::SetInputData(vtkDataObject* input)
{
	this->SetInputData(0, input);
}

void
Rasterize2d::SetInputData(int index, vtkDataObject* input)
{
	this->SetInputDataInternal(index, input);
}

