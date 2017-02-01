#ifndef Rasterize2d_h
#define Rasterize2d_h

#include "vtkAlgorithm.h"
#include "QuadTree.h"

class vtkDataObject;
class vtkPolyData;
class vtkPolygon;
class vtkHyperOctree;

class Rasterize2d : public vtkAlgorithm {

	public:
		static Rasterize2d *New();
		vtkTypeMacro(Rasterize2d,vtkAlgorithm);
		void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

		vtkSetMacro( Depth, int );
		vtkGetMacro( Depth, int );

		vtkHyperOctree* GetOutput();
		vtkHyperOctree* GetOutput(int);

		virtual int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) VTK_OVERRIDE;

		// Does not handle setting up a connection
		void SetInputData(vtkDataObject *);
		void SetInputData(int, vtkDataObject *);

	protected:
		Rasterize2d();
		~Rasterize2d();

		int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);	
		int FillOutputPortInformation(int port, vtkInformation* info) VTK_OVERRIDE;
		int FillInputPortInformation(int port, vtkInformation* info) VTK_OVERRIDE;

		int Depth;

	private:
		void operator=(const Rasterize2d&) VTK_DELETE_FUNCTION;	
		Rasterize2d(const Rasterize2d&) VTK_DELETE_FUNCTION;	
};

#endif


