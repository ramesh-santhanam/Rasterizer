# Rasterizer
Rasterize 2d polygon + Meshes with antialiasing using non-standard Haar Wavelets
# Paper
The implementation is based on the 2011 Eurographics paper "Wavelet Rasterization"  authored by Manson & Schaefer 
# Mathematica
The analytical derivation for the computation of the Haar wavelet coefficient is done using Mathematica. The polygon contours can be defined using linear edges, quadratic Bezier or cubic Bezier. Also enclosed is the Haar coefficients for 3d trianglular meshes. The mathemtica files (.nb) are stored in a Mathematica folder. The key idea is to express the constant function of the enclosed area (volume) as the sum of Haar basis function. The coeffcients of the basis function is ascertained using the divergence theorem by equating the computation of the surface integral by its line integral (2d polygons), the volume integral by the surface integral (3d meshes).In the 2D-contour cases, horizontal and vertical lines will have their respective wavelet (10, 01, 11) coefficients to be 0.0
## Build
git clone  --recursive https://github.com/ramesh-santhanam/Rasterizer.git  
mkdir build  
cd build  
cmake -DVTK_Group_StandAlone:Bool=OFF -DVTK_Group_Rendering:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DVTK_WITHOUT_QT:BOOL=ON -DModule_vtkCommonDataModel:BOOL=ON -DModule_vtkFiltersGeometry:BOOL=ON  ..  
make 
