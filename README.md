# Rasterizer
Rasterize 2d polygon + Meshes with antialiasing using non-standard Haar Wavelets
# Paper
The implementation is based on the 2011 Eurographics paper "Wavelet Rasterization"  authored by Manson & Schaefer 
# Mathematica
The analytical derivation for the computation of the Haar wavelet coefficient is done using Mathematica. The polygon contours can be defined using linear edges, quadratic Bezier or cubic Bezier. To come is the derivation of the 3DHaar coefficient for triangle meshes in 3D. The mathemtica files (.nb) are stored in a Mathematica folder.
## Build
git clone  --recursive https://github.com/ramesh-santhanam/Rasterizer.git  
mkdir build  
cd build  
cmake -DVTK_Group_StandAlone:Bool=OFF -DVTK_Group_Rendering:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DVTK_WITHOUT_QT:BOOL=ON -DModule_vtkCommonDataModel:BOOL=ON -DModule_vtkFiltersGeometry:BOOL=ON  ..  
make 
