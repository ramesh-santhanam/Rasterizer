# Rasterizer
Rasterize 2d polygon + Meshes with antialiasing using non-standard Haar Wavelets
# Paper
The implementation is based on the 2011 Eurographics paper "Wavelet Rasterization"  authored by Manson & Schaefer 
## Build
git clone  -r https://github.com/ramesh-santhanam/Rasterizer.git
mkdir build  
cd build  
cmake -DVTK_Group_StandAlone:Bool=OFF -DVTK_Group_Rendering:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DVTK_WITHOUT_QT:BOOL=ON -DModule_vtkCommonDataModel:BOOL=ON -DModule_vtkFiltersGeometry:BOOL=ON  ..
make 
