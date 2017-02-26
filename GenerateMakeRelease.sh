#!bin/sh
basepath="$(cd "$(dirname "$0")" && pwd)"
echo $basepath
builddir="$basepath/build"
echo $builddir
cmake -E make_directory $builddir
cd $builddir
cmake -DVTK_Group_StandAlone:Bool=OFF -DVTK_Group_Rendering:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DVTK_WITHOUT_QT:BOOL=ON -DModule_vtkCommonDataModel:BOOL=ON -DModule_vtkFiltersGeometry:BOOL=ON  .. 
cd $builddir
make 
