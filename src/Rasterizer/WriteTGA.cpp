#include "WriteTGA.h"

#include <cstddef>
#include <cassert>
#include <iostream>
#include <cmath>

#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctype.h>

using namespace std;

bool
writeTGA(const std::vector< std::vector<double>>& img, const std::string& name )
{
	assert(img.size());

	if(img.size() == 0 )
		return false;

	std::ofstream file;
	file.open(name.c_str(),ios::binary);
	assert(file.is_open());

	if( ! file.is_open() )
		return false;

	short height = img.size();
	short width = img[0].size();

	// write header.
	{
		file.put(0);
		file.put(0);
		file.put(3); // monochrome
		file.put(0); file.put(0);
		file.put(0); file.put(0);
		file.put(0);
		file.put(0); file.put(0); // X - origin
		file.put(0); file.put(0); // Y - origin
		file.put(width&0x00FF);
		file.put((width&0xFF00)/256);
		file.put(height&0x00FF);
		file.put((height&0xFF00)/256);
		file.put(8);	
		file.put(0);	
	}

	auto clamp = [] (float v) -> float {
    	if( v < 0.0 )
            return 0.0;
        if( v > 1.0 )
            return 1.0;
        return v;
    };

	for( int r = 0; r < img.size(); r++ ) {
		assert(img[r].size() == width);
		for( int j = 0; j < img[r].size(); j++ ) {
            float val = img[r][j];
            val = clamp(val);
            char  gv = (int)((1 - val ) * 255.0);
			file << hex << gv;
		}
	}
	file.close();	
	return true;
}

