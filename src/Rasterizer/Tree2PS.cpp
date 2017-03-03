#include "Tree2PS.h"
#include "QuadTree.h"
#include "Rasterer2d.h"

#include <cstddef>
#include <cassert>
#include <iostream>
#include <cmath>

#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctype.h>

#include <set>
#include <map>
#include <vector>
#include <string>

using namespace std;

typedef std::set<double> DblSet;
typedef std::map<double, DblSet> MultiLine;

class ExportPS {

	public:
	ExportPS(std::string& name);
	~ExportPS();

	public:
	bool write(const MultiLine& lines);

	public:
	std::ofstream m_file;
	int m_w;
	int m_h;

	double m_ox;
	double m_oy;
	double m_sx;
	double m_sy;
};

ExportPS::ExportPS(std::string& name)
{
	m_file.open(name.c_str());
	assert(m_file.is_open());

	m_w = 595;
	m_h = 842;

	m_ox = 0.25*m_w;	
	m_oy = 0.25*m_h;	

	m_sx = 0.5*m_w;
	m_sy = 0.5*m_h;	
	if( m_file.is_open() ) {

		m_file << "%!PS-Adobe " << '\n'
        	<< "%%Title: Rasterizer2D" << std::endl;

		m_file << "/m {moveto} bind def"      << std::endl
         << "/l {lineto} bind def"      << std::endl
         << "/s {setrgbcolor} bind def" << std::endl
         << "/sg {setgray} bind def"    << std::endl
         << "/slw {setlinewidth} bind def" << std::endl
         << "/sls {setlinestyle} bind def" << std::endl
         << "/gs {gsave} bind def"      << std::endl
         << "/gr {grestore} bind def"   << std::endl
         << "/pf {lineto closepath gsave 0.9 setgray fill grestore 1.1 setlinewidth 0.25 setgray stroke} bind def" << std::endl
		 << "/lx {lineto closepath stroke} bind def" << std::endl
         << "/lf {lineto closepath fill} bind def"   << std::endl
         << "/ci2 { 2 0 360 arc closepath stroke } bind def" << std::endl
         << "/cif2 { 2 0 360 arc closepath fill } bind def" << std::endl;

		m_file << "%%EndProlog" << '\n';


		// draw box.
		double ox = 0.25*m_w;
		double oy = 0.25*m_h;
		double sx = 0.5 * m_w;
		double sy = 0.5 * m_h;
		m_file << ox-5 << " " << oy-5;
    	m_file << " m " ;
    	m_file << ox+sx+5 << " " << oy-5;
    	m_file << " l " ;
    	m_file << ox+sx+5 << " " << oy+sy+5;
    	m_file << " l " ;
    	m_file << ox-5 << " " << oy+sy+5;
    	m_file << " lx " << '\n';
	}
}

ExportPS::~ExportPS()
{
	m_file <<  std::endl << "showpage" << std::endl;
	m_file.close();	
}
bool
ExportPS::write(const MultiLine& lines)
{
	assert(m_file.is_open());
	if( ! m_file.is_open() )
		return false;


	auto xform = [&] (double& x, double &y) {
		x = m_ox + x * m_sx;
		y = m_oy + y * m_sy;
	};
	auto writeLine = [&] (double x1, double y1, double x2, double y2 ) {
		xform(x1, y1);
		m_file << x1 << " " << y1 << " " << " m " <<  " " ;
		xform(x2, y2);
		m_file << x2 << " " << y2 << " " << " l " ;
		m_file << " stroke " << std::endl;
	};
	// write lines. - x = const.
	m_file << "%const grid lines" << std::endl;
	m_file << "0.45" << " slw " << std::endl;

	for( auto s = lines.begin(); s != lines.end(); s++ ) {
		double v = s->first;
		const DblSet& ords = s->second;

		double x = v;
		double y = *ords.begin(); 
		for( int i = 1; i < ords.size(); i++ ) {
			double yn = *std::next(ords.begin(), i);
			// write line.
			writeLine(x, y, x, yn);	
			y = yn;	
		}
	} 

	return true;
}


///////////////

Tree2PS::Tree2PS()
{
}

Tree2PS::~Tree2PS()
{
}

void
Tree2PS::write( Rasterer2d& rast, std::string& name)
{
	MultiLine lineX; // x - ordinate is key, lines with slope infinity.
	MultiLine lineY; // y -ordinate is key, lines with slope zero.

	auto quadSize = [] (int level ) -> double {
		double size = 1.;
		for( int i = 0; i < level; i++ )
			size = size/2.;
		return size;	
	};

	std::vector<int> qcodes;
	auto quadPosition = [&] (QuadNode* n, double& x, double& y ) {
		x = 0.;
		y = 0.;
		if( n->getParent() == nullptr )
			return;		

		qcodes.resize(n->m_level+1,-1);
		QuadNode* qc = n;
		do {
			int index = -1;
			QuadNode* qp = qc->getParent();
			if( qp == nullptr )
				index = 0;
			else {
				for( int i = 0; i < 4; i++ ) {
					if( qp->getSubQuad(i) == qc ) {
						index = i;
						break;
					}
				}
			}
			qcodes[qc->m_level] = index;
			qc = qp;
		} while( qc != nullptr );

		// transform (0,0) to right location using codes
		// scaling by 0.5 each time along the way.
		double s = 1.0;
	
		// skip the first as quad would be outer one.	
		for( int i = 1; i < qcodes.size(); i++ ) {

			assert(qcodes[i] != -1);
			double sx = (qcodes[i] & 0x01) ? 0.5*s : 0;
			x = x + sx;			
			double sy = (qcodes[i] & 0x02) ? 0.5*s : 0;
			y = y + sy;
			s *= 0.5;
		}
		//cout << "came out::" << qcodes[n->m_level] << ":size:" << s << "<" << x << "," << y << ">" << "\n";
	};

	auto nodeFunc = [&] (QuadNode* n)->void {
		assert(n);
		if( ! n )
			return;
		double  s = quadSize(n->m_level);		
		double x, y;
		quadPosition(n, x, y);
		
		// lineX
		lineX[x].insert(y);
		lineX[x].insert(y+s);
		// lineY
		lineY[y].insert(x);
		lineY[y].insert(x+s);
					
	};

	visitTree(rast.m_root, nodeFunc );

	// write PS.
	{
		ExportPS wtoPS(name);
		wtoPS.write(lineX);
		//wtoPS.write(lineY);
	}
	cout << " X lines: " << lineX.size();
	cout << " Y lines: " << lineY.size();
}

