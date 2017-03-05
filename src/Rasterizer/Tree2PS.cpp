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
typedef std::vector<double> DblValues;
typedef std::map<double, DblValues> MultiLine;

class ExportPS {

	public:
	ExportPS(std::string& name);
	~ExportPS();

	public:
	bool write(const MultiLine& lines, int);

	public:
	std::ofstream m_file;
	double m_w;
	double m_h;

	double m_ox;
	double m_oy;
	double m_hx;
	double m_hy;
};

ExportPS::ExportPS(std::string& name)
{
	m_file.open(name.c_str());
	assert(m_file.is_open());

	double pw = 8.5;
	double ph = 11.0;
	
	m_w = 8; // width in inches A4 size
	m_h = 8; // ht in inches.

	m_ox = 0.5 * pw;
	m_oy = 0.5 * ph; 

	m_hx = 0.5*m_w;
	m_hy = 0.5*m_h;	

	if( m_file.is_open() ) {

		m_file << "%!PS-Adobe " << '\n'
        	<< "%%Title: Rasterizer2D" << std::endl;

		m_file << "/m {moveto} bind def"      << std::endl
         << "/l {lineto} bind def"      << std::endl
         << "/sg {setgray} bind def"    << std::endl
         << "/slw {setlinewidth} bind def" << std::endl
         << "/gs {gsave} bind def"      << std::endl
         << "/gr {grestore} bind def"   << std::endl;
		m_file << "%%EndProlog" << '\n';


		m_file << " 72 72 scale " << endl;
		m_file << m_ox << " " << m_oy << " translate" << endl;
 		m_file << "0.5 72 div " << "slw" << endl;	
	}
}

ExportPS::~ExportPS()
{
	m_file <<  std::endl << "showpage" << std::endl;
	m_file.close();	
}
bool
ExportPS::write(const MultiLine& lines, int kind)
{
	assert(m_file.is_open());
	if( ! m_file.is_open() )
		return false;


	auto xform = [&] (double& x, double &y) {
		x = -m_hx + x * m_w;
		y = -m_hy + y * m_h;
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
	if( kind == 0 ) {
		m_file << "% X = const grid lines" << std::endl;
 		m_file << "0.5 72 div " << "slw" << endl;	
	} else {
		m_file << "% Y = const grid lines" << std::endl;
 		m_file << "0.25 72 div " << "slw" << endl;	
	}

	for( auto s = lines.begin(); s != lines.end(); s++ ) {
		double v = s->first;
		const DblValues& ords = s->second;

		double x = v;	
		
		for( int i = 0; i < ords.size(); i+= 2 ) {
			double y1 = ords[i]; 
			double y2 = ords[i+1]; 
			if( kind == 0 )
				writeLine(x, y1, x, y2);	
			else
				writeLine(y1, x, y2, x);	
		

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
		lineX[x].push_back(y);
		lineX[x].push_back(y+s);
		lineX[x+s].push_back(y);
		lineX[x+s].push_back(y+s);
		// lineY
		lineY[y].push_back(x);
		lineY[y].push_back(x+s);
		lineY[y+s].push_back(x);
		lineY[y+s].push_back(x+s);
					
	};

	visitTree(rast.m_root, nodeFunc );

	// write PS.
	{
		ExportPS wtoPS(name);
		wtoPS.write(lineX,0);
		wtoPS.write(lineY,1);
	}
	//cout << " X lines: " << lineX.size();
	//cout << " Y lines: " << lineY.size();
}

