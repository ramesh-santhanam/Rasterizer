#include "Rasterer2d.h"
#include "QuadTree.h"

#include <cstddef>
#include <cassert>
#include <iostream>

using namespace std;

Rasterer2d::Rasterer2d(int depth)
:m_depth(depth)
{
	m_qtree.reset();
	m_root=nullptr;
}

Rasterer2d::~Rasterer2d()
{
	m_qtree.reset();
	m_root = nullptr;
}

void
Rasterer2d::debug()
{
	cout << "coefficient level 0:" <<  m_qtree.m_haar00 << "\n";
	cout << "coeff level::" << m_root->m_level << ":" << m_root->m_haar[0] << "," << m_root->m_haar[1] << "," << m_root->m_haar[2] << "\n";
	cout << "+++++++++++\n";
}

void
Rasterer2d::render()
{

}

void
Rasterer2d::rasterize(double x0, double y0, double x1, double y1)
{
	// eqvt: area of triangle (x0,y0), (x1,y1) && (0,0)
	auto evalC00 = [] (double x0, double y0, double x1, double y1) ->double {
        	return 0.5 * ( x0 * y1 - x1 * y0 );
    	};

	m_qtree.m_haar00 += evalC00( x0, y0, x1, y1 );	
	int level = 0;

	if( m_depth > 0 ) {
		if( !m_root  ) {
			m_root = m_qtree.newQuad(nullptr);
			//cout << "allocated root: %p " << m_root << "\n";
		}
		assert(m_root);
		QuadNode* r = m_root;
		//cout << "root: %p " << r << "\n";

		//cout << "x0: " << x0 << ":y0:" << y0 << ":x1:" << x1 << ":y1:" << y1 << "\n";
		_rasterize( r, x0, y0, x1, y1, level );	
	}
}

QuadNode*
Rasterer2d::getChild(QuadNode* n, int index)
{
	QuadNode* child = n->getSubQuad(index);
	if( ! child ) {
		//cout << "%p(%d)" << n << n->m_level << ":" << "child: " << index << ": is nil\n";
		child = m_qtree.newQuad(n);
		//cout << "%p(%d)" << n << n->m_level << ":child:" << "(%p)" << child << "n";
		n->setSubQuad(index, child);
	}
	assert(child);
	assert(child->getParent() == n);
	return child;
}

void
Rasterer2d::_rasterize(QuadNode* node, double x0, double y0, double x1, double y1, int level )
{
	auto xform = [] (double &px, double &py, int quad ) {

      	switch(quad) {
		case 0:	
       		       	px *= 2.0;
                	py *= 2.0;
                	break;

            	case 1:
	                px = 2*px - 1;
       		        py *= 2.0;
                	break;

		case 2:
                	px *= 2.0;
                	py = 2.*py -1.;
                	break;

            	default:
			assert(quad==3);
       		        px = 2.*px -1.;
                	py = 2.*py -1.;
                	break;
        	}
    	};

	// classify p(x,y) in [0,1] into one of four quadrants of
	// a unit square split using (0.5,0.5).
	//
	auto pointCode = [] (double px, double py ) ->unsigned int {
		unsigned int code = 0;
		if ( px >= 0.5 )
			code |= 0x01;
		if ( py >= 0.5 )
			code |= 0x02;
		assert(code==0||code==1||code==2||code==3);
		return code;
	};

	// basis function C01 - contribution per quadrant.
	auto evalC01 = [] (double x0, double y0, double x1, double y1, int quad)->double {
    		double c = 1./8.;
       		if( quad == 0 || quad == 1)
       			return c * (x0 - x1 ) * (y0 + y1 );
       		assert(quad == 2 || quad ==3);
       		return c * (2. - (y0 + y1)) * (x0-x1);
    	};

    // basis function C10 - contribution per quadrant.
    auto evalC10 = [&] (double x0, double y0, double x1, double y1, int quad)->double {
    	double c = 1./8.;
       	if( quad == 0 || quad == 2)
          		return c * ( y1 - y0 ) * ( x0 + x1 );
       	assert(quad == 1 || quad == 3 );
       	return c * (2. - (x0+x1) ) * ( y1 - y0 );
    };

	// square basis function C11 - contrbution per quadrant
    auto evalC11 = [] (double x0, double y0, double x1, double y1, int quad)->double {
       	double c = 1./8.;
       	if( quad == 0 )
       		return c * (x0 + x1) * (y1 - y0);
       	if( quad == 1 )
       		return c * (2. - (x0 + x1)) * (y1 - y0);
       	if( quad == 2 )
      		return c * (x0 + x1) * (y0 - y1);
       	assert(quad == 3);
       	return c * (2. - (x0 + x1)) * (y0 - y1);
    };

    // evaluate basis function contribution for line segment.
    auto evalCoefficients = [&] (double *coeff, double px, double py, double qx, double qy, int quad ) ->void {
    	coeff[0] += evalC10(px, py, qx, qy, quad);
       	coeff[1] += evalC01(px, py, qx, qy, quad);
       	coeff[2] += evalC11(px, py, qx, qy, quad);
    };

    auto swap = [] (double &t0, double &t1) {
    	if( t0 > t1 ) {
        	double t = t1;
           	t1 = t0;
           	t0 = t;
       	}
       	assert(t1 >= t0);
    };

	auto interpolate = [] (double px, double qx, double t ) -> double {
		return px + t * (qx - px );
	};

	node->setLevel(level);
	if( level == m_depth )
		node->setLeaf(true);

	unsigned int code0 = pointCode(x0, y0);
	unsigned int code1 = pointCode(x1, y1);

	if( code0 == code1 ) {
		xform( x0, y0, code0 );	
		xform( x1, y1, code1 );	
		evalCoefficients( node->m_haar, x0, y0, x1, y1, code0 );

		node->setSubQuadEdge(code0);

		if( level == m_depth ) 
			return;

		node = getChild(node, code0);
		return _rasterize( node, x0, y0, x1, y1, level+1);
	}

	unsigned int v = code0^code1;
	switch(v) {
		case 1:
            {
				// split in x.
   	            double t = (0.5 - x0 )/ (x1 - x0);
               	double rx[2], ry[2];
               	rx[0] = 0.5;
               	ry[0] = y0 + t * (y1-y0);

               	rx[1] = 0.5;
               	ry[1] = ry[0];

               	xform( x0, y0, code0 );
               	xform( rx[0], ry[0], code0 );
                evalCoefficients( node->m_haar, x0, y0, rx[0], ry[0], code0 );
				node->setSubQuadEdge(code0);
				if( level != m_depth ) {
					QuadNode* c = getChild(node, code0);
					_rasterize( c, x0, y0, rx[0], ry[0], level+1);
				}
                xform( x1, y1, code1 );
                xform( rx[1], ry[1], code1 );
                evalCoefficients( node->m_haar, rx[1], ry[1], x1, y1, code1 );
				node->setSubQuadEdge(code1);
				if( level != m_depth ) { 
					QuadNode* c = getChild(node, code1);
					_rasterize( c, rx[1], ry[1], x1, y1, level+1);
				} 
               	break;
			}
		case 2:
            {
           		// split in y.
           		double t = (0.5 - y0 )/ (y1 - y0);
               	double rx[2], ry[2];
               	rx[0] = x0 + t * (x1-x0);
               	ry[0] = 0.5;
               	rx[1] = rx[0];
               	ry[1] = ry[0];
               	xform( x0, y0, code0);
               	xform( rx[0], ry[0], code0);
                evalCoefficients( node->m_haar, x0, y0, rx[0], ry[0], code0 );
				node->setSubQuadEdge(code0);
				if( level != m_depth ) {
					QuadNode* c = getChild(node, code0);
					_rasterize( c, x0, y0, rx[0], ry[0], level+1);
				}

                xform( x1, y1, code1);
                xform( rx[1], ry[1], code1);
                evalCoefficients( node->m_haar, rx[1], ry[1], x1, y1, code1 );
				node->setSubQuadEdge(code1);
				if( level != m_depth ) {
					QuadNode* c = getChild(node, code1);
					_rasterize( c, rx[1], ry[1], x1, y1, level+1);
				}
                break;
           }	

		default:
		{
			assert(v == 3);

            double rx[4], ry[4];
            // split in x
            double t0 = (0.5 - x0 )/ (x1 - x0);
            // split in y.
            double t1 = (0.5 - y0 )/ (y1 - y0);

            swap(t0, t1);
            rx[0] = interpolate( x0, x1, t0 );
            ry[0] = interpolate( y0, y1, t0 );
            rx[1] = rx[0];
            ry[1] = ry[0];

            rx[2] = interpolate(x0, x1, t1 );
            ry[2] = interpolate(y0, y1, t1 );
            rx[3] = rx[2];
            ry[3] = ry[2];

            xform( x0, y0, code0);
            xform( rx[0], ry[0], code0);
            evalCoefficients( node->m_haar, x0, y0, rx[0], ry[0], code0 );
			node->setSubQuadEdge(code0);
			if( level != m_depth ) {
				QuadNode* c = getChild(node, code0);
				_rasterize( c, x0, y0, rx[0], ry[0], level+1);
			}	
 			if( t1 > t0 ) {
                unsigned int code = pointCode(rx[1], ry[1]);
                xform( rx[1], ry[1], code);
                xform( rx[2], ry[2], code);
                evalCoefficients( node->m_haar, rx[1], ry[1], rx[2], ry[2], code );
				node->setSubQuadEdge(code);
				if( level != m_depth ) {
					QuadNode* c = getChild(node, code);
					_rasterize( c, rx[1], ry[1], rx[2], ry[2], level+1);
				}	
            }
            xform( rx[3], ry[3], code1);
            xform( x1, y1, code1);
            evalCoefficients( node->m_haar, rx[3], ry[3], x1, y1, code1 );		
			node->setSubQuadEdge(code1);
			if( level != m_depth ) {
				QuadNode* c = getChild(node, code1);
				_rasterize( c, rx[3], ry[3], x1, y1, level+1);
			}
			break;
		}
	} // switch	
	
}
