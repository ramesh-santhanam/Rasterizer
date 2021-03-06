#ifndef Rasterer2d_h
#define Rasterer2d_h

#include "QuadTree.h"

class Tree2ps;

// Rasterize 2d polygon one segment at a time.
class Rasterer2d {

	public:
	friend class Tree2ps;

	public:


        Rasterer2d(int depth);
        ~Rasterer2d();

        // rasterize a line
        void rasterize(double x0, double y0, double x1, double y1);

	// render as image 
	void render(std::vector<std::vector<double>>& img);

	// to postscript grey scale.
	void toPostscript(const std::vector<std::vector<double>>& img, const std::string& name);

   public:
        void debug();

   protected:

	void writeImage( std::vector<std::vector<double>>& img,
						 int offset_x, int offset_y,
						 int size,
						 const QuadNode* node,
						 double val );

	void writeImageBlock( std::vector<std::vector<double>>& img,
						int  ox,
						int  oy,
						int  sx,
						int  sy,
						double val );

	void _rasterize(QuadNode* node, double x0, double y0, 
			double x1, double y1, int level);

        QuadNode* getChild(QuadNode* n, int index);

    public:
          int m_depth;
          QuadTree m_qtree;
          QuadNode* m_root;
};
#endif

