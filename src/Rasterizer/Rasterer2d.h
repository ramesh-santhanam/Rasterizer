#ifndef Rasterer2d_h
#define Rasterer2d_h

#include "QuadTree.h"

// Rasterize 2d polygon one segment at a time.
class Rasterer2d {

	public:

        Rasterer2d(int depth);
        ~Rasterer2d();

        // rasterize a line
        void rasterize(double x0, double y0, double x1, double y1);

		// render as image 
		void render(std::vector<std::vector<double>>& img);

   public:
        void debug();

   protected:

		void writeImage( std::vector<std::vector<double>>& img,
						 int offset_x, int offset_y,
						 int size,
						 QuadNode* node, double val );

		void _rasterize(QuadNode* node, double x0, double y0, 
						double x1, double y1, int level);
        QuadNode* getChild(QuadNode* n, int index);

    public:
          int m_depth;
          QuadTree m_qtree;
          QuadNode* m_root;
};
#endif

