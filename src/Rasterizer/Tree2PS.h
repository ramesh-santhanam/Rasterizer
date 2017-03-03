#ifndef Tree2ps_h
#define Tree2ps_h

#include <string>

class Rasterer2d;

class Tree2PS {
public:
        Tree2PS();
        ~Tree2PS();
public:
	void write(Rasterer2d& rast, std::string& name);
};
#endif

