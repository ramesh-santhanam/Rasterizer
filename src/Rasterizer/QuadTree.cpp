#include "QuadTree.h"
#include <cstddef>
#include <assert.h>

QuadTree::QuadTree()
{
	m_haar00 = 0;
}

QuadTree::~QuadTree()
{
	m_quads.reset();
}

QuadNode*
QuadTree::newQuad( QuadNode* parent )
{
	QuadNode* qn = m_quads.allocItem();
	if( parent )
		qn->setParent(parent);
	return qn;
}

void
QuadTree::reset()
{
	m_quads.reset();
}

////////////////////////////

QuadNode::QuadNode()
{
	reset();
}

QuadNode::~QuadNode()
{
}

void
QuadNode::reset()
{
	m_parent = nullptr;
	m_quads[0] = nullptr;
	m_quads[1] = nullptr;
	m_quads[2] = nullptr;
	m_quads[3] = nullptr;
	m_level = -1;
	m_leaf = false;
	m_haar[0] = 0.0;
	m_haar[1] = 0.0;
	m_haar[2] = 0.0;
	m_hasEdge = 0;
}

void
QuadNode::setParent(QuadNode* p)
{
	assert(p);
	m_parent = p;
}

void
QuadNode::setLevel(int l)
{
	m_level = l;
}

void
QuadNode::setLeaf(bool isLeaf)
{
	m_leaf = isLeaf;
}

bool
QuadNode::isLeaf() const
{
	return m_leaf;
}

void
QuadNode::setSubQuad(int index, QuadNode* q)
{
	assert(index >=0 && index < 4);
	m_quads[index] = q;
}

QuadNode*
QuadNode::getSubQuad(int q)
{
	assert(q >=0 && q < 4);
	return m_quads[q];
}

const QuadNode*
QuadNode::getSubQuad(int q) const
{
	assert(q >=0 && q < 4);
	return m_quads[q];
}

QuadNode*
QuadNode::getParent()
{
	return m_parent;
}

void
QuadNode::setSubQuadEdge( int q )
{
	switch( q ) {
		case 0:
			m_hasEdge |=  0x01;
			break;
		case 1:
			m_hasEdge |= 0x02;
			break;

		case 2:
			m_hasEdge |= 0x03;
			break;

		case 3:
			m_hasEdge |= 0x04;
			break;
	}	
}

bool
QuadNode::subQuadEdge(int q) const 
{
	if ( q == 0 )
		return m_hasEdge&0x01;
	if ( q == 1 )
		return m_hasEdge&0x02;
	if ( q == 2 )
		return m_hasEdge&0x03;
	assert(q==3);
	return m_hasEdge&0x04;
}
