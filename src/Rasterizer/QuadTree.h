#ifndef QUAD_ALLOCATOR_H
#define QUAD_ALLOCATOR_H

#include <vector>

template<class T>
class QuadBlockAllocator 
{
	public:
		QuadBlockAllocator()
		: m_activeBlock(0),
		  m_numItemsOut(0)
		{
		}

		~QuadBlockAllocator()
		{
			int num = m_blocks.size();
			for( int i = 0; i < num; i++ )
				delete m_blocks[i];
		}


	public:

		T *allocItem()
		{
			if( m_freeItems.size() == 0 )
			{
				itemBlock *block;
				if( m_activeBlock == m_blocks.size() )
				{
					block = new itemBlock;
					m_blocks.push_back(block);
				} else
					block = m_blocks[m_activeBlock];

				++m_activeBlock;
				block->m_used = true;
				T *items = block->m_items;
				m_freeItems.resize(BlockSize);
				for( int i = 0; i < BlockSize; ++i )
					m_freeItems[i] = &items[i];
			}
			T *item = m_freeItems.back();
			m_freeItems.pop_back();
			++m_numItemsOut;
			return item;
		}

		void reset()
		{
			m_freeItems.clear();
			m_numItemsOut = 0;
			m_activeBlock = 0;

			int n = m_blocks.size();
			for( int b = 0; b < n; b++ )
			{
				itemBlock *block = m_blocks[b];
				if( block->m_used )
				{
					for( int i = 0; i < BlockSize; ++i )
						block->m_items[i].reset();
					block->m_used = false;

				}
			}
		}

		void returnItem( T *item )
		{
			assert(item);
			item->reset();
			m_freeItems.push_back(item);
			--m_numItemsOut;
		}

		int numUsed() const
		{
			return m_numItemsOut;
		}	

	private:
	
		enum {
			BlockSize = 50
		};

		struct itemBlock
		{
			T m_items[BlockSize];
			bool m_used;
		};

		std::vector<itemBlock*> m_blocks;
		std::vector<T *> m_freeItems;
		int	m_activeBlock;
		int	m_numItemsOut;
};

class QuadNode {

	public:
		QuadNode();
		~QuadNode();

		void reset();	
		void setParent(QuadNode*);
		void setLevel(int);
		void setLeaf(bool);
		void setSubQuad(int, QuadNode*);		
		QuadNode* getSubQuad(int);
		const QuadNode* getSubQuad(int) const;
	
		void setSubQuadEdge(int);
		bool subQuadEdge(int) const;	

		QuadNode* getParent();
		bool isLeaf() const;
	public:
		double m_haar[3];
		int 	m_level;
	private:

		QuadNode*	m_quads[4];
		QuadNode*	m_parent;
		int			m_hasEdge;
		bool		m_leaf;
};

class QuadTree
{
	public:
		QuadTree();
		~QuadTree();

	public:
		QuadNode* newQuad( QuadNode* parent = 0);
		void reset();
	public:
		double m_haar00;
	private:
		QuadBlockAllocator<QuadNode> m_quads;		
};
#endif
