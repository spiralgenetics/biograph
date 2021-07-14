
#ifndef __aggregate_map_h__
#define __aggregate_map_h__

#include <stdlib.h>
#include <utility>
#include <boost/iterator/iterator_facade.hpp>

#ifdef AGGREGATE_MAP_DEBUG
#include <stdio.h>
#endif

#include "base/base.h"

template<class Key, class Data>
class aggregate_map
{
private:
	typedef std::pair<Key, Data> value_type;
	struct node_t
	{
		value_type value; // No const because erase does evil things
		Data total;
		node_t* link[2];
		node_t* parent;
		bool red;
		node_t(const value_type& _value)
			: value(_value)
			, total(_value.second)
			, red(true)
		{
			link[0] = NULL;
			link[1] = NULL;
		}
		void recompute()
		{
			total = value.second;
			if (link[0])
				total += link[0]->total;
			if (link[1])
				total += link[1]->total;
		}
		void set_link(int dir, node_t* child)
		{
			link[dir] = child;
			if (child != NULL) child->parent = this;
		}
	};
private:
	bool is_red(node_t* n) { return n != NULL && n->red; }
	node_t* rot_single(node_t* n, int dir)
	{
		node_t *s = n->link[!dir];
		n->set_link(!dir, s->link[dir]);
		s->set_link(dir, n);
		n->red = true;
		s->red = false;
		n->recompute();
		return s;
	}
	node_t *rot_double(node_t* n, int dir)
	{
		n->set_link(!dir, rot_single(n->link[!dir], !dir ));
		return rot_single(n, dir);
	}

	node_t *insert_rec(node_t* n, const value_type& v)
	{
		if (n == NULL)
		{
			n = new node_t(v);
			return n;
		}
		CHECK(n->value.first < v.first || v.first < n->value.first);  // Assert no duplicates

		int dir = n->value.first < v.first;
		n->set_link(dir, insert_rec(n->link[dir], v));
		if (is_red(n->link[dir]))
		{
			if (is_red(n->link[!dir]))
			{
				n->red = true;
				n->link[0]->red = false;
				n->link[1]->red = false;
			}
			else
			{
				if (is_red(n->link[dir]->link[dir]))
					n = rot_single(n, !dir);
				else if (is_red(n->link[dir]->link[!dir]))
					n = rot_double(n, !dir);
			}
		}
		n->recompute();
		return n;
	}

	node_t* erase_rec(node_t* n, bool& done, Key k, bool& removed)
	{
		if (n == NULL)
			return n;  // Ignore missing deletes

		if (n->value.first == k)
		{
			if (n->link[0] == NULL || n->link[1] == NULL)
			{
				node_t* s = n->link[n->link[0] == NULL ? 1 : 0];
				if (is_red(n))
					done = true;
				else if (is_red(s))
				{
					s->red = false;
					done = true;
				}
				delete n;
				removed = true;
				return s;
			}
			node_t* heir = n->link[0];
			while( heir->link[1] != NULL)
				heir = heir->link[1];
			
			n->value = heir->value;
			k = heir->value.first;
		}

		int dir = n->value.first < k;
		n->set_link(dir, erase_rec(n->link[dir], done, k, removed));
		if (!done)
			n = erase_balance(n, dir, done);
		n->recompute();
		return n;
	}

	node_t* erase_balance(node_t* n, int dir, bool& done)
	{
		node_t *p = n;
		node_t *s = n->link[!dir];

		if (is_red(s))
		{
			n = rot_single(n, dir);
			n->recompute();
			s = p->link[!dir];
		}

		if (s != NULL)
		{
			if (!is_red(s->link[0]) && !is_red(s->link[1]))
			{
				if (is_red(p))
					done = true;
				p->red = false;
				s->red = true;
			}
			else
			{
				bool save = p->red;
				bool new_root = (n == p);
				
				if (is_red(s->link[!dir]))
					p = rot_single(p, dir);
				else
					p = rot_double(p, dir);
		
				p->red = save;
				p->link[0]->red = false;
				p->link[1]->red = false;
				p->recompute();
		
				if (new_root)
					n = p;
				else
				{
					n->set_link(dir, p);
					n->recompute();
				}
			
				done = true;
			}
		}
		return n;		
	}

	node_t* find_rec(node_t* n, const Key& k) const
	{
		if (n == NULL)
			return NULL;
		if (n->value.first == k)
			return n;
		if (k < n->value.first)
			return find_rec(n->link[0], k);
		else
			return find_rec(n->link[1], k);
	}

	node_t* lower_bound_rec(node_t* n, const Key& k) const
	{
		if (n->value.first == k)
			return n;
		if (k < n->value.first)
		{
			if (n->link[0] == NULL)
				return n;
			return lower_bound_rec(n->link[0], k);
		}
		else
		{
			if (n->link[1] == NULL)
				return move(n, 1);
			return lower_bound_rec(n->link[1], k);
		}
	}

	Data total_dir_rec(node_t* n, const Key& k, int dir, bool inclusive) const
	{
		if (n == NULL)
			return Data();

		bool key_in_range = (n->value.first == k ? inclusive : int(k < n->value.first) == dir);
		if (!key_in_range)
			return total_dir_rec(n->link[dir], k, dir, inclusive);

		Data tot = n->value.second;
		if (n->link[dir] != NULL)
			tot += n->link[dir]->total;
		tot += total_dir_rec(n->link[!dir], k, dir, inclusive);
		return tot;
	}

	Data total_rec(node_t* n, const Key& begin, const Key end, bool begin_inclusive, bool end_inclusive) const
	{
		// Nothing in range, return empty value
		if (n == NULL)
			return Data();
		
		// If key is outside of range, descend in proper direction
		if (n->value.first < begin || (!begin_inclusive && n->value.first == begin))
			return total_rec(n->link[1], begin, end, begin_inclusive, end_inclusive);
		if (n->value.first > end || (!end_inclusive && n->value.first == end))
			return total_rec(n->link[0], begin, end, begin_inclusive, end_inclusive);

		// If key is in range, do 'directional sums' in both directions
		Data tot = n->value.second;
		tot += total_dir_rec(n->link[0], begin, 1, begin_inclusive);
		tot += total_dir_rec(n->link[1], end, 0, end_inclusive);
		return tot;
	}

	static node_t* move(node_t* n, int dir)
	{
		if (n->link[dir] != NULL)
		{
			n = n->link[dir];
			while(n->link[!dir] != NULL)
				n = n->link[!dir];
			return n;
		}
		while(true)
		{
			if (n->parent == NULL || n == n->parent->link[!dir])
			{
				n = n->parent;
				break;
			}
			n = n->parent;
		}
		return n;
	}

#ifdef AGGREGATE_MAP_DEBUG
	int validate_rec(node_t* n)
	{
		if (n == NULL)
			return 1;
		
		node_t *ln = n->link[0];
		node_t *rn = n->link[1];

		if (is_red(n) &&
			(is_red(ln) || is_red(rn)))
		{
			printf("Red violation");
			return 0;
		}

		if ((ln != NULL && ln->parent != n) ||
		    (rn != NULL && rn->parent != n))
		{
			printf("Parent violation");
			return 0;
		}

		int lh = validate_rec(ln);
		int rh = validate_rec(rn);
	
		if ((ln != NULL && ln->value.first >= n->value.first) ||
		    (rn != NULL && rn->value.first <= n->value.first))
		{
			printf("Binary tree violation");
			return 0;
		}

		if (lh != 0 && rh != 0 && lh != rh)
		{
			printf("Black violation");
			return 0;
		}
		
		Data total = n->value.second;
		if (ln != NULL)
			total += ln->total;
		if (rn != NULL)
			total += rn->total;
		if (total != n->total)
		{
			printf("Total violation");
			return 0;
		}

		if (lh != 0 && rh != 0)
			return is_red(n) ? lh: lh + 1;
		else
			return 0;
	}

	void dump_rec(node_t* n, int depth)
	{
		if (n == NULL)
			return;
		dump_rec(n->link[0], depth+1);
		for(int i = 0; i < depth; i++)
			printf("  ");
		printf("%d, %d, %s\n",
			n->value.first, n->value.second, (n->red ? "red" : "black"));
		dump_rec(n->link[1], depth+1);
	}
#endif

public:
	class const_iterator : public boost::iterator_facade< 
		const_iterator, 
		const value_type, 
		boost::bidirectional_traversal_tag
		>
	{
		friend class aggregate_map;
		friend class boost::iterator_core_access;
	public:
		const_iterator() : m_node(NULL) {}
		
	private:
		const_iterator(node_t* root) : m_root(root), m_node(NULL) {}
		const_iterator(node_t* root, node_t* node) : m_root(root), m_node(node) {}

		void increment() { m_node = move(m_node, 1); }
		void decrement() {
			if (m_node != NULL)
				m_node = move(m_node, 0); 
			else
			{
				m_node = m_root;
				while(m_node && m_node->link[1] != NULL)
					m_node = m_node->link[1];
			}
		}

		bool equal(const_iterator const& other) const
		{
			return m_node == other.m_node;
		}

		const value_type& dereference() const
		{
			return m_node->value;
		}
		node_t* m_root;
		node_t* m_node;
	};
	
	aggregate_map() : m_size(0), m_root(NULL) {}
	~aggregate_map() { clear(); }

	Data total(const Key& begin, const Key& end) const
	{
		CHECK_GE(end, begin);
		return total_rec(m_root, begin, end, true, false);
	}

	Data total(const const_iterator& begin, const const_iterator& end) const
	{
		if (end == const_iterator(m_root))
		{
			if (begin == const_iterator(m_root))
				return Data();
			return total_dir_rec(m_root, begin->first, 1, true);
		}
		return total_rec(m_root, begin->first, end->first, true, false);
	}

	const_iterator begin() const 
	{
		node_t* leftmost = m_root;
		if (leftmost != NULL)
		{
			while(leftmost->link[0] != NULL)
				leftmost = leftmost->link[0];
		}
		return const_iterator(m_root, leftmost);
	}

	const_iterator end() const 
	{
		return const_iterator(m_root);
	}

	const_iterator find(const Key& k) const
	{
		node_t * r = find_rec(m_root, k);
		return const_iterator(m_root, r);
	}

	const_iterator lower_bound(const Key& k) const
	{
		node_t * r = NULL;
		if (m_root != NULL)
			r = lower_bound_rec(m_root, k);
		return const_iterator(m_root, r);
	}

	void insert(const value_type& val)
	{
		m_root = insert_rec(m_root, val);
		m_root->red = false;
		m_root->parent = NULL;
		m_size++;
	}

	// Note, iterators are NOT stable during erase
	void erase(const Key& k)
	{
		bool done = false;
		bool removed = false;
		m_root = erase_rec(m_root, done, k, removed);
		if (removed) m_size--;
		if (m_root != NULL)
		{
			m_root->red = false;
			m_root->parent = NULL;
		}
	}

	void erase(const const_iterator& it)
	{
		erase(it->first);
	}

	void update(const_iterator it, const Data& data)
	{
		node_t* node = it.m_node;
		node->value.second = data;
		while(node)
		{
			node->recompute();
			node = node->parent;
		}
	}
	
	size_t size() const { return m_size; }

private:
	void clear_rec(node_t* n)
	{
		if (n == NULL) return;
		clear_rec(n->link[0]);
		clear_rec(n->link[1]);
		delete n;
	}
public:
	void clear() 
	{
		clear_rec(m_root);
		m_root = NULL;
		m_size = 0;
	}

#ifdef AGGREGATE_MAP_DEBUG
	bool validate()
	{
		int r = validate_rec(m_root);
		return r;
	}

	void dump()
	{
		printf("-------------------------------\n");
		dump_rec(m_root, 0);
	}
#endif

public:
	// Special 'tree' interface for tree walkers (really just interval_tree)
	void* get_root() const { return m_root; }
	void* get_node(const const_iterator& it) const { return it.m_node; }
	const Key& get_key(void* node) const { return ((node_t*) node)->value.first; }
	const Data& get_data(void* node) const { return ((node_t*) node)->value.second; }
	const Data& get_total(void* node) const { return ((node_t*) node)->total; }
	void* get_down(void* node, int dir) const { return ((node_t*) node)->link[dir]; }
	void* get_up(void* node) const { return ((node_t*) node)->parent; }

private:	
	size_t m_size;
	node_t* m_root;
};

#endif

