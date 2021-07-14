#pragma once

#include <vector>
#include <cassert>

#include "base/base.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/io.h"

// This class can be returned by mmap_buffer and represents a very limited implementation
// of std::vector that operates on a mmap_buffer.
//
// Only the bare minimum needed for our use is implemented, but feel free to add more
// methods as needed...
//
// NOTE:  This "vector" does NOT reallocate.  You can't grow it past its original capacity!

template <typename T>
class mmap_vector
{
public:
	typedef typename std::vector<T>::reference reference;
	typedef typename std::vector<T>::const_reference const_reference;
	typedef typename std::vector<T>::size_type size_type;
	typedef typename std::vector<T>::value_type value_type;
	
public:
	size_type size() const;
	
	reference operator[] (size_type n);
	const_reference operator[] (size_type n) const;
	
	void push_back(const value_type& value);
	void push_back(value_type&& rvalue);
	
	// Use this constructor when you want a new mmap.  Be sure and pass it
	// the capacity as the number of elements you want in the vector, not the
	// mmap size in bytes.
	mmap_vector(size_type vector_capacity)
		: m_size(0)
		, m_capacity(vector_capacity)
	{}	

	void resize(size_type size) { CHECK_LE(size, m_capacity); m_size = size; }
	size_t capacity() { return m_capacity; }

	mmap_buffer& get_buffer() { return m_mmap_buffer; }
	size_t buffer_size() { return m_capacity * sizeof(T); }

	T* begin() { return &(*this)[0]; }
	T* end() { return &(*this)[m_size]; }
	const T* begin() const { return &(*this)[0]; }
	const T* end() const { return &(*this)[m_size]; }
	
	void sync() { m_mmap_buffer.sync(); }
	void reserve(size_type n) { throw io_exception("mmap_vector reserve called, but you may not change its capacity!"); }

private:
	mmap_buffer m_mmap_buffer;
	size_type m_size;
	const size_type m_capacity;
};


template <typename T>
inline typename mmap_vector<T>::size_type mmap_vector<T>::size() const
{
	return m_size;
}


template <typename T>
inline typename mmap_vector<T>::reference mmap_vector<T>::operator[] (mmap_vector<T>::size_type n)
{
	return *(reinterpret_cast<T*>(m_mmap_buffer.buffer()) + n);
}


template <typename T>
inline typename mmap_vector<T>::const_reference mmap_vector<T>::operator[] (mmap_vector::size_type n) const
{
	return *(reinterpret_cast<const T*>(m_mmap_buffer.buffer()) + n);
}


template <typename T>
inline void mmap_vector<T>::push_back(const mmap_vector<T>::value_type& value)
{
	CHECK_LT(m_size, m_capacity);
	operator[](m_size++) = value;
}


template <typename T>
inline void mmap_vector<T>::push_back(mmap_vector<T>::value_type&& rvalue)
{
	CHECK_LT(m_size, m_capacity);
	std::swap(operator[](m_size++), rvalue);
}
