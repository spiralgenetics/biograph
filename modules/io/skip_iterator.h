
#ifndef __skip_iterator_h__
#define __skip_iterator_h__

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/transform_iterator.hpp>

template<class It>
class skip_iterator : public boost::iterator_facade<
        skip_iterator<It>,
        typename std::iterator_traits<It>::value_type,
        std::random_access_iterator_tag,
        typename std::iterator_traits<It>::reference,
        typename std::iterator_traits<It>::difference_type>
{
        friend class boost::iterator_core_access;
public:
        skip_iterator() : m_mult(0) {}
        skip_iterator(const It& it, int mult) : m_it(it), m_mult(mult) {}

private:
        typename std::iterator_traits<It>::reference dereference() const
                { return *m_it; }
        bool equal(const skip_iterator& rhs) const
                { return m_it == rhs.m_it; }
        void increment() { m_it += m_mult; }
        void decrement() { m_it -= m_mult; }
        void advance(typename std::iterator_traits<It>::difference_type n)
                { m_it += n * m_mult; }
        typename std::iterator_traits<It>::difference_type distance_to(const skip_iterator& rhs) const
                { return (rhs.m_it - m_it) / m_mult; }
        It m_it;
        int m_mult;
};

#endif

