#include <algorithm>
#include <vector>
template<class Iter, class T>
Iter binary_find(Iter begin, Iter end, T val)
{
    // Finds the lower bound in at most log(last - first) + 1 comparisons
	Iter i = std::lower_bound(begin, end, val);

	if (i != end && !(val < *i))
		return i;
	else
		return end; 
}

//Return the index of an element in sorted vector
template<class Iter, class T> 
int binary_find_index(Iter begin, Iter end, T val)
{
	auto i = binary_find(begin, end, val);
	//think this may need an i+1 if i==end...
  	return std::distance(begin, i);
}
