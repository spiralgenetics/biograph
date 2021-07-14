#pragma once

#include <memory>

// TODO: When we upgrade to C++14, switch to std::make_unique.
template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
	return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}
