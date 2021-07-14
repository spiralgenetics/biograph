#ifndef __progress_tracker_types_h__
#define __progress_tracker_types_h__

#include <functional> 
#include <cstddef>

// update function called by the zip_* classes.
//
// total_input_read is how many bytes we have read so far
// total_output_written is how many bytes we have written so far
// 
// Returns a <modulo> which can be used by the caller of the update function
// to delay the next call to it until <modulo> bytes have been seen.
//
// This allows the provider of the update function to dynamically redefine the
// speed at which things should be updated.
// 
// It may not make sense to call the update function for every byte of an input
// that is GB in size.
//
typedef const std::function< size_t (size_t total_input_read, size_t total_output_written) > progress_t;

// A trivial update function that does nothing.
// Used a default input parameter.
//
inline size_t no_update(size_t in, size_t out){ return (size_t)1;};

#endif
