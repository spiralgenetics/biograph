#pragma once

#include "modules/io/io.h"
#include <seqan/bam_io.h>

class seqan_readable_adapter
{
public:
	seqan_readable_adapter(readable & source)
		: m_source(source)
    {}
	
	bool at_eof() const 
    { 
        return m_at_eof; 
    }

	size_t get_position() const 
    { 
        return m_position; 
    }
	
	size_t read(char* buf, size_t len)
	{
		if (len == 0) {
            return 0;
        }
		size_t amount_read = m_source.read(buf, len);
		if (amount_read < len) {
            m_at_eof = true;
        }
		m_position += amount_read;
		return amount_read;
	}

private:
	readable& m_source;
	size_t m_position = 0;
	bool m_at_eof = false;
};


namespace seqan {

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<seqan_readable_adapter &, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<seqan_readable_adapter &, IsOutput>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<seqan_readable_adapter &, HasPeek>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<seqan_readable_adapter &, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature<seqan_readable_adapter &, Seek<TSpec> >
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<seqan_readable_adapter &, Tell>
{
    typedef True Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function streamPeek()
// ----------------------------------------------------------------------------

// Not implemented.

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, seqan_readable_adapter & stream)
{
    size_t bytes_read = stream.read(&c, 1);
    if (bytes_read == 0) {
        return EOF;
    }
    
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(seqan_readable_adapter & stream)
{
    return stream.at_eof();
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(seqan_readable_adapter & stream)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
streamReadBlock(char * target, seqan_readable_adapter & stream, size_t maxLen)
{
    return stream.read(target, maxLen);
}


// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

inline int
streamFlush(seqan_readable_adapter & stream)
{
    // Does nothing.
    return 0;
}


// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline size_t
streamTell(seqan_readable_adapter & stream)
{
    return stream.get_position();
}


}  // namespace seqan
