
#ifndef __inverse_kvread_h__
#define __inverse_kvread_h__

#include "modules/io/io.h"
#include "modules/io/loop_io.h"
#include "modules/io/keyvalue.h"

#include <deque>

class inverse_kvread : public reset_readable
{
public:
	inverse_kvread(reset_kv_source& source);
	size_t read(char* buf, size_t len) override;
	void reset() override;

private:
	reset_kv_source& m_source;
	loop_io m_loop; 
	kv_writer m_loop_write;
};

#endif 
