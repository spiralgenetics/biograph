
#pragma once

#include "modules/mapred/reducer.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/resource_manager.h"
#include "modules/io/transfer_object.h"
#include "modules/bio_base/seq_position.h"
#include "modules/bio_base/struct_var.h"
#include "modules/bio_base/sv_call.h"
#include "modules/io/bitcount.h"
#include "modules/bio_base/reference.h"

struct sv_call_params
{
	std::string reference;
	manifest coverage;
	size_t read_size;
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(reference);
		FIELD(coverage);
		FIELD(read_size);
	}
};

class sv_call_reducer : public typed_reducer<sv_call_reducer, seq_position, struct_var, seq_position, sv_call>
{
public:
	sv_call_reducer(const std::string& params);

	void set_watchdog(const std::function<void()>& watchdog) override { m_watchdog = watchdog; }
	void typed_start(const seq_position& key);
	void typed_add_value(const seq_position& key, const struct_var& value);
	void typed_end();

private:
	void dump_current();

	sv_call_params m_params;
	reference m_ref;
	mmap_buffer m_buf;
	std::unique_ptr<bitcount> m_bc_uniq;
	std::unique_ptr<bitcount> m_bc_guess;
	seq_position m_pos_start;
	seq_position m_pos_end;
	sv_call m_current;
	std::function<void()> m_watchdog;
};
