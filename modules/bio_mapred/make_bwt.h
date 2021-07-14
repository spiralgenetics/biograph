#pragma once

#include "modules/bio_base/bwt_file.h"
#include "modules/bio_base/flat_ref.h"
#include "modules/io/file_io.h"
#include "modules/io/mmap_vector.h"
#include "modules/io/progress.h"
#include "modules/io/version.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/task.h"

struct bwt_flyweight {
		bwt_flyweight() = default;

		bwt_flyweight(uint32_t extent_, uint32_t offset_)
				: extent(extent_)
				, offset(offset_)
		{}

		uint32_t extent;		// supercontig ID
		uint32_t offset;		// offset within supercontig
};

// Takes a .ref and builds a .bwt
class make_bwt_task : public task_impl<make_bwt_task>
{
public:
	static std::string s_type() { return "make_bwt_task"; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input_ref, TF_STRICT);
		FIELD(output_bwt, TF_STRICT);
		FIELD(cent_mod, TF_STRICT);
	}

	void run();
	void validate();

	task_requirements get_requirements() override
	{
		return task_requirements {
			.profile = "himem",
			.cpu_minutes = 60,
		};
	}

	std::string input_ref;
	std::string output_bwt;
	size_t cent_mod = 64;
private:
	std::unique_ptr<flat_ref> m_ref;
};
