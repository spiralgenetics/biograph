
#include "modules/mapred/manifest.h"
#include "modules/mapred/task.h"
#include "modules/io/progress.h"
#include "modules/io/mmap_vector.h"
#include "modules/io/version.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/reference.h"
#include "modules/io/file_io.h"
#include "modules/io/track_mem.h"

struct flyweight {
        flyweight() : non_ref(false), empty(false) {}
        flyweight(uint64_t start_, uint16_t length_, bool flipped_, bool non_ref_ = false, bool empty_ = false)
                : start(start_)
                , length(length_)
                , flipped(flipped_)
                , non_ref(non_ref_)
                , empty(empty_)
        {}

        flyweight rev_comp() {
                if (flipped) {
                        return flyweight(start - (length - 1), length, false, non_ref, empty);
                } else {
                        return flyweight(start + (length - 1), length, true, non_ref, empty);
                }
        }

        bool valid() { return (! non_ref) && (! empty); }

        uint64_t start : 48;
        unsigned int length  : 13;
        unsigned int flipped : 1; // True when flyweight is in the reverse direction
        unsigned int non_ref : 1; // True when flyweight is a place holder for a non-ref read.
        unsigned int empty   : 1; // True when flyweight is place holder for non-existent mate pair.
};

class mem_seqset_task : public task_impl<mem_seqset_task>
{
public:
	static std::string s_type() { return "mem_seqset_task"; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(read_size, TF_STRICT);
		FIELD(input, TF_STRICT);
		FIELD(coverage, TF_STRICT);
		FIELD(num_threads, TF_STRICT);
		FIELD(max_mem, TF_STRICT);
		FIELD(is_paired, TF_STRICT);
		FIELD(run_tests, TF_STRICT);
		FIELD(ref_name, TF_ALLOW_NULL);
		FIELD(m_write_flat, TF_ALLOW_NULL);
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
	manifest build_seqset();

	std::string ref_name;
	manifest input;
	manifest coverage;
	size_t num_threads = 0;
	size_t max_mem = 0;
	uint8_t read_size = 0;
	bool is_paired = true;
	bool run_tests = false;
	bool m_write_flat = false; 
	
	size_t flywt_index(size_t record_id) const
	{
		return record_id * (is_paired ? 4 : 2);
	}

	static size_t bases_to_data_size(size_t base_count) { return (base_count + 3) / 4; }

private:
	size_t m_max_buf_size;
	mmap_buffer m_repo_mmap;
	const char* m_repo;
	// m_originals contains the flyweights in this format:
	//    m_originals[n]     == normal flyweight
	//    m_originals[n + 1] == reverse complement flyweight
	// if paired:
	//    m_originals[n + 2] == mate flyweight
	//    m_originals[n + 3] == mate reverse complement flyweight
	std::shared_ptr<mmap_vector<flyweight>> m_originals;
	flyweight m_worst_ever;
	std::atomic<uint64_t> m_next_read;
	size_t m_base_pos[5];

	size_t expand_one_read(tracked_vector<flyweight>& output, flyweight read);
	std::string one_expand_pass(const progress_handler_t& progress);
	size_t do_merge(file_writer& output, std::vector<file_reader>& inputs);
	void load_repo(const reference& ref, const progress_handler_t& progress);
	manifest output_seqset(const progress_handler_t& progress, const std::string& output_name, size_t tot_size, bool write_flat);
	manifest do_mem_seqset(const progress_handler_t& progress);

	size_t count_flyweights() const;
};

