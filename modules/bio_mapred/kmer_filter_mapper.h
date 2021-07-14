
#ifndef __kmer_filter_mapper_h__
#define __kmer_filter_mapper_h__

#include "modules/io/transfer_object.h"
#include "modules/mapred/mapper.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/kmer.h"
#include "modules/bio_base/overrep.h"

struct kmer_filter_params
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(min_count, TF_STRICT);
		FIELD(sys_err_thresh, float(0.0));
		FIELD(rnd_err_thresh, float(0.0));
		FIELD(kmer_size, TF_STRICT);
		FIELD(overrep);
	}
	size_t min_count;
	float sys_err_thresh = 0.0; // .1 
	float rnd_err_thresh = 0.0; // .005 
	size_t kmer_size;
	manifest overrep;
	void validate();
};

class kmer_filter_mapper : public typed_mapper<kmer_filter_mapper,
	kmer_t, kcount_pair,
	kmer_t, kcount_pair>
{
public:
	kmer_filter_mapper(const std::string& params);
	void set_watchdog(const std::function<void()>& watchdog) override { m_watchdog = watchdog; }
        void setup() override;

	~kmer_filter_mapper() {}

	void typed_map(kmer_t key, const kcount_pair& value);

private:
	kmer_filter_params m_params;
	std::function<void()> m_watchdog;
	overrep_map m_overrep;
};

#endif
