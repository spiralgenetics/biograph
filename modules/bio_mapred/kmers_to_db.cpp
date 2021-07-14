#include "modules/bio_base/kmer.h"
#include "modules/bio_mapred/kmers_to_db.h"
#include "modules/bio_mapred/align_kmer.h"
#include "modules/bio_base/reference.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/resource_manager.h"

#include <boost/bind.hpp>

REGISTER_TASK(kmers_to_db_task);

void kmers_to_db_task::run()
{
	// Setup progess handling
	progress_handler_t my_prog = boost::bind(&kmers_to_db_task::void_progress, this, _1);

	unsigned kmer_size = input.metadata().get<unsigned>(meta::ns::readonly, "kmer_size");

	// Make kmer set
	subprogress sp1(my_prog, 0.0, (ref_name == "" ? 0.8 : 0.5));
	SPLOG("kmers_to_db_task::run> Loading k-mers into memory");
	manifest_reader mr(input);
	auto num_records = input.get_num_records();
	// TODO: Using kmer_trace only for callback management is lame
	kmer_set kdb(mr, num_records, kmer_size, [&](
		size_t index,
		const kmer_t& kmer,
		size_t kmer_size,
		const std::string& value)
	{
		sp1(double(index) / num_records);
	});

	// Write DB out
	subprogress sp2(my_prog, 0.8, 1.0);
	SPLOG("kmers_to_db_task::run> Writing k-mers out");
	std::string kmer_db = kdb.save(get_root(), sp2);

	manifest out;
	out.metadata().set(meta::ns::internal, "kmer_db", kmer_db);
	set_output(out);
}
