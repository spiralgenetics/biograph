
#include "modules/bio_mapred/compute_coverage.h"
#include "modules/io/bitcount.h"
#include "modules/mapred/resource_manager.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/coverage_record.h"
#include "modules/bio_base/seq_position.h"
#include <boost/bind.hpp>

REGISTER_TASK(compute_coverage_task);

void compute_coverage_task::run()
{
	progress_handler_t my_prog = boost::bind(&compute_coverage_task::void_progress, this, _1);
        subprogress sp1(my_prog, 0.0, 0.7);
        subprogress sp2(my_prog, 0.7, 0.75);
        subprogress sp3(my_prog, 0.75, 0.8);
        subprogress sp4(my_prog, 0.8, 1.0);

	::reference ref(reference);
	SPLOG("compute_coverage_task::run> Start");

	size_t one_size = bitcount::compute_size(ref.size());
	mmap_buffer bit_buf;
	resource_manager rm;
	rm.create_resource(bit_buf, one_size*2);

	bitcount bc_uniq(bit_buf.buffer(), ref.size());
	bitcount bc_guess(bit_buf.buffer() + one_size, ref.size());
	bc_uniq.init();
	bc_guess.init();

	manifest_reader mr(input);
	seq_position sp;
	coverage_record cr;
	double tot_recs = input.get_num_records();
	size_t i = 0;
	while(mr.read_msgpack(sp, cr)) {
		sp1(double(i) / tot_recs);
		size_t flat_pos = ref.flatten(sp);
		if (cr.match_count == 1)
			bc_uniq.set(flat_pos, true);
		else
			bc_guess.set(flat_pos, true);
	}
	bc_uniq.finalize(sp2);
	bc_guess.finalize(sp3);

	manifest out;
	rm.write_resource(out, bit_buf, get_root(), "cov", sp4);
	out.metadata().set(meta::ns::readonly, "ref_size", ref.size());
	set_output(out);
	SPLOG("compute_coverage_task::run> Done");
}
