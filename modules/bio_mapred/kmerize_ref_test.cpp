#include "modules/test/test_utils.h"
#include "modules/test/local_context.h"
#include "modules/bio_mapred/kmerize_ref.h"
#include "modules/mapred/histogram_export.h"
#include "gtest/gtest.h"

TEST(kmerize_reference, DISABLED_basic)
{
	task_mgr_local ltm;

	std::unique_ptr<kmerize_ref_task> t = make_unique<kmerize_ref_task>();

	json_deserialize(t->params, R"|(
			{
				"kmer_size" : 20,
				"reference" : "e_coli_k12_ASM584v1"
			}
		)|"
	);
	t->params.validate();

	manifest kmers;
	ltm.run_task(kmers, path(make_path("kmerize/bits")), std::move(t));

	std::unique_ptr<map_reduce_task> t2 = make_unique<map_reduce_task>();
        t2->input = kmers;
        t2->map = "value_count";
        t2->sort = "uint64";
        t2->reduce = "sum";
        t2->is_summary = true;
        t2->use_sort = true;
	manifest histogram;
	ltm.run_task(histogram, path(make_path("kmerize/bits")), std::move(t2));

	auto out_path = make_path("kmerize/results");
	simple_export<histogram_exporter>(out_path, histogram);

	ASSERT_TRUE(diff(out_path, "golden/kmerize.txt"));
}

