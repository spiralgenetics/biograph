#include "modules/test/build_ref.h"
#include "modules/pipeline/primitives.h"
#include "modules/bio_format/make_vars.h"
#include <gtest/gtest.h>

TEST(make_vars, runs_without_exceptions)
{
	perform_build_ref("yeast", "datasets/fasta/saccharomyces_cerevisiae_EF4.fasta");

	make_vars mv("yeast", 400, 100, "yout");
	mv.snp("snp");
	mv.random_insert("ins", 512);
	mv.repeat_insert("rep", 512);
	mv.random_delete("del", 512);
	mv.transpose("trn", 512);
}
