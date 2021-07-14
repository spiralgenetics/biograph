#include "modules/test/build_ref.h"
#include "modules/pipeline/primitives.h"
#include <gtest/gtest.h>

class build_ref : public ::testing::Test
{
public:
	build_ref() {}

protected:
	void SetUp() override {
		impl.setup();
	}

protected:
	build_ref_impl impl;
};

TEST_F(build_ref, yeast)
{
	impl.run_task("yeast", "datasets/fasta/saccharomyces_cerevisiae_EF4.fasta");
}

TEST_F(build_ref, e_coli)
{
	impl.run_task("e_coli", "datasets/fasta/e_coli_k12.ASM584v1.fasta");
}

