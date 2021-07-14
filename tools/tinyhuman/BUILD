# Runs an evaluation of biograph's variant calling capabilities based
# on a ~1% random sampling of high confidence regions of HG001.
#
# This makes it easy to evaluate the effect of changes to the variant
# caller in minutes rather than hours.
#
# It requires a number of utilities to already be installed:
#   vcftools, vct-sort, bgzip, tabix, rtg
#
# To run:  bazel run //tools/tinyhuman:compare-tinyhuman-rand
sh_binary(
    name = "compare-tinyhuman-rand",
    srcs = ["compare-tinyhuman-rand.sh"],
    data = ["//modules/biograph:bgbinary"],
)
