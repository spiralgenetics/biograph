#!/bin/bash


#rm /home/ubuntu/notebook/results/candidates_test_two_parent.*
#rm /home/ubuntu/notebook/results/candidates_test_single_parent.*
# test against NA12878
python kmer_analysis.py --ref "/reference/human_g1k_v37" --proband "/biographs/juvpsych/0847-01_sorted_new" --variant_file "/home/ubuntu/notebook/tests/test_variant_types.vcf" --parent "/biographs/juvpsych/0847-02_map_new" --parent "/biographs/juvpsych/0847-03_sorted_new" --output "/home/ubuntu/customer_integration/output/regression/candidates_test_two_parent" > /dev/null

if [ "$1" == "nosnps" ]; then
	diff /home/ubuntu/notebook/results/candidates_base_two_parent_nosnps.vcf /home/ubuntu/customer_integration/output/regression/candidates_test_two_parent.vcf > /dev/null
	echo "two parent vcf = "$?
	diff /home/ubuntu/notebook/results/candidates_base_two_parent_nosnps.csv /home/ubuntu/customer_integration/output/regression/candidates_test_two_parent.csv > /dev/null
	echo "two parent csv = "$?
else
	diff /home/ubuntu/notebook/results/candidates_base_two_parent.vcf /home/ubuntu/customer_integration/output/regression/candidates_test_two_parent.vcf > /dev/null
	echo "two parent vcf = "$?
	diff /home/ubuntu/notebook/results/candidates_base_two_parent.csv /home/ubuntu/customer_integration/output/regression/candidates_test_two_parent.csv > /dev/null
	echo "two parent csv = "$?
fi

python kmer_analysis.py --ref "/reference/human_g1k_v37" --proband "/biographs/juvpsych/0847-01_sorted_new" --variant_file "/home/ubuntu/notebook/tests/test_variant_types.vcf" --parent "/biographs/juvpsych/0847-02_map_new" --output "/home/ubuntu/customer_integration/output/regression/candidates_test_single_parent" > /dev/null

if [ "$1" == "nosnps" ]; then
	diff /home/ubuntu/notebook/results/candidates_base_single_parent_nosnps.vcf /home/ubuntu/customer_integration/output/regression/candidates_test_single_parent.vcf > /dev/null
	echo "two parent vcf = "$?
	diff /home/ubuntu/notebook/results/candidates_base_single_parent_nosnps.csv /home/ubuntu/customer_integration/output/regression/candidates_test_single_parent.csv > /dev/null
	echo "two parent csv = "$?
else
	diff /home/ubuntu/notebook/results/candidates_base_single_parent.vcf /home/ubuntu/customer_integration/output/regression/candidates_test_single_parent.vcf > /dev/null
	echo "two parent vcf = "$?
	diff /home/ubuntu/notebook/results/candidates_base_single_parent.csv /home/ubuntu/customer_integration/output/regression/candidates_test_single_parent.csv > /dev/null
	echo "two parent csv = "$?
fi

echo "Test done."
