#### Step 1 : Identify variants #### 

/scratch/biograph-5.0.0-rc6/biograph variants --in benchmark/${i}_lambda.bg --ref toy_reference --out ${i}_variants.vcf --assemblies-out ${i}_assemblies.csv

#### Step 2 : Run pcmp #### 

bgzip ${i}_variants.vcf
tabix -p vcf ${i}_variants.vcf.gz

bgtools pcmp -b benchmark/${i}_lambda.bg -v ${i}_variants.vcf.gz -r toy_reference/ -o ${i}.pcmp.vcf

bgzip ${i}.pcmp.vcf
tabix -p vcf ${i}.pcmp.vcf.gz

#### Step 3 : Run truvari and rtg #### 
base='/scratch/lambda/lambdaToyData/ml_lambda/lambda_base.vcf.gz'
comp=${i}'.pcmp.vcf.gz'
truvari --reference /scratch/reference/source.fasta --base $base --comp $comp --out disc_truvari --multimatch --passonly
/Project/rtg-tools-3.10.1/rtg vcfeval -b $base -c $comp -t ../ml_lambda/lambda.sdf -o disc_rtg

# bgzip and tabix truvari output vcfs.
cd disc_truvari
for j in *.vcf; do bgzip $j; tabix ${j}.gz; done
cd ..

#### Step 4 : Build ml_tables #### 

#python3 /scratch/lambda/scripts/bgvar2table.py -a ${i}_assemblies.csv -v ${i}.pcmp.vcf.gz  -o ${i}.bench.jl  --mode bench  --trudir disc_truvari --rtgdir disc_rtg

python3 /Project/random_forest_July16/scripts/bgvar2table_sm.py -a ${i}_assemblies.csv -v ${i}.pcmp.vcf.gz  -o ${i}.bench.mod.jl  --mode bench --trudir disc_truvari --rtgdir disc_rtg

####  Step 5 : Run classifier #### 

# For SV's
/Project/random_forest_July16/scripts/run_classifier_rfc_july24.py -c ${i}.bench.mod.jl -v ${i}.pcmp.vcf -m /Project/random_forest_July16/models/model_sv_23july2019.ml -o ${i}.filter.rfc.sv.vcf  -l 0.335 --svonly

# For snps and indels
/Project/random_forest_July16/scripts/run_classifier_rfc_july24.py -c ${i}.bench.mod.jl -v ${i}.pcmp.vcf -m /Project/random_forest_July16/models/model_ao_23july2019.ml -o ${i}.filter.rfc.ao.vcf  -l 0.25 

####  Step 6 : Test filtered results Truvari/RTG #### 
bgzip ${i}.filter.rfc.sv.vcf 
tabix -p vcf ${i}.filter.rfc.sv.vcf.gz

base='/scratch/lambda/lambdaToyData/ml_lambda/lambda_base.vcf.gz'
sdf='/scratch/lambda/lambdaToyData/ml_lambda/lambda.sdf'
comp=${i}'.filter.rfc.sv.vcf.gz'

truvari --reference /scratch/reference/source.fasta --base $base --comp $comp --out ml_truvari --multimatch --passonly


bgzip ${i}.filter.rfc.ao.vcf 
tabix -p vcf ${i}.filter.rfc.ao.vcf.gz
comp=${i}'.filter.rfc.ao.vcf.gz'
/Project/rtg-tools-3.10.1/rtg vcfeval -b $base -c $comp -t $sdf -o ml_rtg
