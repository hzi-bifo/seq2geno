#### Step 1. check if Conda is installed
#### Step 2. check where seq2geno is installed
It might look like `/YOUR/HOME/bin/seq2geno`
#### Step 3. edit the value of SEQ2GENO and add the lines below in your ~/.profile
```
export SEQ2GENO_HOME=$HOME/bin/seq2geno
export PATH=$(realpath $SEQ2GENO_HOME)/main:$PATH
```
If this package was acquired as Seq2Geno2Pheno was installed, these variable will be set up by Seq2Geno2Pheno.

#### Step 4. install the core environment
```
conda create --file snakemake_env.yml
```
#### Step 5. test if the core environment 
```
source activate snakemake_env
source deactivate
```
Although `source activate` is functionally interchangeable with `conda activate`, we still recommend to test with the former. It will be required by snakemake, so these commands could also help to ensure.  
