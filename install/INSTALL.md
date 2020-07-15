#### Step 1. check if Conda is installed
#### Step 2. check where seq2geno is installed
It might look like `/YOUR/HOME/bin/seq2geno`
#### Step 3. edit the environment variables and add the lines below in your ~/.profile
```
export SEQ2GENO_HOME=$HOME/bin/
export PATH=$( realpath $SEQ2GENO_HOME )/main:$PATH
```
If this package was acquired as Seq2Geno2Pheno was installed, these variable will be set up by Seq2Geno2Pheno.

#### Step 4. install the core environment
```
conda env create --file=snakemake_env.yml
```
#### Step 5. test if the core environment 
```
source activate snakemake_env
source deactivate
```
Although `source activate` is functionally interchangeable with `conda activate`, the former will be required by snakemake so we still recommend to test with it. 

#### Step 6. install dependencies of Raory 
Roary already has its own script for installing dependencies, so we can simply use it:

```
source activate snakemake_env
cd $SEQ2GENO_HOME/denovo/lib/Roary
./install_dependencies.sh
source deactivate
```
These softwares should now be available under `$SEQ2GENO_HOME/denovo/lib/Roary/build`.
