#### Step 1. check if Conda is installed
#### Step 2. check where seq2geno is installed
It might look like `/YOUR/HOME/bin/seq2geno`
#### Step 3. edit the environment variables and add the lines below in your ~/.profile
```
export SEQ2GENO_HOME=/YOUR/HOME/bin/seq2geno
export PATH=$( realpath $SEQ2GENO_HOME )/main:$PATH
```

#### Step 4. install the core environment
create the environment with the commands:
```
conda env create -n snakemake_env --file=snakemake_env.yml
```
We use the name "snakemake_env", but
this might need changes if any previous environment is already named after it
#### Step 5. test if the core environment 
```
source activate snakemake_env
seq2geno -h
source deactivate
```
Your conda might ask you to replace `source activate` with `conda activate` 

#### Step 6. install dependencies of Raory 
Roary already has its own script for installing dependencies, so we can simply use it:

```
source activate snakemake_env
cd $SEQ2GENO_HOME/denovo/lib/Roary
./install_dependencies.sh
source deactivate
```
These softwares should now be available under `$SEQ2GENO_HOME/denovo/lib/Roary/build`.
