'''
Purpose:
    Automatically create symbolic cpoies of the data and generate the input
    markdown of Geno2Pheno.
'''


gp_output_dir= 
# All the tables need to have the same type (binary/continuous)?
# Absolute paths only?
gp_input_tables_dir= 
gp_input_tables_type= 
gp_input_table_transpose=
gp_input_scf_dir=
gp_input_phenotype=
gp_input_tr=
gp_kmer_size=

gp_input_template= '''
<project output=\"{}\" name=\"{}\">
    <genotype>
        <tables path=\"{}\"
        normalization="{}" 
        transpose=\"{}\">

        {project_name}

        </tables>
        <sequence
        path=\"{}\"
        kmer=\'{}\'>
        </sequence>
    </genotype>

    <phenotype
    path=\"{}\">
    </phenotype>

    <phylogentictree path=\"{}\">
    </phylogentictree>
</project>
'''
print(gp_input_template)
