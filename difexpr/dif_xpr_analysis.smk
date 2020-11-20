#' Purpose:
#' - Determine differential expressions
#' Materials:
#' - expression levels  
#' - phenotypes
#' Methods:
#' - differnetial expression analyzed by DESeq2
#' Output:
#' - list of differntially expressed genes
#
import yaml
DIF_XPR_OUT=config['out_dir']
PHE_TABLE_F=config['pheno_tab']
expr_config_f= config['expr_config']
expr_config= yaml.load(open(expr_config_f, 'r'))
EXPR_OUT= os.path.join(os.path.dirname(expr_config_f), expr_config['out_table'])
DIFXPR_ALPHA=config['alpha']
DIFXPR_LFC=config['lfc']

rule all:
    input:
        dif_xpr_out=directory(DIF_XPR_OUT)
        
rule dif_xpr_analysis:
    input:
        phe_table= PHE_TABLE_F,
        expr_table= EXPR_OUT
    output:
        dif_xpr_out=directory(DIF_XPR_OUT)
    params:
        alpha_cutoff=DIFXPR_ALPHA,
        lfc_cutoff= DIFXPR_LFC
    threads: 12 
    conda: 'dif_expr_env.yml'
    script: 'diffEpr_analysis.R'
        
