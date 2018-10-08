import re
import os

def divide_smk(smk_f, output_dir):
    rule_name= ''
    rule_dir= os.path.join(output_dir, rule_name)
    if not os.path.exists(rule_dir):
        os.makedirs(rule_dir)
    rule_f= os.path.join(rule_dir, '{}.smk'.format(rule_name))
    rule_fh= open(rule_f, 'w')
    for l in open(smk_f, 'r'):
        title_check= re.search('^rule\s*(\w+)\s*:', l)
        if l.startswith('#'):
            continue
        # detect rule name
        elif not (title_check is None):
            rule_fh.close()
            rule_name= title_check.group(1)
            rule_dir= os.path.join(output_dir, rule_name)
            if not os.path.exists(rule_dir):
                os.makedirs(rule_dir)
            rule_f= os.path.join(rule_dir, '{}.smk'.format(rule_name))
            rule_fh= open(rule_f, 'w')
        rule_fh.write(l)

    rule_fh.close()

# list the smk files
old_smk_dir= 'smk'
old_smk_files= [os.path.join(old_smk_dir, f) for f in os.listdir(old_smk_dir)
        if not( re.search('\.smk$', f) is None)]
# read the smk file
for smk_f in old_smk_files:
    divide_smk(smk_f, 'smk/rules')
