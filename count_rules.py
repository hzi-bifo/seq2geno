import re
import os

def divide_smk(smk_f):
    titles= []
    for l in open(smk_f, 'r'):
        title_check= re.search('^rule\s*(\w+)\s*:', l)
        if l.startswith('#'):
            continue
        # detect rule name
        elif not (title_check is None):
            rule_name= title_check.group(1)
            if rule_name == 'all':
                continue
            else:
                titles.append(rule_name)

    return(titles)
# list the smk files
old_smk_dir= 'smk/old'
old_smk_files= [os.path.join(old_smk_dir, f) for f in os.listdir(old_smk_dir)
        if not( re.search('\.smk$', f) is None)]
# read the smk file
rule_counts= {}
for smk_f in old_smk_files:
    rules= divide_smk(smk_f)
    for r in rules:
        rule_counts[r]= (rule_counts[r]+1 if r in rule_counts else 1)
rules= sorted(rule_counts.keys(), key= lambda x: rule_counts[x], reverse= True)

for r in rules:
    print('{}\t{}'.format(rule_counts[r], r))

