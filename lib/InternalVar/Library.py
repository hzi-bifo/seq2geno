class Library:
    def __init__(self, user_opt):
        import os
        self.rulespath= os.path.join(user_opt['seq2geno_smk_dir'], 
        'rules')
        self.recipespath= os.path.join(user_opt['seq2geno_smk_dir'], 
        'recipes')
        self.libpath= user_opt['seq2geno_lib_dir']
    
    def __str__(self):
        return(
        "Rules@{}\nRecipes@{}".format(self.rulespath, 
            self.recipespath))

    def find_rules(self):
        return(self.rulespath)

    def find_recipes(self):
        return(self.recipespath)
