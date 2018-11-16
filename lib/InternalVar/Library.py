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
        "Rules@{}\nRecipes@{}\nLib@{}".format(self.rulespath, 
            self.recipespath, self.libpath))

    def call_rules_path(self):
        return(self.rulespath)

    def call_recipes_path(self):
        return(self.recipespath)

    def call_lib_path(self):
        return(self.libpath)
