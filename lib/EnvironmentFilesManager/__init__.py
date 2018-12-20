class Pool:
    def __init__(self, env_files_home):
        self.env_files_home= env_files_home
    def find_yaml(self, env_name):
        import os
        # for activating the environment
        p= os.path.join(self.env_files_home, env_name+'.yaml')
        return(p)
    def update_variables_cmd(self, var_name, var_value, check_existence= False):
        import os
        # for updating the environment variables
        try: 
            if check_existence:
               self._test_env_var(var_name)
            cmd= 'export {}={}'.format(var_name, var_value)
            return(cmd)
        except Exception as e:
            import sys
            sys.exit(str(e))
    def _test_env_var(self, var_name):
        import os
        if os.environ.get(var_name) is None:
            raise Exception('ERROR: the environment variable "{}" was not found'.format(var_name))
    def call_pool_location(self):
        return(self.env_files_home)

