class EnvironmentManager:
    def __init__(self, env_home):
        self.env_home= env_home
    def activate_env_cmd(self, env_name):
        import os
        # for activating the environment
        cmd= 'source activate {}'.format(os.path.join(self.env_home, env_name))
        return(cmd)
    def update_variables_cmd(self, var_name, var_value):
        import os
        # for updating the environment variables
        if environ.get(var_name) is not None:
            cmd= 'export {}={}'.format(var_name, var_value)
            return(cmd)
        else:
            raise Exception('Variable "{}" not found'.format(var_name))
            
