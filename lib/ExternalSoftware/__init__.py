class SoftwarePool:
    def __init__(self, env_dir):
        '''
        By turning the sensitivity on, the environmental paths 
        will be updaated every time when searching for sortwares
        '''
        import os
        self.env_dir= env_dir
    def __str__(self):
        return(
        'home of environments: {}'.format(
            self.env_dir))
    def find_software(self, software_name, target_dir= None, 
            target_env= None) -> str:
        '''
        priorty: target_dir -> target_env -> PATH (under seq2geno/ng_seq2geno) 
        '''
        import os

        libpaths= []
        # the target path, if specified
        if type(target_dir) is str:
            libpaths= [target_dir]
        # whether it's installed under the env
        elif type(target_env) is str:
            libpaths= [os.path.join(self.env_dir, target_env, 'bin')]
            
        # the others
        other_paths_str= os.environ.get('PATH')
        other_paths= other_paths_str.split(':')
        libpaths= libpaths+other_paths

        # filter by directory
        libpaths=filter(lambda x: os.path.isdir(x), libpaths) 
        software_candidates= [os.path.join(p, software_name) for p in libpaths]
        # filter by file
        software_candidates= [s for s in software_candidates if os.path.isfile(s)]

        if len(software_candidates) == 0:
            raise Exception('cannot find {}'.format(software_name))
        else:
            return(software_candidates[0])
