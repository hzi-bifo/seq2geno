class SoftwarePool:
    def __init__(self, env_sensitive= True):
        '''
        By turning the sensitivity on, the environmental paths 
        will be updaated every time when searching for sortwares
        '''
        import os
        self.env_sen= env_sensitive
        self.env_paths=os.environ.get('PATH') 
    def __str__(self):
        return(
        '{} environmental paths\ninitial environmental paths: {}'.format(
            'dynamic' if self.env_sen else 'fixed', self.env_paths))
    def find_software(self, software_name, target_dir= None, 
            include_env= True) -> str:
        '''
        priorty: lib -> PATH (by the environment) 
        '''
        import os

        libpaths= []
        # the target path, if specified
        if type(target_dir) is str:
            libpaths= [target_dir]
        # the environmental paths
        ## update the environement path?
        if include_env:
            env_path= self.env_paths if not self.env_sen else os.environ.get('PATH')
            env_paths= env_path.split(':')
            libpaths= libpaths+env_paths

        # filter by directory
        libpaths=filter(lambda x: os.path.isdir(x), libpaths) 
        software_candidates= [os.path.join(p, software_name) for p in libpaths]
        # filter by file
        software_candidates= [s for s in software_candidates if os.path.isfile(s)]

        if len(software_candidates) == 0:
            raise Exception('cannot find {}'.format(software_name))
        else:
            return(software_candidates[0])
