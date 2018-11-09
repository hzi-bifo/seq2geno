class SoftwarePool:
    def find_software(self, software_name, target_dir= None):
        '''
        priorty: lib -> PATH (by the environment) -> the others
        '''
        import os

        libpaths= []
        if target_dir is None :
            libpaths= os.environ.get('PATH') .split(':')
        elif type(target_dir) is str:
            libpaths= [target_dir]
        elif type(target_dir) is list:
            libpaths= target_dir
        else:
            raise Exception(
            'illegal path to search for software')
        print('search \"{}\" in \"{}\"'.format(software_name, ':'.join(libpaths)))

        libpaths=filter(lambda x: os.path.isdir(x), libpaths) 
        software_candidates= [os.path.join(p, software_name) for p in libpaths if
                software_name in os.listdir(p)]

        if len(software_candidates) == 0:
            raise Exception('cannot find {}'.format(software_name))
        else:
            return(software_candidates[0])
