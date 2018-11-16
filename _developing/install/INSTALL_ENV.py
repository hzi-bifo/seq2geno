import os
import subprocess
import sys
home=os.path.dirname(os.path.realpath(__file__))
# the environments
env_yml_dir= os.path.join(home, 'env', 'yaml')
env_home= os.path.join(home, 'env')
env_list_f= 'ENV_LIST'
env_names= [l.strip() for l in open(env_list_f, 'r')]
print('environments to install:\n\t{}'.format(
    ' '.join(env_names)))

'''
for env_name in env_names[1:]:
    print('installing {}...'.format(env_name))
    env_dir= os.path.join(env_home, env_name)
    if not os.path.exists(env_dir):
        os.makedirs(env_dir)

    print(env_dir)
    print(os.path.exists(env_dir))
    conda_cmd= ['conda', 'env', 'create', '-f', 
        os.path.join(env_yml_dir, env_name+'.yaml'),
        '--prefix', 
        os.path.join(env_home, env_name)]
    print(' '.join(conda_cmd))
    try:
        subprocess.run(conda_cmd)
    except AttributeError:
        subprocess.call(conda_cmd)
    except Exception as e:
        sys.exit('Errors in the installation of env "{}":\n\t{}'.format(env_name, str(e)))
'''
