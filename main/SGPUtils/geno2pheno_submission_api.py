'''
To resolve:
    Where the Clinet is?
    key acuired by function
'''

import os
import time
from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.tools.inputs import inputs
import sys
import logging
# the validator
from .Geno2PhenoClient import validator
import .create_genyml

def make_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(filename)s %(asctime)s %(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return(logger)

def setup_key():
    '''
    This current one is registered under my email
    Remove before sending the script
    '''
    key= 'a118891c05aa8244b1f3b06cbfe8d60a'
    return(key)

def validate_and_compress(config, log, gp_zip):
    assert os.path.isfile(config), logger.error(
        '{} not found'.format(config))
    try:
    # run the validator
        val= validator.ValidateGenML(config, log)
        val.create_zip_file(gp_zip)
        print('Compression done')
    except IOError as e:
        logger.error('Validation or compression failed. Exit')
        sys.exit()


def submit_gp_zip( key, geno2pheno_input_zip, email= ''):
    '''upload the validated zip but not yet start the process'''
    url= 'https://galaxy.bifo.helmholtz-hzi.de/galaxy/'
    gi = GalaxyInstance(url=url, key=key, verify= False)
    history_id= gi.histories.get_histories()[0]['id']
    job_obj= ''
    job_obj=gi.tools.upload_file(geno2pheno_input_zip, history_id)

    email_regex = '^(\w|\.|\_|\-)+[@](\w|\_|\-|\.)+[.]\w{2,3}$'
    job_id= job_obj['jobs'][0]['id']
    myinputs= ''
    if email != '' and re.search(email, email_regex):
        myinputs = inputs().set_dataset_param("input", job_id, src="hda").set_param(
            "email", email)
    elif email != '' and re.search(email, email_regex) is None: 
        raise TypeError('Incorrect email format')
    else:
        myinputs = inputs().set_dataset_param("input", job_id, src="hda")
    history_obj= gi.histories.show_history(history_id, contents=False)
    return(gi, history_obj, myinputs)

project_name= 'sgp'
gp_config= './kp_config.yml'
seq2geno_outdir=('/net/metagenomics/data/from_moni/old.tzuhao/'
                 'seq2geno/test_Kpneumoniae/v4/seq2geno_test/')
log= 'geno2phenoclient.log'
geno2pheno_input_zip= './kp.zip'
geno2pheno_out= 'test_geno2phenoclient_out.zip' 
email= ''

####
##  SG output -> genml
genml_creator_args_array= ['--seq2geno', seq2geno_outdir,
                         '--yaml', gp_config, 
                         '--proj', project_name]
genml_creator_args_parser= create_genyml.make_parser()
genml_creator_args= genml_creator_args_parser.parse_args(genml_creator_args_array)
try:
    create_genyml.make_genyml(genml_creator_args) 
except IOError as e:
    print('Creation of {} failed'.format(gp_config))

####
##  genml -> GP zip
key= setup_key()
# create the zip file
validate_and_compress(gp_config, log, geno2pheno_input_zip)
try :
    assert os.path.isfile(geno2pheno_input_zip)
except AssertionError as e:
    print('{} not found'.format(geno2pheno_input_zip))

if geno2pheno_out != '':
    assert not os.path.isfile(geno2pheno_out), 'Output path exists!'
else:
    os.makedirs(os.path.dirname(geno2pheno_out))

####
##  GP zip submission to server
gi, history_obj, myinputs= submit_gp_zip(key=key,
              geno2pheno_input_zip= geno2pheno_input_zip,
              email=email) 
# ensuring that the submission is done 
while history_obj['state'] != 'ok':
    assert history_obj['state'] == 'error', sys.exit('Data submission failed')
    time.sleep(100)
print('Data submission done')

####
##  GP run with submitted data
gi.tools.run_tool(history_obj['id'], 'genopheno', myinputs)
# check status and download the result once finished
history_obj= gi.histories.show_history(history_id, contents=False)
while history_obj['state'] != 'ok':
    # check it every 100s
    time.sleep(100)
    assert history_obj['state'] == 'error', (
        'Error during modeling')
    history_obj= gi.histories.show_history(history_id, contents=False)
print('Predictive modeling finished')

####
## download the GP result
datasets_obj= gi.datasets.get_datasets(history_id=history_id)
gp_out= [d for d in gi.datasets.get_datasets() 
         if d['name'] == 'GenoPheno output'].pop()
gi.datasets.download_dataset(gp_out['dataset_id'], 
                             file_path=geno2pheno_out, 
                             use_default_filename=False)
