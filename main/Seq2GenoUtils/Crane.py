import os
import time
from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.tools.inputs import inputs
import sys
import logging
import argparse

class Seq2Geno_Crane():
    def __init__(self,logger= '', email= '',
                 url= 'https://galaxy.bifo.helmholtz-hzi.de/galaxy/'):
        self.url=url
        self.email= ''
        self.logger= logger

    def _setup_key(self):
        '''
        This current one is registered under my email
        Remove before sending the script
        '''
        self.key= 'a118891c05aa8244b1f3b06cbfe8d60a'

    def _fetch_galaxy(self):
        '''upload the validated zip but not yet start the process'''
        assert hasattr(self, 'key'), 'Sumission key not set'
        self.gi = GalaxyInstance(url=self.url, key=self.key, verify= False)

    def _submit(self):
        self.history_id= self.gi.histories.get_histories()[0]['id']
        job_obj= ''
        job_obj=self.gi.tools.upload_file(self.seq2geno_input_zip, self.history_id)

        email_regex = '^(\w|\.|\_|\-)+[@](\w|\_|\-|\.)+[.]\w{2,3}$'
        job_id= job_obj['jobs'][0]['id']
        self.myinputs = inputs().set_dataset_param("input", job_id, src="hda")
        if self.email != '' and not (re.search(self.email, email_regex) is None):
            self.myinputs.set_param("email", self.email)
        else:
            print('Empty or invalid email address. Setup skipped')

    def launch(self, seq2geno_outdir, seq2geno_input_zip, 
             output_zip_f= 'seq2geno_results.zip',
             tool= 'seq2geno'):

        ####
        ##  SG zip submission to server
        self.seq2geno_input_zip= seq2geno_input_zip
        self._setup_key()
        self._fetch_galaxy() 
        if self.logger is logging.Logger:
            self.logger.info('Submitting {} to the server'.format(
                self.seq2geno_input_zip))
        self._submit() 
        # ensuring that the submission is done 
        history_obj= self.gi.histories.show_history(self.history_id, contents=False)
        if self.logger is logging.Logger:
            self.logger.info(
                'Submission state {}'.format(history_obj['state']))

        while history_obj['state'] != 'ok':
            assert history_obj['state'] != 'error', sys.exit('Data submission failed')
            time.sleep(100)
            history_obj= self.gi.histories.show_history(self.history_id, contents=False)
        if self.logger is logging.Logger:
            self.logger.info('Data submission done')

#        ####
#        ##  SG run with submitted data
#        self.gi.tools.run_tool(history_obj['id'],tool, self.myinputs)
#        # check status and download the result once finished
#        history_obj= self.gi.histories.show_history(self.history_id, contents=False)
#        while history_obj['state'] != 'ok':
#            # check it every 100s
#            time.sleep(100)
#            assert history_obj['state'] == 'error', (
#                'Error during modeling')
#            history_obj= self.gi.histories.show_history(self.history_id, contents=False)
#        if self.logger is logging.Logger:
#            self.logger.info('Genotyping finished')
#
#        ####
#        ## download the GP result
#        data_name= 'GenoPheno output'
#        datasets_obj= self.gi.datasets.get_datasets(history_id=self.history_id)
#        sg_out= [d for d in self.gi.datasets.get_datasets() 
#                 if d['name'] == data_name].pop()
#        self.gi.datasets.download_dataset(sg_out['dataset_id'], 
#                                     file_path=os.path.join(seq2geno_outdir,output_zip_f), 
#                                     use_default_filename=False)
#
#        self.logger.info('Data downloaded')
