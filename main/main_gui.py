#!/usr/bin/env python3
#' Role: Manager 
#' Purpose: 
#' Parse the arguments through the GUI

import tkinter as tk
from tkinter import scrolledtext
from tkinter import filedialog
from functools import partial
from tkinter import ttk
import sys
import os
import yaml
import re
from pprint import pprint

import UserOptions 
#' >>>>
#' Appearance attributes
#fieldname_font=("Helvetica",12,"bold") 
#def fieldname_cols(is_optional):
#    '''
#    field names colored by optional and mandatory fields
#    '''
#    return('red' if not is_optional else 'black')

#' >>>>
#' Target values to determine

config_dict= dict()
func_dict= dict()
primary_dict= dict()

#' >>>>
#' functions

def browseDirs(field): 
    '''
    allow the user to select a directory using the file browser and
    auto-fill the textbox
    '''
    d_name= filedialog.askdirectory(title= 'select directory', 
                                initialdir= '.')
    #' update the value
    config_dict[field].set(d_name)


def browse_file_and_update(field): 
    '''
    allow the user to select a file using the file browser and
    auto-fill the textbox
    '''
    f_name= filedialog.askopenfilename(title= 'select file', 
                                initialdir= '.',
                                filetypes=[('all', '*'),
                                        ('.tsv', '*.tsv'),
                                        ('.fa', '*.fa')])
    #' update the value
    config_dict[field].set(f_name)

def browse_primary_file_and_update(field): 
    '''
    allow the user to select a file using the file browser and
    auto-fill the textbox
    '''
    f_name= filedialog.asksaveasfilename(title= 'set filename', 
                                initialdir= '.')
    #' update the value
    primary_dict[field].set(f_name)

def make_primary_file_field(root, field, is_optional= False):
    #' the row
    row = ttk.Frame(root)
    row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=3)
    #' name of this field
    lab = ttk.Label(row, width=15, text=field)
    lab.pack(side=tk.LEFT)
    #' the text
    #' updated when filename selected or typed
    primary_dict[field]= tk.StringVar(row)
    ent = ttk.Entry(row, textvariable= primary_dict[field])
    ent.configure(font= ('nimbus mono l', 12))
    ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
    if is_optional:
        optout_but= ttk.Button(row, text= 'skip', width= 6, 
                               command=partial(field_opt_out, field))
        optout_but.pack(side=tk.RIGHT, padx=5, pady=5)
    but= ttk.Button(row, text= 'browse', width= 10, command= partial(browse_primary_file_and_update, field))
    but.pack(side=tk.RIGHT, padx=5, pady=5)

def makeform_primary(root, args_dict):
    '''
    create the form for determining input data
    '''
    for field in args_dict:
        if args_dict[field]['class'] == 'file':
            #' arguments of filenames
            make_primary_file_field(
                root, field, 
                (True if len(args_dict[field]['pattern'])==0 else False))

def field_opt_out(field):
    '''
    For the optional arguments, update the value with "-"
    '''
    config_dict[field].set('-')

def make_file_field_shared(root, field, is_optional= False):
    '''
    the field name and the textbox for all fields in common
    '''
    #' the row
    row = ttk.Frame(root)
    row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
    #' name of this field
    lab = ttk.Label(row, width=15, text=field)
    lab.pack(side=tk.LEFT)
    #' the text
    #' updated when filename selected or typed
    config_dict[field]= tk.StringVar(row)
    ent = ttk.Entry(row, textvariable= config_dict[field])
    ent.configure(font= ('nimbus mono l', 12))
    ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
    
    return(row)

def make_file_field(root, field, is_optional= False):
    '''
    the fields where the value should be a file
    '''
    row= make_file_field_shared(root, field, is_optional)
    if is_optional:
        optout_but= ttk.Button(row, text= 'skip', width= 6, 
                               command=partial(field_opt_out, field))
        optout_but.pack(side=tk.RIGHT, padx=5, pady=5)
    but= ttk.Button(row, text= 'browse', width= 10, command= partial(browse_file_and_update, field))
    but.pack(side=tk.RIGHT, padx=5, pady=5)

def make_dir_field(root, field):
    '''
    the fields where the value should be a directory
    '''
    row= make_file_field_shared(root, field)
    but= ttk.Button(row, text= 'browse', width= 10, command= partial(browseDirs, field))
    but.pack(side=tk.RIGHT, padx=5, pady=5)

def make_bool_field(root, field):
    row = ttk.Frame(root)
    row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
    lab = ttk.Label(row, width=15, text=field)
    lab.pack(side=tk.LEFT)
    config_dict[field]= tk.StringVar(row)
    opt_but= ttk.OptionMenu(row, config_dict[field],
                                'N', *['N', 'Y'])
    opt_but.pack(side=tk.LEFT)


def makeform_general(root, args_dict):
    '''
    create the form for determining input data
    '''
    for field in args_dict:
        if args_dict[field]['class'] == 'file':
            #' arguments of filenames
            make_file_field(
                root, field, 
                (True if len(args_dict[field]['pattern'])==0 else False))
        elif args_dict[field]['class'] == 'dir':
            #' arguments of directories
            make_dir_field(root, field)
        elif args_dict[field]['class'] == 'bool':
            #' arguments of directories
            make_bool_field(root, field)
        else:
            #' other types of arguments
            make_file_field_shared(root, field)

def makeform_functions(root, func_options):
    '''
    options of functions
    '''
    for func in func_options:
        row = ttk.Frame(root)
        row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        lab = ttk.Label(row, width=10, text=func, anchor='w')
        lab.pack(side=tk.LEFT)
        func_dict[func]= tk.StringVar(row)
        #' the options of ttk.OptionMenu differ from those of tk.OptionMenu 
        opt_but= ttk.OptionMenu(row, func_dict[func],
                                func_options[func][0], *func_options[func])
        opt_but.pack(side=tk.LEFT)


#def func_name_adaptor(d):
#    '''
#    Becasue the function names used in the main script is abbreviated,
#    this interface describe thw workflows in a more understoodable way.
#    To allow the main script recognize them,
#    this adaptor convert the names and values
#    '''
#    #' ensure the datatype first
#    assert isinstance(d, dict)
#    #' the dictionary for the conversion
#    func_name_adaptor={ 'mode': 'dryrun',
#            'snps': 'snps',
#            'denovo': 'denovo',
#            'expr':'expr',
#            'phylo':'phylo',
#            'ar':'ar',
#            'de': 'de'  }
#    func_val_adaptor= { 'include': 'Y', 'skip': 'N',
#                       'dryrun': 'Y','execute': 'N'}
#    new_d= dict()
#    for k in func_name_adaptor:
#        old_val= d[k]
#        new_val= func_val_adaptor[old_val]
#        new_key= func_name_adaptor[k]
#        new_d[new_key]= new_val
#    return(new_d)


def load_theme(root):
    parent_d=os.path.dirname(__file__)
    awtheme_d=os.path.join(parent_d, 'GUIutils', 'theme', 'awthemes-10.0.0')
    root.tk.call('lappend', 'auto_path', awtheme_d)
    root.tk.call('package', 'require', 'awdark')
    s= ttk.Style()
    s.theme_use('awdark')
    s.configure('.', font=('Helvetica', 12))


def make_arguments_for_main(func_dict, config_dict, argspace):
    ''' 
    prepare the argument object for Seq2Geno
    '''
    #' encode the boolean variables
    #' features section
    func_plainstr_dict= {k: func_dict[k].get() for k in func_dict}
    func_dict_for_main= func_plainstr_dict
    pprint(func_dict_for_main)
    #' general section
    config_plainstr_dict= {k: config_dict[k].get() for k in config_dict}
    config_dict_for_main= config_plainstr_dict
    pprint(config_plainstr_dict)

    #' create the arguments object
    args= UserOptions.arguments()
    try:
        args.add_opt(**func_dict_for_main)
        args.add_opt(**config_plainstr_dict)
    except KeyError as e:
        sys.exit('ERROR: {} not found in the input file'.format(str(e)))
    else:
        print(args.__dict__)
        args.check_args()
        return(args)

def read_arguments_space(as_f):
    as_fh= open(as_f, 'r')
    args= yaml.safe_load(as_fh)
    as_fh.close()
    return(args)

def load_old_yaml():
    '''
    Allow the user to select the old yaml file and parse the arguments
    '''
    yml_f= filedialog.askopenfilename(title= 'select file', 
                                initialdir= '.',
                                filetypes=[('yml', '*.yml'),
                                        ('.yaml', '*.yaml'),
                                        ('all', '*')])
    reset_args_with_yml(yml_f)

def reset_args_with_yml(yml_f):
    if len(yml_f) > 0 :
        primary_dict['yml_f'].set(yml_f)
        #' in case the selection is canceled or accidents
        old_args= UserOptions.parse_arg_yaml(yml_f)
        for func in func_dict:
            if hasattr(old_args, func):
                arg_val= ('Y' if getattr(old_args, func)=='Y' else 'N')
                #print('{}: {} '.format(func, func_dict[func].get()))
                func_dict[func].set(arg_val)
                #print('---> {}'.format(func_dict[func].get()))
        for k in config_dict:
            if hasattr(old_args, k):
                arg_val= getattr(old_args, k)
                #print('{}: {} '.format(k, config_dict[k].get()))
                config_dict[k].set(arg_val)
                #print('---> {}'.format(config_dict[k].get()))

def write_yaml(func_dict, config_dict):
    '''
    Save the setting 
    '''
    import yaml
    func_plainstr_dict= {k: func_dict[k].get() for k in func_dict}
    config_plainstr_dict= {k: config_dict[k].get() for k in config_dict}
    arg_dict={'features': func_plainstr_dict, 'general': config_plainstr_dict}
    #' determine the filename
    yml_f= filedialog.asksaveasfilename(title= 'save as...', 
                                initialdir= '.')
    if len(yml_f) > 0 :
        primary_dict['yml_f'].set(yml_f)
        with open(yml_f, 'w') as yml_fh:
            yaml.safe_dump(arg_dict, yml_fh)
        return(yml_f)

def load_old_log(out):
    '''
    Read the old log file to allow the user to know the current status
    '''
    #' open the log file
    import UserOptions 
    import os
    log_f= filedialog.askopenfilename(title= 'select file', 
                                initialdir= '.',
                                filetypes=[('log', '*.log'),
                                        ('all', '*')])
    primary_dict['log_f'].set(log_f)

    #' print the information in the log file
    yml_f= ''
    out.delete('1.0', tk.END)
    with open(log_f, 'r') as log_fh:
        for l in log_fh.readlines():
            is_config_line= re.search('#CONFIGFILE:(.+)', l.strip())
            if not is_config_line is None:
                yml_f= is_config_line.group(1)
            else:
                out.insert(tk.END, l)
    
    #' set the arguments same as the config file that generated the log
    assert os.path.isfile(yml_f), (
        'Config file "{}" described in the log not found or broken'.format(yml_f))
    reset_args_with_yml(yml_f)


class seq2geno_gui:
    def __init__(self, root):
        self.win_root= root
    def show(self):
        win_root= self.win_root
        parent_d=os.path.dirname(__file__)

        #' primary arguments that the user would set with the commandline interface
        #' read the arguments space
        p_as_f=os.path.join(parent_d, 'GUIutils', 'PrimaryArgSpace.yml')
        self.p_argspace= read_arguments_space(p_as_f)
        panel_primary=ttk.Frame(win_root, width= 1000, borderwidth= 10)
        makeform_primary(panel_primary, self.p_argspace)


        #' arguments listed in the input yaml file
        #' read the arguments space
        as_f=os.path.join(parent_d, 'GUIutils', 'ArgSpace.yml')
        self.argspace= read_arguments_space(as_f)
        win_mainframe= ttk.Notebook(win_root, width= 1000)
        #' group the arguments
        #' options of workflows 
        panel_functions= ttk.Frame(win_mainframe)
        win_mainframe.add(panel_functions, text= 'features')
        makeform_functions(panel_functions, self.argspace['features'])

        #' IO panel 
        panel_general= ttk.Frame(win_mainframe)
        win_mainframe.add(panel_general, text= 'general')
        makeform_general(panel_general, self.argspace['general'])

        #' the log display 
        panel_log= ttk.Frame(win_mainframe)
        win_mainframe.add(panel_log, text= 'log')
        v_box= scrolledtext.ScrolledText(panel_log, height=200)
        v_box.configure(font=('nimbus mono l', 12))
        v_box.pack()

        #' the menu
        #' file
        win_menubar= tk.Menu(win_root)
        filemenu= tk.Menu(win_menubar, tearoff= False)
        filemenu.add_command(label= 'Load yaml', command=load_old_yaml)
        filemenu.add_command(
            label= 'Save yaml',
            command=partial(write_yaml,func_dict, config_dict))
        filemenu.add_command(label= 'Load log', 
                             command=partial(load_old_log, v_box))
        filemenu.add_command(label= 'Run', command=self.exec)
        filemenu.add_command(label= 'Exit', command=win_root.quit)
        win_menubar.add_cascade(menu= filemenu, label= 'File')

#        panel_verbose= ttk.Frame(win_root, height= 200)
#        v_box= scrolledtext.ScrolledText(panel_verbose, height=200)
#        v_box.configure(font=('nimbus mono l', 12))
#        v_box.pack()
#        redir= RedirectText(v_box)
#        #' recognizable choices
#        sys.stdout= redir
#        sys.stderr= redir
#
        #' theming
        panel_primary.pack(fill= tk.X)
        win_mainframe.pack(fill= tk.X)
        win_root.config(menu= win_menubar)
#        panel_verbose.pack(fill=tk.X)
        load_theme(win_root)
        win_root.mainloop()
    def extract_args(self):
        #' collect the arguments
        args_for_main= make_arguments_for_main(func_dict,
                                               config_dict,
                                               self.argspace)
        self.args= args_for_main
        return(self.args)
    def exec(self):
        #' Instead of replacing the command line interface,
        #' the GUI is a an adaptor between the user and the command line interface.
        #' Therefore, launching seq2geno from this GUI simply is passing the needed 
        #' primary arguments: the yaml filename and the log filename
        yml_f= write_yaml(func_dict, config_dict)
        #log_f= 

#        import seq2geno
#        self.extract_args()
#        seq2geno.main(self.args)


if __name__ == '__main__':
    win_root = tk.Tk()
    win_root.configure(background='grey')
    win_root.title('Seq2Geno')
    seq2geno_gui= seq2geno_gui(win_root)
    seq2geno_gui.show()
    print(seq2geno_gui.extract_args().__dict__)
#    seq2geno_gui.exec()
