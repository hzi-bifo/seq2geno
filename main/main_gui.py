#!/usr/bin/env python3
#' Role: Manager 
#' Purpose: 
#' Parse the arguments through the GUI

import tkinter as tk
from functools import partial
from tkinter import ttk
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
config_plainstr_dict= dict()
func_dict= dict()
func_plainstr_dict= dict()

#' >>>>
#' functions

def browseDirs(field): 
    '''
    allow the user to select a directory using the file browser and
    auto-fill the textbox
    '''
    from tkinter import filedialog
    d_name= filedialog.askdirectory(title= 'select directory', 
                                initialdir= '.')
    #' update the value
    config_dict[field].set(d_name)

def browseFiles(field): 
    '''
    allow the user to select a file using the file browser and
    auto-fill the textbox
    '''
    from tkinter import filedialog
    f_name= filedialog.askopenfilename(title= 'select file', 
                                initialdir= '.',
                                filetypes=[('all', '*'),
                                        ('.tsv', '*.tsv'),
                                        ('.fa', '*.fa')])
    #' update the value
    config_dict[field].set(f_name)

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
    but= ttk.Button(row, text= 'browse', width= 10, command= partial(browseFiles, field))
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
    import os
    parent_d=os.path.dirname(__file__)
    awtheme_d=os.path.join(parent_d, 'GUIutils', 'theme', 'awthemes-10.0.0')
    root.tk.call('lappend', 'auto_path', awtheme_d)
    root.tk.call('package', 'require', 'awdark')
    s= ttk.Style()
    s.theme_use('awdark')

def make_arguments_for_main(func_dict, config_dict, argspace):
    ''' 
    prepare the argument object for Seq2Geno
    '''
    from pprint import pprint
    #' encode the boolean variables
    #' features section
    func_plainstr_dict= {k: func_dict[k].get() for k in func_dict}
    func_dict_for_main= func_plainstr_dict
    for k in func_dict_for_main:
        func_dict_for_main[k]= (True if func_plainstr_dict[k] == 'Y'
                                      else False)
    pprint(func_dict_for_main)
    #' general section
    config_plainstr_dict= {k: config_dict[k].get() for k in config_dict}
    config_dict_for_main= config_plainstr_dict
    for k in config_dict_for_main:
        if argspace['general'][k]['class'] == 'bool' :
            config_dict_for_main[k]= (True if config_plainstr_dict[k] == 'Y'
                                      else False)
    pprint(config_plainstr_dict)
    import UserOptions 
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

def read_arguments_space():
    import yaml
    import os
    parent_d=os.path.dirname(__file__)
    as_f=os.path.join(parent_d, 'GUIutils', 'ArgSpace.yml')
    as_fh= open(as_f, 'r')
    args= yaml.safe_load(as_fh)
    as_fh.close()
    return(args)

def load_old_yaml():
    '''
    Allow the user to select the old yaml file and parse the arguments
    '''
    from tkinter import filedialog
    import UserOptions 
    import os
    yml_f= filedialog.askopenfilename(title= 'select file', 
                                initialdir= '.',
                                filetypes=[('yml', '*.yml'),
                                        ('.yaml', '*.yaml'),
                                        ('all', '*')])
    if len(yml_f) > 0 :
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

class seq2geno_gui:
    def __init__(self, root):
        self.win_root= root
    def show(self):
        win_root= self.win_root
        #' read the arguments space
        self.argspace= read_arguments_space()

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

        #' the menu
        #' file
        win_menubar= tk.Menu(win_root)
        filemenu= tk.Menu(win_menubar, tearoff= False)
        filemenu.add_command(label= 'Load yaml', command=load_old_yaml)
        filemenu.add_command(label= 'Exit', command=win_root.quit)
        win_menubar.add_cascade(menu= filemenu, label= 'File')
        #' theming
        win_mainframe.pack()
        win_root.config(menu= win_menubar)
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
        import seq2geno
        seq2geno.main(self.args)


if __name__ == '__main__':
    win_root = tk.Tk()
    win_root.configure(background='grey')
    win_root.title('Seq2Geno')
    seq2geno_gui= seq2geno_gui(win_root)
    seq2geno_gui.show()
    print(seq2geno_gui.extract_args().__dict__)
    seq2geno_gui.exec()
