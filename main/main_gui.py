#!/usr/bin/env python3
#' Role: Manager 
#' Purpose: 
#' Determine, initiate and launch the workflows based on user-defined
#' arguments

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
    d_name= filedialog.askdirectory(title= 'select firectory', 
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
        optout_but= ttk.Button(row, text= 'skip', width= 10, command=
                          partial(field_opt_out, field))
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

def make_plain_field(root, field):
    row= make_file_field_shared(root, field)

def makeform_general(root):
    '''
    create the form for determining input data
    '''
    is_optional_dict= {'old_config': True, 'dna_reads': False, 
                       'ref_fa': False, 'ref_gbk': False,'ref_gff': False,  
                       'rna_reads': True,
                       'adaptor': True}
    for field in is_optional_dict:
        make_file_field(root, field, is_optional_dict[field])
    dir_fields= ['wd']
    for field in dir_fields:
        make_dir_field(root, field)
    #' fields that do not need the file browser
    make_file_field_shared(root, 'cores')

def makeform_functions(root):
    '''
    options of functions
    '''
    func_options= {'snps': ['skip', 'include'],
        'denovo': ['skip','include'],
        'expr': ['skip','include'],
        'phylo': ['skip','include'],
        'ar': ['skip','include'],
        'de': ['skip','include'],
        'mode': ['dryrun', 'execute']}
    for func in func_options:
        row = ttk.Frame(root)
        row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        lab = ttk.Label(row, width=15, text=func, anchor='w')
        lab.pack(side=tk.LEFT)
        func_dict[func]= tk.StringVar(row)
        func_dict[func].set(func_options[func][0])
        opt_but= ttk.OptionMenu(row, func_dict[func], *func_options[func])
        opt_but.pack(side=tk.LEFT)

def func_name_adaptor(d):
    '''
    Becasue the function names used in the main script is abbreviated,
    this interface describe thw workflows in a more understoodable way.
    To allow the main script recognize them,
    this adaptor convert the names and values
    '''
    #' ensure the datatype first
    assert isinstance(d, dict)
    #' the dictionary for the conversion
    func_name_adaptor={ 'mode': 'dryrun',
            'snps': 'snps',
            'denovo': 'denovo',
            'expr':'expr',
            'phylo':'phylo',
            'ar':'ar',
            'de': 'de'  }
    func_val_adaptor= { 'include': 'Y', 'skip': 'N',
                       'dryrun': 'Y','execute': 'N'}
    new_d= dict()
    for k in func_name_adaptor:
        old_val= d[k]
        new_val= func_val_adaptor[old_val]
        new_key= func_name_adaptor[k]
        new_d[new_key]= new_val
    return(new_d)


#class ToggledFrame(ttk.Frame):
#    def __init__(self, root, name):
#        self.show = tk.IntVar()
#        self.show.set(0)
#        from tkinter import ttk
#        self.main_frame= ttk.Frame(root)
#        lab = ttk.Label(self.main_frame, width=15, text=name, anchor='w')
#        lab.pack(side=tk.LEFT)
#        toggle_button = ttk.Checkbutton(self.main_frame, width=2, text='+', command=self.toggle,
#                                            variable=self.show, style='Toolbutton')
#        toggle_button.pack(side="left")
#        self.sub_frame = ttk.Frame(self, relief="sunken", borderwidth=1)
#    def toggle(self):
#        if bool(self.show.get()):
#            self.sub_frame.pack(fill="x", expand=1)
#            self.toggle_button.configure(text='-')
#        else:
#            self.sub_frame.forget()
#            self.toggle_button.configure(text='+')

def load_theme(root):
    import os
    parent_d=os.path.dirname(__file__)
    awtheme_d=os.path.join(parent_d, 'theme', 'awthemes-10.0.0')
    root.tk.call('lappend', 'auto_path', awtheme_d)
    root.tk.call('package', 'require', 'awdark')

def make_arguments_for_main(func_dict, config_dict):
    ''' 
    prepare the argument object for Seq2Geno
    '''
    from pprint import pprint
    config_plainstr_dict= {k: config_dict[k].get() for k in config_dict}
    pprint(config_plainstr_dict)
    func_plainstr_dict= {k: func_dict[k].get() for k in func_dict}
    func_plainstr_dict_for_main= func_name_adaptor(func_plainstr_dict)
    pprint(func_plainstr_dict_for_main)
    import UserOptions 
    args= UserOptions.arguments()
    try:
        args.add_opt(**func_plainstr_dict_for_main)
        args.add_opt(**config_plainstr_dict)
    except KeyError as e:
        sys.exit('ERROR: {} not found in the input file'.format(str(e)))
    else:
        print(args.__dict__)
        args.check_args()
        return(args)

if __name__ == '__main__':
    win_root = tk.Tk()
    win_root.configure(background='black')
    win_root.title('Seq2Geno')
    win_root.geometry('1200x1000')

    #' functions
    panel_functions= ttk.Frame(win_root).pack(side=tk.TOP)
    makeform_functions(panel_functions)
    #' the configurable arguments
    panel_general= ttk.Frame(win_root).pack(side= tk.TOP)
    makeform_general(panel_general)
    save_but = ttk.Button(win_root, text='Save', command=win_root.quit)
    save_but.pack(side=tk.TOP, padx=5, pady=5)
    #' theming
    load_theme(win_root)
#    win_root.tk.call('lappend', 'auto_path',
#                     '/home/thkuo/projects/test_tkinter/awthemes-10.0.0')
#    win_root.tk.call('package', 'require', 'awdark')
    s= ttk.Style()
    s.theme_use('awdark')
    win_root.mainloop()
    #' collect the arguments
    args_for_main= make_arguments_for_main(func_dict, config_dict)
    print(args_for_main)
