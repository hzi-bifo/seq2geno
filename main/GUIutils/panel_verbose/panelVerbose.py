
class RedirectText(object):
    import tkinter as tk
    ''' 
    Class to extend the scrolledText
    '''
    def __init__(self, text_ctrl):
        self.output = text_ctrl
    def write(self, string):
        self.output.insert(tk.END, string)

def create_verbose_frame(win_root, outlevel, width= 1000, height= 100):
    '''
    make a panel to display the verbose output (stdout)
    '''
    from tkinter import scrolledtext
    from tkinter import ttk
    import sys
    panel_verbose= ttk.Frame(win_root, width= width, height= height)
    v_box= scrolledtext.ScrolledText(win_root, width= width)
    v_box.pack()
    redir= RedirectText(v_box)
    #' recognizable choices
    all_outlevels= ['stdout', 'stderr']
    outlevel= [l for l in all_outlevels if l in all_outlevels]
    if 'stdout' in outlevel:
        sys.stdout= redir
    if 'stderr' in outlevel:
        sys.stderr= redir
    return(v_box)
