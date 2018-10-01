#! /usr/bin/env python
import os

if __name__== '__main__':
    # create the config file
    import UserOptions
    args, config_f= UserOptions.main()    
    print(args)
    print(config_f)
