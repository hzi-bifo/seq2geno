def random_filename(n=12):
    '''
    For generating an unique filename that fills in the rules, although it won't
be really created. 
    '''
    f= ''.join(random.choice(string.ascii_uppercase + string.digits) for x in
range(n))
    return(f)
