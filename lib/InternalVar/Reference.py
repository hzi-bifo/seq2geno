class Reference:
    def __init__(self, user_opt):
        self.REF_FA=user_opt['ref_fa']
        self.REF_GBK=user_opt['ref_gbk']
    
    def __str__(self):
        return(
        "Reference sequence: {}\nReference annotation: {}".format(self.REF_FA, 
            self.REF_GBK))

    def call_ref_seq_file(self):
        return(self.REF_FA)

    def call_ref_anno_file(self):
        return(self.REF_GBK)
