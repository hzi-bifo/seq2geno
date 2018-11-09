class reader():
    def __init__(self):
        self= self

    def load_materials(self, user_opt):
        from .Samples import Samples
        self.materials= Samples(user_opt)

class reader_for_main(reader):
    def __init__(self):
        self= self

