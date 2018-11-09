class reader():
    def __init__(self, user_opt):
        self= self
        self.user_opt= user_opt

    def load_materials(self):
        from .Samples import Samples
        self.materials= Samples(self.user_opt)

    def load_reference(self):
        from .Reference import Reference
        self.reference= Reference(self.user_opt)

    def load_library(self):
        from .Library import Library
        self.library= Library(self.user_opt)


class reader_for_main(reader):
    def __init__(self, user_opt):
        self= self
        self.user_opt= user_opt

