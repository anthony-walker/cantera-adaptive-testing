class Surface(object):
    def __init__(self):
        self.surface = None
        self.sphase = None
        self.fuel = ""


class PlatinumSmall(Surface):
    def __init__(self, *args, **kwargs):
        super(PlatinumSmall, self).__init__(*args, **kwargs)
        self.surface = 'PT(S):1.0'
        self.fuel = "CH4:1.0"
        self.sphase = 'surface-small'


class PlatinumMedium(Surface):
    def __init__(self, *args, **kwargs):
        super(PlatinumMedium, self).__init__(*args, **kwargs)
        self.surface = 'Pt(9):1.0'
        self.fuel = "CH4:1.0"
        self.sphase = 'surface-medium'


class PlatinumLarge(Surface):
    def __init__(self, *args, **kwargs):
        super(PlatinumLarge, self).__init__(*args, **kwargs)
        self.surface = 'Pt(9):0.1'
        self.fuel = "CH4:1.0"
        self.sphase = 'surface-large'
