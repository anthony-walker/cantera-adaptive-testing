from .model_base import ModelBase
import cantera as ct


class GAS_DEF(ModelBase):
    def __init__(self, *args, **kwargs):
        super(GAS_DEF, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('gas-def.yaml')
        self.fuel = 'CH4:1.0'


class Hydrogen(ModelBase):
    def __init__(self, *args, **kwargs):
        super(Hydrogen, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('hydrogen-10-28.yaml')
        self.fuel = 'H2:1.0'


class GRI(ModelBase):
    def __init__(self, *args, **kwargs):
        super(GRI, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('gri-mech-55-325.yaml')
        self.fuel = 'CH4:1.0'

class A1(ModelBase):
    def __init__(self, *args, **kwargs):
        super(A1, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('a1-119-841.yaml')
        self.fuel = 'POSF10264:1.0'

class A2(ModelBase):
    def __init__(self, *args, **kwargs):
        super(A2, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('a2-119-841.yaml')
        self.fuel = 'POSF10325:1.0'

class A3(ModelBase):
    def __init__(self, *args, **kwargs):
        super(A3, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('a3-119-841.yaml')
        self.fuel = 'POSF10289:1.0'

class C1(ModelBase):
    def __init__(self, *args, **kwargs):
        super(C1, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('c1-119-841.yaml')
        self.fuel = 'POSF11498:1.0'

class C5(ModelBase):
    def __init__(self, *args, **kwargs):
        super(C5, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('c5-119-843.yaml')
        self.fuel = 'POSF12345:1.0'

class JP10(ModelBase):
    def __init__(self, *args, **kwargs):
        super(JP10, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('jp10-120-841.yaml')
        self.fuel = 'C10H16:1.0'
class R1(ModelBase):
    def __init__(self, *args, **kwargs):
        super(R1, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('r1-119-841.yaml')
        self.fuel = 'POSF7688:1.0'

class R2(ModelBase):
    def __init__(self, *args, **kwargs):
        super(R2, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('r2-119-841.yaml')
        self.fuel = 'POSF5433:1.0'

class DME(ModelBase):
    def __init__(self, *args, **kwargs):
        super(DME, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('dme-propane-122-711.yaml')
        self.fuel = 'CH3OCH3:1.0'

class JetA(ModelBase):
    def __init__(self, *args, **kwargs):
        super(JetA, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('jetA-detailed-NOx-203-1589.yaml')
        self.fuel = 'POSF10325:1.0'


class Butane(ModelBase):
    def __init__(self, *args, **kwargs):
        super(Butane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('butane-230-2461.yaml')
        self.fuel = 'C4H10:1.0'


class TwoButonane(ModelBase):
    def __init__(self, *args, **kwargs):
        super(TwoButonane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path(
            '2-butonane-ch3coch2ch3-315-1803.yaml')
        self.fuel = 'C4H8O1-2:1.0'


class NHexanal(ModelBase):
    def __init__(self, *args, **kwargs):
        super(NHexanal, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('n-hexanal-482-5182.yaml')
        self.fuel = 'C6H14:1.0'


class Aramco(ModelBase):
    def __init__(self, *args, **kwargs):
        super(PlatinumAramco, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('aramco-493-2716.yaml')
        self.fuel = 'CH4:1.0, C3H8:1.0, C2H6:1.0'


class IsoButene(ModelBase):
    def __init__(self, *args, **kwargs):
        super(IsoButene, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('isobutene-ic4h8-493-2716.yaml')
        self.fuel = 'IC4H8:1.0'


class IsoPentanol(ModelBase):
    def __init__(self, *args, **kwargs):
        super(IsoPentanol, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('i-pentanol-511-3010.yaml')
        self.fuel = 'ic5h10oh:1.0'


class T124MCH(ModelBase):
    def __init__(self, *args, **kwargs):
        super(T124MCH, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('t124mch-533-3193.yaml')
        self.fuel = 'T124MCH:1.0'


class OneTwoDME(ModelBase):
    def __init__(self, *args, **kwargs):
        super(OneTwoDME, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('one-two-dme-570-2960.yaml')
        self.fuel = 'CH3OCH3:1.0'


class PropylAcetate(ModelBase):
    def __init__(self, *args, **kwargs):
        super(PropylAcetate, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('propyl-acetate-628-4182.yaml')
        self.fuel = 'pa:1.0'


class NHeptane(ModelBase):
    def __init__(self, *args, **kwargs):
        super(NHeptane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('n-heptane-c7h16-654-4846.yaml')
        self.fuel = 'NC7H16:1.0'


class DEE(ModelBase):
    def __init__(self, *args, **kwargs):
        super(DEE, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('dee-746-3555.yaml')
        self.fuel = 'DEE:1.0'


class IsoOctane(ModelBase):  # Iso-Octane
    def __init__(self, *args, **kwargs):
        super(IsoOctane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('ic8-874-6864.yaml')
        self.fuel = 'IC8H18:1.0'


class GasolineSurrogate(ModelBase):
    def __init__(self, *args, **kwargs):
        super(GasolineSurrogate, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path("gasoline-surrogate-1389.yaml")
        self.fuel = "C6H5CH3:1.0, NC7H16:1.0, IC8H18:1.0"


class ThreeMethylHeptane(ModelBase):
    """“Kinetic modeling of gasoline surrogate components and mixtures
under engine conditions | Elsevier Enhanced Reader.”
https://reader.elsevier.com/reader/sd/pii/S1540748910000787?token=7EE5B546D255AEA89A00CA8C8015063F3A98600E157AFBC19AB73BD38F1A5A81EE4B8F923AE3DF3629DC10661C6C81D2&originRegion=us-east-1&originCreation=20211027184033
(accessed Oct. 27, 2021).
    """
    def __init__(self, *args, **kwargs):
        super(ThreeMethylHeptane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path(
            '3-methylheptane-c8h18-3-1378-8143.yaml')
        self.fuel = 'c8h18-3:1.0'
        # default temperature and pressure for a fuel
        self.options["T"] = kwargs.get("T", 1300)
        self.options["P"] = kwargs.get("P", ct.one_atm)


class NHexadecane(ModelBase):
    """
    LNLL n-hexadecane
    """
    def __init__(self, *args, **kwargs):
        super(NHexadecane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('nhexadecane-2115-13341.yaml')
        self.fuel = 'nc16h34:1.0'
        # default temperature and pressure for a fuel
        self.options["T"] = kwargs.get("T", 1300)
        self.options["P"] = kwargs.get("P", ct.one_atm)

class C8_C18_Blends(ModelBase):
    """
    LNLL n-hexadecane
    """
    def __init__(self, *args, **kwargs):
        super(C8_C18_Blends, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('nhexadecane-2115-13341.yaml')
        self.fuel = 'nc16h34:1.0, nc8h18:1.0, nc12h26:1.0, nc10h22:1.0' # NHexadecane, NOctane, NDodecane
        # default temperature and pressure for a fuel
        self.options["T"] = kwargs.get("T", 1300)
        self.options["P"] = kwargs.get("P", ct.one_atm)


class MethylFiveDeconate(ModelBase):
    """
    LNLL methyl-5-deconate
    """
    def __init__(self, *args, **kwargs):
        super(MethylFiveDeconate, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('md5d-2649-10487.yaml')
        self.fuel = 'md5d:1.0'
        # default temperature and pressure for a fuel
        self.options["T"] = kwargs.get("T", 1250)
        self.options["P"] = kwargs.get("P", ct.one_atm)


class IsoOctaneDetailed(ModelBase):  # Iso-Octane
    def __init__(self, *args, **kwargs):
        super(IsoOctaneDetailed, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('ic8-detailed-2768-11850.yaml')
        self.fuel = 'IC8H18:1.0'
        # default temperature and pressure for a fuel
        self.options["T"] = kwargs.get("T", 1250)
        self.options["P"] = kwargs.get("P", ct.one_atm)

class MethylNineDeconate(ModelBase):
    """
    LNLL methyl-9-deconate
    """
    def __init__(self, *args, **kwargs):
        super(MethylNineDeconate, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('md9d-3298.yaml')
        self.fuel = 'md9d:1.0'
        # default temperature and pressure for a fuel
        self.options["T"] = kwargs.get("T", 1250)
        self.options["P"] = kwargs.get("P", ct.one_atm)


class MethylDeconateNHeptane(ModelBase):
    """
    LNLL methyldeconate with nheptane
    """

    def __init__(self, *args, **kwargs):
        super(MethylDeconateNHeptane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('md-nc7-3787-10264.yaml')
        self.fuel = 'md:1.0, nc7h16:1.0'
        # default temperature and pressure for a fuel
        self.options["T"] = kwargs.get("T", 1250)
        self.options["P"] = kwargs.get("P", ct.one_atm)


class TwoMethylnonadecane(ModelBase):
    """
    S.M. Sarathy, M. Mehl, C. K. Westbrook, W. J. Pitz,C. Togbe, P.
    Dagaut, H. Wang, M.A. Oehlschlaeger, U. Niemann, K. Seshadri,P.S.
    Veloo, C. Ji, F.N. Egolfopoulos, T. Lu Comprehensive chemical
    kinetic modeling of the oxidation of 2-methylalkanes from C7 to C20
    Combustion and Flame 2011
    """
    def __init__(self, *args, **kwargs):
        super(TwoMethylnonadecane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('mmc5-7171-38324.yaml')
        self.fuel = 'c20h42-2:1.0'
        # default temperature and pressure for a fuel
        self.options["T"] = kwargs.get("T", 1250)
        self.options["P"] = kwargs.get("P", ct.one_atm)

