from .model_base import ModelBase
import cantera as ct


class GAS_DEF(ModelBase):
    def __init__(self, *args, **kwargs):
        super(GAS_DEF, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('gas-def.yaml')
        self.fuel = 'CH4:1.0'
        self.skip_database_build = True

class Hydrogen(ModelBase):
    def __init__(self, *args, **kwargs):
        super(Hydrogen, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('hydrogen-10-28.yaml')
        self.fuel = 'H2:1.0'
        self.skip_database_build = False

class PlatinumSmallHydrogen(ModelBase):
    def __init__(self, *args, **kwargs):
        super(PlatinumSmallHydrogen, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('hydrogen-10-28.yaml')
        self.fuel = 'H2:1.0'
        self.surface = 'PT(S):10.0'
        self.sphase = 'surface-small'
        self.skip_database_build = False

class PlatinumMediumHydrogen(ModelBase):
    def __init__(self, *args, **kwargs):
        super(PlatinumMediumHydrogen, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('hydrogen-10-28.yaml')
        self.fuel = 'H2:1.0'
        self.surface = 'Pt(9):10.0'
        self.sphase = 'surface-medium'
        self.skip_database_build = False

class PlatinumLargeHydrogen(ModelBase):
    def __init__(self, *args, **kwargs):
        super(PlatinumLargeHydrogen, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('hydrogen-10-28.yaml')
        self.fuel = 'H2:1.0'
        self.surface = 'Pt(9):10.0'
        self.sphase = 'surface-large'
        self.skip_database_build = False

class PlatinumSmallAramco(ModelBase):  # Hydrogen with more species
    def __init__(self, *args, **kwargs):
        super(PlatinumSmallAramco, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('aramco-493-2716.yaml')
        self.fuel = 'CH4:1.0, C3H8:1.0, C2H6:1.0'
        self.surface = 'PT(S):10.0'
        self.sphase = 'surface-small'

class PlatinumMediumAramco(ModelBase):  # Hydrogen with more species
    def __init__(self, *args, **kwargs):
        super(PlatinumMediumAramco, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('aramco-493-2716.yaml')
        self.fuel = 'CH4:1.0, C3H8:1.0, C2H6:1.0'
        self.surface = 'Pt(9):10.0'
        self.sphase = 'surface-medium'
        self.skip_database_build = False

class PlatinumLargeAramco(ModelBase):  # Hydrogen with more species
    def __init__(self, *args, **kwargs):
        super(PlatinumLargeAramco, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('aramco-493-2716.yaml')
        self.fuel = 'CH4:1.0, C3H8:1.0, C2H6:1.0'
        self.surface = 'Pt(9):10.0'
        self.sphase = 'surface-large'
        self.skip_database_build = False

class PlatinumSmallGRI(ModelBase):  # Hydrogen with more species
    def __init__(self, *args, **kwargs):
        super(PlatinumSmallGRI, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('gri-mech-55-325.yaml')
        self.fuel = 'CH4:1.0'
        self.surface = 'PT(S):10.0'
        self.gphase = 'gas'
        self.sphase = 'surface-small'
        self.skip_database_build = False

class PlatinumMediumGRI(ModelBase):  # Hydrogen with more species
    def __init__(self, *args, **kwargs):
        super(PlatinumMediumGRI, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('gri-mech-55-325.yaml')
        self.fuel = 'CH4:1.0'
        self.surface = 'Pt(9):10.0'
        self.sphase = 'surface-medium'
        self.skip_database_build = False

class PlatinumLargeGRI(ModelBase):  # Hydrogen with more species
    def __init__(self, *args, **kwargs):
        super(PlatinumLargeGRI, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('gri-mech-55-325.yaml')
        self.fuel = 'CH4:1.0'
        self.surface = 'Pt(9):10.0'
        self.sphase = 'surface-large'
        self.skip_database_build = False

class DME(ModelBase):
    def __init__(self, *args, **kwargs):
        super(DME, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('dme-propane-122-711.yaml')
        self.fuel = 'CH3OCH3:1.0, C3H8:1.0'
        self.skip_database_build = False


class JetA(ModelBase):
    def __init__(self, *args, **kwargs):
        super(JetA, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('jetA-detailed-NOx-203-1589.yaml')
        self.fuel = 'POSF10325:1.0'
        self.skip_database_build = False


class Butane(ModelBase):
    def __init__(self, *args, **kwargs):
        super(Butane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('butane-230-2461.yaml')
        self.fuel = 'C4H10:1.0'
        self.skip_database_build = False


class TwoButonane(ModelBase):
    def __init__(self, *args, **kwargs):
        super(TwoButonane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path(
            '2-butonane-ch3coch2ch3-315-1803.yaml')
        self.fuel = 'C4H8O1-2:1.0'
        self.skip_database_build = False


class NHexanal(ModelBase):
    def __init__(self, *args, **kwargs):
        super(NHexanal, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('n-hexanal-482-5182.yaml')
        self.fuel = 'C6H14:1.0'
        self.skip_database_build = False


class IsoButene(ModelBase):
    def __init__(self, *args, **kwargs):
        super(IsoButene, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('isobutene-ic4h8-493-2716.yaml')
        self.fuel = 'IC4H8:1.0'
        self.skip_database_build = False


class IsoPentanol(ModelBase):
    def __init__(self, *args, **kwargs):
        super(IsoPentanol, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('i-pentanol-511-3010.yaml')
        self.fuel = 'ic5h10oh:1.0'
        self.skip_database_build = False


class T124MCH(ModelBase):
    def __init__(self, *args, **kwargs):
        super(T124MCH, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('t124mch-533-3193.yaml')
        self.fuel = 'T124MCH:1.0'
        self.skip_database_build = False


class OneTwoDME(ModelBase):
    def __init__(self, *args, **kwargs):
        super(OneTwoDME, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('one-two-dme-570-2960.yaml')
        self.fuel = 'CH3OCH3:1.0'
        self.skip_database_build = False


class PropylAcetate(ModelBase):
    def __init__(self, *args, **kwargs):
        super(PropylAcetate, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('propyl-acetate-628-4182.yaml')
        self.fuel = 'pa:1.0'
        self.skip_database_build = False


class NHeptane(ModelBase):
    def __init__(self, *args, **kwargs):
        super(NHeptane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('n-heptane-c7h16-654-4846.yaml')
        self.fuel = 'NC7H16:1.0'
        self.skip_database_build = False


class DEE(ModelBase):
    def __init__(self, *args, **kwargs):
        super(DEE, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('dee-746-3555.yaml')
        self.fuel = 'DEE:1.0'
        self.skip_database_build = False


class IsoOctane(ModelBase):  # Iso-Octane
    def __init__(self, *args, **kwargs):
        super(IsoOctane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('ic8-874-6864.yaml')
        self.fuel = 'IC8H18:1.0'
        self.skip_database_build = False


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
        self.skip_database_build = False


class NHexadecane(ModelBase):
    """
    LNLL n-hexadecane
    """
    def __init__(self, *args, **kwargs):
        super(NHexadecane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('nhexadecane-2115-13341.yaml')
        self.fuel = 'nc16h34:1.0'
        self.skip_database_build = False


class MethylFiveDeconate(ModelBase):
    """
    LNLL methyl-5-deconate
    """
    def __init__(self, *args, **kwargs):
        super(MethylFiveDeconate, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('md5d-2649-10487.yaml')
        self.fuel = 'md5d:1.0'
        self.skip_database_build = False


class IsoOctaneDetailed(ModelBase):  # Iso-Octane
    def __init__(self, *args, **kwargs):
        super(IsoOctaneDetailed, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('ic8-detailed-2768-11850.yaml')
        self.fuel = 'IC8H18:1.0'
        self.skip_database_build = False


class MethylNineDeconate(ModelBase):
    """
    LNLL methyl-9-deconate
    """
    def __init__(self, *args, **kwargs):
        super(MethylNineDeconate, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('md9d-3298.yaml')
        self.fuel = 'md9d:1.0'
        self.skip_database_build = False


class MethylDeconateNHeptane(ModelBase):
    """
    LNLL methyldeconate with nheptane
    """

    def __init__(self, *args, **kwargs):
        super(MethylDeconateNHeptane, self).__init__(*args, **kwargs)
        self.model = self.get_test_set_path('md-nc7-3787-10264.yaml')
        self.fuel = 'md:1.0, nc7h16:1.0'
        self.skip_database_build = False


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
        self.skip_database_build = False
