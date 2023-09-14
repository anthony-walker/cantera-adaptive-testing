import math

def KRO2NO(T):
    return 2.7e-12 * math.exp(360 / T)

def KRO2HO2(T):
    return 2.91e-13 * math.exp(1300 / T)

def KAPHO2(T):
    return 5.2e-13 * math.exp(980 / T)

def KAPNO(T):
    return 7.5e-12 * math.exp(290 / T)

def KRO2NO3(T):
    return 2.3e-12

def KNO3AL(T):
    return 1.4e-12 * math.exp(-1860 / T)

def KDEC(T):
    return 1.0e6

def KROPRIM(T):
    return 2.5e-14 * math.exp(-300 / T)

def KROSEC(T):
    return 2.5e-14 * math.exp(-300 / T)

def KCH3O2(T):
    return 1.03e-13 * math.exp(365 / T)

def K298CH3O2(T):
    return 3.5e-13

def K14ISOM1(T):
    return 3.00e7 * math.exp(-5300 / T)

def KD0(T, M):
    return 1.10e-05 * M * math.exp(-10100 / T)

def KDI(T):
    return 1.90e17 * math.exp(-14100 / T)

def KRD(T, M):
    return KD0(T, M) / KDI(T)

def FCD(T):
    return 0.30

def NCD(T):
    return 0.75 - 1.27 * (math.log10(FCD(T)))

def FD(T):
    return 10 ** (math.log10(FCD(T)) / (1 + (math.log10(KRD(T)) / NCD(T)) ** 2))

def KBPAN(T, M):
    return (KD0(T, M) * KDI(T) * FD(T)) / (KD0(T, M) + KDI(T))

def KC0(T, M):
    return 3.28e-28 * M * (T / 300) ** -6.87

def KCI(T):
    return 1.125e-11 * (T / 300) ** -1.105

def KRC(T, M):
    return KC0(T, M) / KCI(T)

def FCC(T):
    return 0.30

def NC(T):
    return 0.75 - 1.27 * (math.log10(FCC(T)))

def FC(T):
    return 10 ** (math.log10(FCC(T)) / (1 + (math.log10(KRC(T)) / NC(T)) ** 2))

def KFPAN(T, M):
    return (KC0(T, M) * KCI(T) * FC(T)) / (KC0(T, M) + KCI(T))

def K10(T, M):
    return 1.0e-31 * M * (T / 300) ** -1.6

def K1I(T):
    return 5.0e-11 * (T / 300) ** -0.3

def KR1(T, M):
    return K10(T, M) / K1I(T)

def FC1(T):
    return 0.85

def NC1(T):
    return 0.75 - 1.27 * (math.log10(FC1(T)))

def F1(T):
    return 10 ** (math.log10(FC1(T)) / (1 + (math.log10(KR1(T)) / NC1(T)) ** 2))

def KMT01(T, M):
    return (K10(T, M) * K1I(T) * F1(T)) / (K10(T, M) + K1I(T))

def K20(T, M):
    return 1.3e-31 * M * (T / 300) ** -1.5

def K2I(T):
    return 2.3e-11 * (T / 300) ** 0.24

def KR2(T, M):
    return K20(T, M) / K2I(T)

def FC2(T):
    return 0.6

def NC2(T):
    return 0.75 - 1.27 * (math.log10(FC2(T)))

def F2(T):
    return 10 ** (math.log10(FC2(T)) / (1 + (math.log10(KR2(T)) / NC2(T)) ** 2))

def KMT02(T, M):
    return (K20(T, M) * K2I(T) * F2(T)) / (K20(T, M) + K2I(T))

def K30(T, M):
    return 3.6e-30 * M * (T / 300) ** -4.1

def K3I(T):
    return 1.9e-12 * (T / 300) ** 0.2

def KR3(T, M):
    return K30(T, M) / K3I(T)

def FC3(T):
    return 0.35

def NC3(T):
    return 0.75 - 1.27 * (math.log10(FC3(T)))

def F3(T):
    return 10 ** (math.log10(FC3(T)) / (1 + (math.log10(KR3(T)) / NC3(T)) ** 2))

def KMT03(T, M):
    return (K30(T, M) * K3I(T) * F3(T)) / (K30(T, M) + K3I(T))

def K40(T, M):
    return 1.3e-3 * M * (T / 300) ** -3.5 * math.exp(-11000 / T)

def K4I(T):
    return 9.7e14 * (T / 300) ** 0.1 * math.exp(-11080 / T)

def KR4(T, M):
    return K40(T, M) / K4I(T)

def FC4(T):
    return 0.35

def NC4(T):
    return 0.75 - 1.27 * (math.log10(FC4(T)))

def F4(T):
    return 10 ** (math.log10(FC4(T)) / (1 + (math.log10(KR4(T)) / NC4(T)) ** 2))

def KMT04(T, M):
    return (K40(T, M) * K4I(T) * F4(T)) / (K40(T, M) + K4I(T))

def KMT05(M):
    return 1.44e-13 * (1 + (M / 4.2e19))

def KMT06(T, H2O):
    return 1 + (1.40e-21 * math.exp(2200 / T) * H2O)

def K70(T, M):
    return 7.4e-31 * M * (T / 300) ** -2.4

def K7I(T):
    return 3.3e-11 * (T / 300) ** -0.3

def KR7(T, M):
    return K70(T, M) / K7I(T)

def FC7(T):
    return 0.81

def NC7(T):
    return 0.75 - 1.27 * (math.log10(FC7(T)))

def F7(T):
    return 10 ** (math.log10(FC7(T)) / (1 + (math.log10(KR7(T)) / NC7(T)) ** 2))

def KMT07(T, M):
    return (K70(T, M) * K7I(T) * F7(T)) / (K70(T, M) + K7I(T))

def K80(T, M):
    return 3.2e-30 * M * (T / 300) ** -4.5

def K8I(T):
    return 3.0e-11

def KR8(T, M):
    return K80(T, M) / K8I(T)

def FC8(T):
    return 0.41

def NC8(T):
    return 0.75 - 1.27 * (math.log10(FC8(T)))

def F8(T):
    return 10 ** (math.log10(FC8(T)) / (1 + (math.log10(KR8(T)) / NC8(T)) ** 2))

def KMT08(T, M):
    return (K80(T, M) * K8I(T) * F8(T)) / (K80(T, M) + K8I(T))

def K90(T, M):
    return 1.4e-31 * M * (T / 300) ** -3.1

def K9I(T):
    return 4.0e-12

def KR9(T, M):
    return K90(T, M) / K9I(T)

def FC9(T):
    return 0.4

def NC9(T):
    return 0.75 - 1.27 * (math.log10(FC9(T)))

def F9(T):
    return 10 ** (math.log10(FC9(T)) / (1 + (math.log10(KR9(T)) / NC9(T)) ** 2))

def KMT09(T, M):
    return (K90(T, M) * K9I(T) * F9(T)) / (K90(T, M) + K9I(T))

def K100(T, M):
    return 4.10e-05 * M * math.exp(-10650 / T)

def K10I(T):
    return 6.0e15 * math.exp(-11170 / T)

def KR10(T, M):
    return K100(T, M) / K10I(T)

def FC10(T):
    return 0.4

def NC10(T):
    return 0.75 - 1.27 * (math.log10(FC10(T)))

def F10(T):
    return 10 ** (math.log10(FC10(T)) / (1 + (math.log10(KR10(T)) / NC10(T)) ** 2))

def KMT10(T, M):
    return (K100(T, M) * K10I(T) * F10(T)) / (K100(T, M) + K10I(T))

def K1(T):
    return 2.40e-14 * math.exp(460 / T)

def K3(T):
    return 6.50e-34 * math.exp(1335 / T)

def K4(T):
    return 2.70e-17 * math.exp(2199 / T)

def K2(T, M):
    return (K3(T) * M) / (1 + (K3(T) * M) / K4(T))

def KMT11(T, M):
    return K1(T) + K2(T, M)

def K120(T, M):
    return 2.5e-31 * M * (T / 300) ** -2.6

def K12I(T):
    return 2.0e-12

def KR12(T, M):
    return K120(T, M) / K12I(T)

def FC12(T):
    return 0.53

def NC12(T):
    return 0.75 - 1.27 * (math.log10(FC12(T)))

def F12(T):
    return 10 ** (math.log10(FC12(T)) / (1 + (math.log10(KR12(T)) / NC12(T)) ** 2))

def KMT12(T, M):
    return (K120(T, M) * K12I(T) * F12(T)) / (K120(T, M) + K12I(T))

def K130(T, M):
    return 2.5e-30 * M * (T / 300) ** -5.5

def K13I(T):
    return 1.8e-11

def KR13(T, M):
    return K130(T, M) / K13I(T)

def FC13(T):
    return 0.36

def NC13(T):
    return 0.75 - 1.27 * (math.log10(FC13(T)))

def F13(T):
    return 10 ** (math.log10(FC13(T)) / (1 + (math.log10(KR13(T)) / NC13(T)) ** 2))

def KMT13(T, M):
    return (K130(T, M) * K13I(T) * F13(T)) / (K130(T, M) + K13I(T))

def K140(T, M):
    return 9.0e-5 * math.exp(-9690 / T) * M

def K14I(T):
    return 1.1e16 * math.exp(-10560 / T)

def KR14(T, M):
    return K140(T, M) / K14I(T)

def FC14(T):
    return 0.36

def NC14(T):
    return 0.75 - 1.27 * (math.log10(FC14(T)))

def F14(T):
    return 10 ** (math.log10(FC14(T)) / (1 + (math.log10(KR14(T)) / NC14(T)) ** 2))

def KMT14(T, M):
    return (K140(T, M) * K14I(T) * F14(T)) / (K140(T, M) + K14I(T))

def K150(T, M):
    return 8.6e-29 * M * (T / 300) ** -3.1

def K15I(T):
    return 9.0e-12 * (T / 300) ** -0.85

def KR15(T, M):
    return K150(T, M) / K15I(T)

def FC15(T):
    return 0.48

def NC15(T):
    return 0.75 - 1.27 * (math.log10(FC15(T)))

def F15(T):
    return 10 ** (math.log10(FC15(T)) / (1 + (math.log10(KR15(T)) / NC15(T)) ** 2))

def KMT15(T, M):
    return (K150(T, M) * K15I(T) * F15(T)) / (K150(T, M) + K15I(T))

def K160(T, M):
    return 8.0e-27 * M * (T / 300) ** -3.5

def K16I(T):
    return 3.0e-11 * (T / 300) ** -1

def KR16(T, M):
    return K160(T, M) / K16I(T)

def FC16(T):
    return 0.5

def NC16(T):
    return 0.75 - 1.27 * (math.log10(FC16(T)))

def F16(T):
    return 10 ** (math.log10(FC16(T)) / (1 + (math.log10(KR16(T)) / NC16(T)) ** 2))

def KMT16(T, M):
    return (K160(T, M) * K16I(T) * F16(T)) / (K160(T, M) + K16I(T))

def K170(T, M):
    return 5.0e-30 * M * (T / 300) ** -1.5

def K17I(T):
    return 1.0e-12

def KR17(T, M):
    return K170(T, M) / K17I(T)

def FC17(T):
    return 0.17 * math.exp(-51 / T) + math.exp(-T / 204)

def NC17(T):
    return 0.75 - 1.27 * (math.log10(FC17(T)))

def F17(T):
    return 10 ** (math.log10(FC17(T)) / (1 + (math.log10(KR17(T)) / NC17(T)) ** 2))

def KMT17(T, M):
    return (K170(T, M) * K17I(T) * F17(T)) / (K170(T, M) + K17I(T))

def KMT18(T, O2):
    return 9.5e-39 * O2 * math.exp(5270 / T) / (1 + 7.5e-29 * O2 * math.exp(5610 / T))

def KPPN0(T, M):
    return 1.7e-03 * math.exp(-11280 / T) * M

def KPPNI(T):
    return 8.3e16 * math.exp(-13940 / T)

def KRPPN(T, M):
    return KPPN0(T, M) / KPPNI(T)

def FCPPN(T):
    return 0.36

def NCPPN(T):
    return 0.75 - 1.27 * (math.log10(FCPPN(T)))

def FPPN(T):
    return 10 ** (math.log10(FCPPN(T)) / (1 + (math.log10(KRPPN(T)) / NCPPN(T)) ** 2))

def KBPPN(T, M):
    return (KPPN0(T, M) * KPPNI(T) * FPPN(T)) / (KPPN0(T, M) + KPPNI(T))

def KAWQ1(T):
    return KCH3O2() * (1-7.18*math.exp(-885/T))

def KAWQ2(T):
    return 8.8e-12*math.exp(-1320/T) + 1.7e-14*math.exp(423/T)

def KAWQ3(T):
    return (1-1/(1+498*math.exp(-1160/T)))

def KAWQ4(T):
    return (1/(1+498*math.exp(-1160/T)))

def KAWQ5(T):
    return (1-math.exp(-550/T))
