"""Metadata for RFMIP species.
"name" and "ordinary_formula" are based on the Hitran molecule headers.
See https://hitran.org/suppl/xsec/molecule_headers/
"""
XSEC_SPECIES_INFO = {
    # Alcohols, ethers and other oxygenated hydrocarbons
    # Bromocarbons, Hydrobromocarbons, and Halons
    "Halon-1211": {
        "name": "Bromochlorodifluoromethane",
        "ordinary_formula": "CBrClF2",
        "rfmip": "halon1211_GM",
    },
    "Halon-1301": {
        "name": "Bromotrifluoromethane",
        "ordinary_formula": "CBrF3",
        "rfmip": "halon1301_GM",
    },
    "Halon-2402": {
        "name": "1,2-Dibromo-1,1,2,2-tetrafluoroethane",
        "ordinary_formula": "CBrF2CBrF2",
        "rfmip": "halon2402_GM",
    },
    # Chlorocarbons and Hydrochlorocarbons
    "CCl4": {
        "name": "Carbon Tetrachloride",
        "ordinary_formula": "CCl4",
        "rfmip": "carbon_tetrachloride_GM",
    },
    "CH2Cl2": {
        "name": "Dichloromethane",
        "ordinary_formula": "CH2Cl2",
        "rfmip": "ch2cl2_GM",
    },
    "CH3CCl3": {
        "name": "1,1,1-Trichloroethane",
        "ordinary_formula": "CH3CCl3",
        "rfmip": "ch3ccl3_GM",
    },
    "CHCl3": {
        "name": "Trichloromethane",
        "ordinary_formula": "CHCl3",
        "rfmip": "chcl3_GM",
    },
    # Chlorofluorocarbons (CFCs)
    "CFC-11": {
        "name": "CFC-11",
        "ordinary_formula": "CCl3F",
        "rfmip": "cfc11_GM",
    },
    "CFC-12": {
        "name": "CFC-12",
        "ordinary_formula": "CCl2F2",
        "rfmip": "cfc12_GM",
    },
    "CFC-113": {
        "name": "CFC-113",
        "ordinary_formula": "C2Cl3F3",
        "rfmip": "cfc113_GM",
    },
    "CFC-114": {
        "name": "CFC-114",
        "ordinary_formula": "C2Cl2F4",
        "rfmip": "cfc114_GM",
    },
    "CFC-115": {
        "name": "CFC-115",
        "ordinary_formula": "C2ClF5",
        "rfmip": "cfc115_GM",
    },
    # Fully Fluorinated Species
    "C2F6": {
        "name": "Hexafluoroethane",
        "ordinary_formula": "C2F6",
        "rfmip": "c2f6_GM",
    },
    "C3F8": {
        "name": "Perfluoropropane",
        "ordinary_formula": "C3F8",
        "rfmip": "c3f8_GM",
    },
    "C4F10": {
        "name": "Perfluorobutane",
        "ordinary_formula": "C4F10",
        "rfmip": "c4f10_GM",
    },
    "C5F12": {
        "name": "Perfluoropentane",
        "ordinary_formula": "n-C5F12",
        "rfmip": "c5f12_GM",
    },
    "C6F14": {
        "name": "Perfluorohexane",
        "ordinary_formula": "n-C6F14",
        "rfmip": "c6f14_GM",
    },
    "C8F18": {
        "name": "Perfluorooctane",
        "ordinary_formula": "C8F18",
        "rfmip": "c8f18_GM",
    },
    "c-C4F8": {
        "name": "Perfluorocyclobutane",
        "ordinary_formula": "c-C4F8",
        "rfmip": "c_c4f8_GM",
    },
    "CF4": {
        "name": "Carbon Tetrafluoride",
        "ordinary_formula": "CF4",
        "rfmip": "cf4_GM",
    },
    "NF3": {
        "name": "Nitrogen trifluoride",
        "ordinary_formula": "NF3",
        "rfmip": "nf3_GM",
    },
    "SF6": {
        "name": "Sulfur Hexafluoride",
        "ordinary_formula": "SF6",
        "rfmip": "sf6_GM",
    },
    "SO2F2": {
        "name": "Sulfuryl fluoride",
        "ordinary_formula": "SO2F2",
        "rfmip": "so2f2_GM",
    },
    # Halogenated Alcohols and Ethers
    # Hydrocarbons
    # Hydrochlorofluorocarbons (HCFCs)
    "HCFC-141b": {
        "name": "HCFC-141b",
        "ordinary_formula": "CH3CCl2F",
        "rfmip": "hcfc141b_GM",
    },
    "HCFC-142b": {
        "name": "HCFC-142b",
        "ordinary_formula": "CH3CClF2",
        "rfmip": "hcfc142b_GM",
    },
    "HCFC-22": {
        "name": "HCFC-22",
        "ordinary_formula": "CHClF2",
        "rfmip": "hcfc22_GM",
    },
    # Hydrofluorocarbons (HFCs)
    "HFC-125": {
        "name": "HFC-125",
        "ordinary_formula": "CHF2CF3",
        "rfmip": "hfc125_GM",
    },
    "HFC-134a": {
        "name": "HFC-134a",
        "ordinary_formula": "CFH2CF3",
        "rfmip": "hfc134a_GM",
    },
    "HFC-143a": {
        "name": "HFC-143a",
        "ordinary_formula": "CF3CH3",
        "rfmip": "hfc143a_GM",
    },
    "HFC-152a": {
        "name": "HFC-152a",
        "ordinary_formula": "CH3CHF2",
        "rfmip": "hfc152a_GM",
    },
    "HFC-227ea": {
        "name": "1,1,1,2,3,3,3-Heptafluoropropane",
        "ordinary_formula": "CF3CHFCF3",
        "rfmip": "hfc227ea_GM",
    },
    "HFC-236fa": {
        "name": "1,1,1,3,3,3-Hexafluoropropane",
        "ordinary_formula": "CF3CH2CF3",
        "rfmip": "hfc236fa_GM",
    },
    "HFC-23": {
        "name": "Trifluoromethane",
        "ordinary_formula": "CHF3",
        "rfmip": "hfc23_GM",
    },
    "HFC-245fa": {
        "name": "1,1,1,3,3-Pentafluoropropane",
        "ordinary_formula": "CHF2CH2CF3",
        "rfmip": "hfc245fa_GM",
    },
    "HFC-32": {
        "name": "HFC-32",
        "ordinary_formula": "CH2F2",
        "rfmip": "hfc32_GM",
    },
    "HFC-365mfc": {
        "name": "1,1,1,3,3-Pentafluorobutane",
        "ordinary_formula": "CH3CF2CH2CF3",
        "rfmip": "hfc365mfc_GM",
    },
    "HFC-43-10mee": {
        "name": "1,1,1,2,2,3,4,5,5,5-Decafluoropentane",
        "ordinary_formula": "CF3CHFCHFCF2CF3",
        "altname": "CF3CHFCHFCF2CF3",
        "rfmip": "hfc4310mee_GM",
    },
    # Iodocarbons and hydroiodocarbons
    # Nitriles, amines and other nitrogenated hydrocarbons
    # Other molecules
    # Sulfur-containing species
}

SPECIES_GROUPS = {
    "reference": [
        "CCl4",
        "CF4",
        "CFC-11",
        "CFC-12",
        "HFC-134a",
        "HFC-23",
    ],
    "rfmip-names": [
        "c2f6_GM",
        "c3f8_GM",
        "c4f10_GM",
        "c5f12_GM",
        "c6f14_GM",
        "c7f16_GM",
        "c8f18_GM",
        "c_c4f8_GM",
        "carbon_dioxide_GM",
        "carbon_tetrachloride_GM",
        "cf4_GM",
        "cfc113_GM",
        "cfc114_GM",
        "cfc115_GM",
        "cfc11_GM",
        "cfc11eq_GM",
        "cfc12_GM",
        "cfc12eq_GM",
        "ch2cl2_GM",
        "ch3ccl3_GM",
        "chcl3_GM",
        "halon1211_GM",
        "halon1301_GM",
        "halon2402_GM",
        "hcfc141b_GM",
        "hcfc142b_GM",
        "hcfc22_GM",
        "hfc125_GM",
        "hfc134a_GM",
        "hfc134aeq_GM",
        "hfc143a_GM",
        "hfc152a_GM",
        "hfc227ea_GM",
        "hfc236fa_GM",
        "hfc23_GM",
        "hfc245fa_GM",
        "hfc32_GM",
        "hfc365mfc_GM",
        "hfc4310mee_GM",
        "methane_GM",
        "methyl_bromide_GM",
        "methyl_chloride_GM",
        "nf3_GM",
        "sf6_GM",
        "so2f2_GM",
    ],
}

RFMIPMAP = {v["rfmip"]: k for k, v in XSEC_SPECIES_INFO.items() if "rfmip" in v.keys()}
SPECIES_GROUPS["rfmip"] = [
    RFMIPMAP[k] for k in SPECIES_GROUPS["rfmip-names"] if k in RFMIPMAP.keys()
]
