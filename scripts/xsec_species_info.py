"""Metadata for RFMIP species."""

XSEC_SPECIES_INFO = {
    # Alcohols, ethers and other oxygenated hydrocarbons
    # Bromocarbons, Hydrobromocarbons, and Halons
    "Halon-1211": {  # no air broadening, only N2
        "altname": "CBrClF2",
        "rfmip": "halon1211_GM",
        "reftemp": 300,
    },
    "Halon-1301": {  # no air broadening, only N2
        "altname": "CBrF3",
        "rfmip": "halon1301_GM",
        "reftemp": 300,
    },
    "Halon-2402": {"altname": "CBrF2CBrF2", "rfmip": "halon2402_GM", "reftemp": 300,},
    # Chlorocarbons and Hydrochlorocarbons
    "CCl4": {  # +++++ fit ok +++++, use only band 700-860
        "active": True,
        "arts_bands": ((700, 860),),
        "rfmip": "carbon_tetrachloride_GM",
        "reftemp": 250,
    },
    "CH2Cl2": {"rfmip": "ch2cl2_GM", "reftemp": 300,},  # no air broadening, only N2
    "CH3CCl3": {"rfmip": "ch3ccl3_GM", "reftemp": 300,},  # no air broadening, only N2
    "CHCl3": {  # not available in Hitran
        "arts_bands": ((580, 7200),),
        "rfmip": "chcl3_GM",
        "reftemp": 300,
    },
    # Chlorofluorocarbons (CFCs)
    "CFC-11": {"active": True, "rfmip": "cfc11_GM",
               #"reftemp": 220,
               },  # +++++ fit ok +++++
    "CFC-12": {"active": True, "rfmip": "cfc12_GM",
               #"reftemp": 220,
               },  # +++++ fit ok +++++
    "CFC-113": {"active": True, "rfmip": "cfc113_GM",},  # only data for 0 torr
    "CFC-114": {"rfmip": "cfc114_GM",},  # only data for 0 torr
    "CFC-115": {"rfmip": "cfc115_GM",},  # only data for 0 torr
    # Fully Fluorinated Species
    "C2F6": {  # !!!!! bad fit !!!!! no high pressure differences available
        "use_average": True,
        "arts_bands": ((500, 6500),),
        "rfmip": "c2f6_GM",
        "reftemp": 300,
    },
    "C3F8": {"rfmip": "c3f8_GM", "reftemp": 300,},  # no air broadening, only N2
    "C4F10": {"rfmip": "c4f10_GM", "reftemp": 300,},  # no air broadening, only N2
    "C5F12": {  # no air broadening, only N2
        "altname": "n-C5F12",
        "arts_bands": ((500, 6500),),
        "rfmip": "c5f12_GM",
        "reftemp": 300,
    },
    "C6F14": {"altname": "n-C6F14", "rfmip": "c6f14_GM",},  # no air broadening, only N2
    "C8F18": {"rfmip": "c8f18_GM",},  # no air broadening, only N2 at 0 Torr
    "c-C4F8": {  # only data for 0 Torr
        "arts_bands": ((550, 6500),),
        "rfmip": "c_c4f8_GM",
        "reftemp": 300,
    },
    "CF4": {  # +++++ fit ok +++++
        "active": True,
        "arts_bands": ((1250, 1290),),
        "rfmip": "cf4_GM",
        #"reftemp": 240,
    },
    "NF3": {"rfmip": "nf3_GM", "reftemp": 300,},  # no air broadening, only N2
    "SF6": {
        "use_average": True,
        "arts_bands": ((560, 6500),),
        "rfmip": "sf6_GM",
        "reftemp": 300,
    },
    "SO2F2": {"rfmip": "so2f2_GM", "reftemp": 300,},  # no air broadening, only N2
    # Halogenated Alcohols and Ethers
    # Hydrocarbons
    # Hydrochlorofluorocarbons (HCFCs)
    "HCFC-141b": {"rfmip": "hcfc141b_GM",},  # only data for 0 torr
    "HCFC-142b": {  # only data for 0 torr
        "arts_bands": ((650, 1500),),
        "rfmip": "hcfc142b_GM",
    },
    "HCFC-22": {  # !!!!! bad fit !!!!! no high pressure differences available
        "use_average": True,
        "arts_bands": ((760, 860), (1060, 1210), (1275, 1380),),
        "rfmip": "hcfc22_GM",
        "reftemp": 210,
    },
    # Hydrofluorocarbons (HFCs)
    "HFC-125": {
        "use_average": True,
        "arts_bands": ((495, 1504),),
        "rfmip": "hfc125_GM",
    },
    "HFC-134a": {  # +++++ fit ok +++++. Use band 750-1600.
        "altname": "CFH2CF3",
        "active": True,
        "arts_bands": ((750, 1600),),
        "rfmip": "hfc134a_GM",
        "reftemp": 220,
    },
    "HFC-143a": {  # not enough xsecs available
        "use_average": True,
        "arts_bands": ((580, 630), (694, 1504),),
        "rfmip": "hfc143a_GM",
    },
    "HFC-152a": {"rfmip": "hfc152a_GM",},
    "HFC-227ea": {
        "altname": "CF3CHFCF3",
        "arts_bands": ((500, 6500),),  # Check first band for overlap
        "rfmip": "hfc227ea_GM",
    },
    "HFC-236fa": {  # No usable input files
        "altname": "CF3CH2CF3",
        "rfmip": "hfc236fa_GM",
    },
    "HFC-23": {"altname": "CHF3", "rfmip": "hfc23_GM", "reftemp": 230,},
    "HFC-245fa": {"altname": "CHF2CH2CF3", "rfmip": "hfc245fa_GM",},  # Only one profile
    "HFC-32": {  # !!!!! bad fit !!!!! not enough xsecs available
        "use_average": True,
        "rfmip": "hfc32_GM",
    },
    "HFC-365mfc": {  # Only one profile ("air")
        "altname": "CH3CF2CH2CF3",
        "rfmip": "hfc365mfc_GM",
    },
    "HFC-43-10mee": {  # not available in Hitran
        "altname": "CF3CHFCHFCF2CF3",
        "rfmip": "hfc4310mee_GM",
    },
    # Iodocarbons and hydroiodocarbons
    # Nitriles, amines and other nitrogenated hydrocarbons
    # Other molecules
    "N2O": {"rfmip": "nitrous_oxide_GM",},
    # Sulfur-containing species
}

SPECIES_GROUPS = {
    "reference": ["CCl4", "CF4", "CFC-11", "CFC-12", "HFC-134a", "HFC-23",],
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
