#!/usr/bin/env python3
import os
import json
import logging
from gzip import GzipFile

from hitran_xsec import (
    XsecFileIndex,
    XsecFile,
    set_default_logging_format,
)
from hitran_molecule_headers import HitranMoleculeHeaders
from xsec_species_info import SPECIES_GROUPS, XSEC_SPECIES_INFO

set_default_logging_format(level=logging.INFO,
                           include_timestamp=True,
                           include_function=True)
logger = logging.getLogger(__name__)

script_path = os.path.dirname(os.path.realpath(__file__))

INPUTDIR = os.path.join(script_path, "../data/HitranXsec")
OUTPUTDIR = os.path.join(script_path, "../data/HitranXsecJson")
species_list = SPECIES_GROUPS["rfmip"]


def xsc_to_json(species):
    MOLECULEDIR = os.path.join(script_path, "../data/HitranMoleculeHeaders")
    hmh = HitranMoleculeHeaders(MOLECULEDIR)

    outfile = os.path.join(OUTPUTDIR, f"{species}.xsc.json.gz")
    xfi = XsecFileIndex(INPUTDIR, species, molecule_headers=hmh)
    bands = xfi.cluster_by_band()

    bands2 = []
    for b in bands:
        band2 = []
        x: XsecFile
        for x in b:
            band2.append({
                "species": x.species,
                "xscfile": x.filename,
                "wmin": x.wmin,
                "wmax": x.wmax,
                "fmin": x.fmin,
                "fmax": x.fmax,
                "pressure": x.pressure,
                "temperature": x.temperature,
                "xsec": list(x.data),
            })
        bands2.append(band2)

    if len(bands2) == 0:
        logger.error(f"No data found for species {species}")
    else:
        logger.info(f"Writing json file {outfile}")
        with GzipFile(outfile, "w") as f:
            f.write(json.dumps(bands2).encode("utf-8"))


for s in species_list:
    xsc_to_json(s)

print()
print("Unavailable RFMIP species:")
for rfmip_species in SPECIES_GROUPS["rfmip-names"]:
    found = False
    for v in XSEC_SPECIES_INFO.values():
        if v["rfmip"] == rfmip_species:
            found = True
            break
    if not found:
        print(rfmip_species)

# Example for reading the data
# import json
# from gzip import GzipFile
#
# with GzipFile(f"CF4-xsc.json.gz") as f:
#     cf4 = json.loads(f.read().decode("utf-8"))
