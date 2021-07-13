import os
import json
from glob import glob


class HitranMoleculeHeaders:
    def __init__(self, path):
        self.headers = []
        for filename in glob(os.path.join(path, "*.json")):
            with open(filename) as f:
                self.headers.append(json.load(f))
        self.aggregate_names()

    def aggregate_names(self):
        """Store all possible aliases in one list"""
        names = [
            "short_alias", "common_name", "ordinary_formula", "stoichiometric_formula",
            "inchi", "inchikey"
        ]
        for mol in self.headers:
            mol["all_names"] = [mol[name] for name in names]
            mol["all_names"] += [alias["alias"] for alias in mol["aliases"]]
            # Remove duplicates
            mol["all_names"] = list(dict.fromkeys(mol["all_names"]))

    def find(self, name):
        """Search molecule list for match"""
        ret = [mol for mol in self.headers if name in mol["all_names"]]
        return ret

    def is_same(self, name1, name2):
        """Check if aliases belong to the same molecule"""
        match = self.find(name1)
        if len(match) == 0:
            raise RuntimeError(f"{name1} not found")
        if len(match) > 1:
            raise RuntimeError(f"Found more than 1 match: {match}")
        return name2 in match[0]["all_names"]


if __name__ == "__main__":
    script_path = os.path.dirname(os.path.realpath(__file__))
    INPUTDIR = os.path.join(script_path, "../data/HitranMoleculeHeaders")

    hmh = HitranMoleculeHeaders(INPUTDIR)
    assert hmh.is_same("HFC-125", "f125")

    molecules = [mol["short_alias"] for mol in hmh.headers]
    molecules.sort()
    print("\n".join(molecules))
