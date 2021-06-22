import itertools
import logging
import os
import re
from copy import deepcopy
from glob import glob
from typing import Iterable

import numpy as np
from scipy.constants import speed_of_light

from xsec_species_info import XSEC_SPECIES_INFO

logger = logging.getLogger(__name__)


def set_default_logging_format(level=None,
                               include_timestamp=True,
                               include_function=True):
    """Generate decently looking logging format string."""

    if level is None:
        level = logging.INFO

    color = "\033[1;%dm"
    reset = "\033[0m"
    black, red, green, yellow, blue, magenta, cyan, white = [
        color % (30 + i) for i in range(8)
    ]
    logformat = "["
    if include_timestamp:
        logformat += f"{red}%(asctime)s.%(msecs)03d{reset}:"
    logformat += f"{yellow}%(filename)s{reset}" f":{blue}%(lineno)s{reset}"
    if include_function:
        logformat += f":{green}%(funcName)s{reset}"
    logformat += f"] %(message)s"

    logging.basicConfig(format=logformat, level=level, datefmt="%H:%M:%S")


def wavenumber2frequency(wavenumber):
    """Convert wavenumber to frequency.

    Parameters:
        wavenumber (float or ndarray): Wavenumber [m^-1].
    Returns:
        float or ndarray: Frequency [Hz].

    """
    return speed_of_light * wavenumber


class XsecError(RuntimeError):
    """Cross section related RuntimeError."""

    pass


class XsecFile:
    """HITRAN cross section file."""
    def __init__(self, filename):
        """Lazy-load cross section file."""
        self.filename = filename
        # noinspection PyUnusedLocal
        rnum = r"[0-9]+\.?[0-9]*"
        m = re.search(
            f"(?P<species>[^_]*)_(?P<T>{rnum})K?[-_](?P<P>{rnum})(Torr|K)?[-_]"
            f"(?P<wmin>{rnum})[-_](?P<wmax>{rnum})(?P<extra>_.*)?\.xsc",
            os.path.basename(self.filename),
        )
        try:
            self.species = m.group("species")
            self.temperature = float(m.group("T"))
            self.torr = float(m.group("P"))
            self.pressure = torr_to_pascal(self.torr)
            self.wmin = float(m.group("wmin"))
            self.wmax = float(m.group("wmax"))
            self.fmin = wavenumber2frequency(self.wmin * 100)
            self.fmax = wavenumber2frequency(self.wmax * 100)
            self.extra = m.group("extra")
            self._header = None
            self._data = None
            self._nfreq = None
        except AttributeError:
            raise XsecError(f"Error parsing filename {filename}")

    def __repr__(self):
        return "XsecFile:" + self.filename

    def __hash__(self):
        return hash(f"{self.species}{self.pressure}{self.temperature}"
                    f"{self.wmin}{self.wmax}")

    def __eq__(self, x):
        return (self.species == x.species and self.pressure == x.pressure
                and self.temperature == x.temperature and self.wmin == x.wmin
                and self.wmax == x.wmax)

    def to_dict(self):
        return {
            "species": self.species,
            "wmin": self.wmin,
            "wmax": self.wmax,
            "fmin": self.fmin,
            "fmax": self.fmax,
            "pressure": self.pressure,
            "temperature": self.temperature,
            "nfreq": self.nfreq,
            "data": deepcopy(self.data),
        }

    def read_hitran_xsec(self):
        """Read HITRAN cross section data file."""
        if self._data is not None:
            return

        logger.info(f"Reading {self.filename}")
        with open(self.filename) as f:
            header = f.readline()
            data = np.hstack(
                list(map(lambda l: list(map(float, l.split())),
                         f.readlines())))

        self._header = header
        self._data = data
        self._nfreq = len(data)

    @property
    def nfreq(self):
        if self._nfreq is None:
            self.read_hitran_xsec()
        return self._nfreq

    @property
    def header(self):
        if self._header is None:
            self.read_hitran_xsec()
        return self._header

    @property
    def data(self):
        if self._data is None:
            self.read_hitran_xsec()
        return self._data

    @data.setter
    def data(self, val):
        self._data = val


class XsecFileIndex:
    """Database of HITRAN cross section files."""
    def __init__(self, directory=None, species=None, ignore=None):
        self.files = []
        self.ignored_files = []
        self.failed_files = []
        if directory is not None and species is not None:
            if "altname" in XSEC_SPECIES_INFO[species]:
                speciesname = XSEC_SPECIES_INFO[species]["altname"]
            else:
                speciesname = species

            for f in glob(os.path.join(directory, "*.xsc")):
                try:
                    xsec_file = XsecFile(f)
                    if xsec_file.species != speciesname:
                        pass
                    elif ignore is not None and re.match(
                            ignore, xsec_file.extra):
                        self.ignored_files.append(f)
                    else:
                        self.files.append(xsec_file)
                        if species != speciesname:
                            xsec_file.species = species
                except XsecError:
                    self.failed_files.append(f)
        self.uniquify()

    @classmethod
    def from_list(cls, xsec_file_list):
        obj = cls()
        obj.files = xsec_file_list
        return obj

    def __repr__(self):
        return "\n".join([f.filename for f in self.files])

    def uniquify(self):
        nfiles = len(self.files)
        checked = {}
        uniqfiles = []
        for item in self.files:
            marker = item
            if marker in checked:
                continue
            checked[marker] = 1
            uniqfiles.append(item)
        nuniqfiles = len(uniqfiles)
        if nuniqfiles < nfiles:
            logger.info(f"Removed {nfiles - nuniqfiles} duplicate data files "
                        f"for {self.files[0].species}")
            self.files = uniqfiles

    def find_file(self, filename):
        ret = [x for x in self.files if x.filename == filename]
        return ret if len(ret) > 1 else ret[0]

    def find(self, wmin=None, wmax=None, temperature=None, pressure=None):
        """Find cross sections that match the criteria."""
        return [
            x for x in self.files
            if (not wmin or x.wmin == wmin) and (not wmax or x.wmax == wmax)
            and (not temperature or x.temperature == temperature) and (
                not pressure or x.torr == pressure)
        ]

    def cluster_by_band(self, wgap=1):
        """Combine files for each band in a list."""
        return _cluster2(self.files,
                         wgap,
                         key=lambda x: x.wmin,
                         key2=lambda x: x.wmax)

    def cluster_by_temperature(self, tgap=3):
        """Combine files for each temperature in a list."""
        return _cluster2(self.files, tgap, key=lambda x: x.temperature)

    def cluster_by_band_and_pressure(self, wgap=1, pgap=100):
        """Combine files for each band and pressure in a nested list."""
        return (_cluster2(
            l, pgap, key=lambda x: x.pressure) for l in _cluster2(
                self.files, wgap, key=lambda x: x.wmin, key2=lambda x: x.wmax))

    def cluster_by_band_and_temperature(self, wgap=1, tgap=3):
        """Combine files for each band and temperature in a nested list."""
        return (_cluster2(
            l, tgap, key=lambda x: x.temperature) for l in _cluster2(
                self.files, wgap, key=lambda x: x.wmin, key2=lambda x: x.wmax))


def torr_to_pascal(torr):
    """Convert Torr to Pascal."""
    return torr * 101325.0 / 760.0


def pascal_to_torr(pascal):
    """Convert Pascal to Torr."""
    return pascal / 101325.0 * 760.0


def _pairify(it):
    """Build pairs."""
    it0, it1 = itertools.tee(it, 2)
    first = next(it0)
    return zip(itertools.chain([first, first], it0), it1)


def _cluster2(iterable: Iterable, maxgap, key=lambda x: x, key2=None):
    """Cluster sequence elements by distance."""
    prev = None
    group = []
    for item in sorted(iterable,
                       key=lambda x: (key(x), key2(x))
                       if key2 is not None else key(x)):
        if not prev or (key(item) - key(prev) <= maxgap and
                        (not key2 or key2(item) - key2(prev) <= maxgap)):
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group
