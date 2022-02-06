try:
    import importlib.metadata as meta  # Python >=3.8
except ModuleNotFoundError:
    import importlib_metadata as meta  # Python 3.7

import pathlib

PACKAGE = "fretlabel"

MODULE_DIR = pathlib.Path(__file__).parent

try:
    metadata = meta.metadata(PACKAGE)
except meta.PackageNotFoundError:
    print(f"{PACKAGE} is not installed yet.")


def _get_urls():
    __urls__ = {}
    try:
        for _url in metadata.get_all("Project-URL"):
            __urls__.update(dict([_url.replace(" ", "").split(",")]))
    except TypeError:
        __urls__["Documentation"] = None
        __urls__["Repository"] = None
    return __urls__


__version__ = metadata["Version"]
__author__ = metadata["Author"]
__keywords__ = metadata["Keywords"]
__license__ = metadata["License"]
__urls__ = _get_urls()


from fretlabel import ff
