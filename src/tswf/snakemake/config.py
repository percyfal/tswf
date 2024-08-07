import contextlib
import os

import numpy as np
import pandas as pd


try:
    from snakemake.logging import logger
except ImportError:
    from snakemake.utils import logger

from tswf._version import version
from tswf.config import PKG_DIR
from tswf.config import PropertyDict
from tswf.config import get_schema


SNAKEMAKE_ROOT = PKG_DIR / "workflow"


def wildcards_or(items: list[str], empty: bool = False):
    """Convert list of items to wildcards string separated by '|'."""
    items = list(set(items))
    if empty:
        items = [""] + items
    return f'({"|".join(items)})'


@contextlib.contextmanager
def cd(path, logger):
    """Change directory and return to original on exit (context manager for cd)."""
    cwd = os.getcwd()
    logger.info("Changing directory from %s to %s", cwd, path)

    os.chdir(path)
    try:
        yield
    except Exception as e:
        logger.warning("Exception caught: %s", e)
    finally:
        logger.info("Changing directory back to %s", cwd)
        os.chdir(cwd)


def add_gitinfo(config, workflow):
    """Add gitinfo to config dict."""
    config["__workflow_basedir__"] = workflow.basedir
    config["__workflow_workdir__"] = os.getcwd()
    config["__version__"] = version
    config["__workflow_commit_link__"] = "https://github.com/percyfal/wg-genealogy-smk"


class Data:
    """Generic data structure with schema and index."""

    _index: None | list[str] = None
    _schema = None

    def __init__(self, *args):
        if len(args) == 1:
            args = args[0]
            if isinstance(args, Data):
                self._schemafile = args._schemafile
                self._data = args.data
            elif isinstance(args, str):
                self._read_tsv(args)
            elif isinstance(args, pd.DataFrame):
                self._data = args
        elif len(args) > 1:
            assert all(isinstance(x, Data) for x in args), logger.error(
                "all instances must be Data"
            )
            self._schemafile = args[0]._schemafile
            self._data = pd.concat(x.data for x in args)
        else:
            raise TypeError
        try:
            self._data.set_index(self._index, drop=True, inplace=True)
        except Exception:
            pass
        if isinstance(self.data, pd.DataFrame):
            for _, record in enumerate(self.data.reset_index().to_dict("records")):
                self.schema.validate(record)
        else:
            self.schema.validate(self.data.reset_index())

    def _read_tsv(self, infile):
        sep = ","
        if infile.endswith(".tsv"):
            sep = "\t"
        self._data = pd.read_csv(infile, sep=sep, comment="#").set_index(
            self._index, drop=False
        )
        self._data = self._data.replace({np.nan: None})
        self._data.index.names = self._index

    @property
    def schema(self):
        return self._schema

    @property
    def data(self):
        return self._data

    @property
    def populations(self):
        return self.data.population

    def subset(self, invert=False, **kw):
        keep = self._data.transpose().all()
        for k, v in kw.items():
            if isinstance(v, set):
                v = list(v)
            if not isinstance(v, list):
                v = [v]
            if len(v) == 0:
                continue
            try:
                if k in self.data.index.names:
                    keep = keep & self.data.index.isin(v)
                else:
                    keep = keep & self.data[k].isin(v)
            except KeyError as e:
                print(e)
                raise
        if invert:
            keep = ~keep
        cls = type(self)
        new = cls(self)
        new._data = new._data[keep]
        return new

    def merge(self, other, **kw):
        self._data = pd.merge(self.data, other.data, **kw)


class SampleData(Data):
    """Data structure to hold samples."""

    _index = ["SM"]
    _schema = get_schema("SAMPLES_SCHEMA")

    def __init__(self, *args):
        super().__init__(*args)

    @property
    def samples(self):
        return self.data.index

    @property
    def unique_samples(self):
        return self.data.SM[list(set(self.data.SM.tolist()))]

    def subset(self, invert=False, **kw):
        if "samples" in kw.keys():
            kw["SM"] = kw["samples"]
            kw.pop("samples")
        return super().subset(invert, **kw)


class PopulationData(Data):
    _index = ["population"]
    _schema = get_schema("POPULATIONS_SCHEMA")

    def __init__(self, *args):
        super().__init__(*args)

    @property
    def populations(self):
        return self.data.index


class Config(PropertyDict):
    def __init__(self, conf, samples):
        super().__init__(conf)
        self["samples"] = samples

    @property
    def genome(self):
        """Return short genome name"""
        return os.path.splitext(os.path.basename(self.ref))[0]
