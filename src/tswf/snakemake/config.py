import contextlib
import os

import numpy as np
import pandas as pd
from snakemake.utils import logger
from tswf._version import version
from tswf.config import get_schema
from tswf.config import PKG_DIR
from tswf.config import PropertyDict


SNAKEMAKE_ROOT = PKG_DIR / "workflow"


def wildcards_or(items, empty=False):
    items = list(set(items))
    if empty:
        items = [""] + items
    return f'({"|".join(items)})'


# context manager for cd
@contextlib.contextmanager
def cd(path, logger):
    CWD = os.getcwd()
    logger.info(f"Changing directory from {CWD} to {path}")

    os.chdir(path)
    try:
        yield
    except Exception as e:
        logger.warning(f"Exception caught: {e}")
    finally:
        logger.info(f"Changing directory back to {CWD}")
        os.chdir(CWD)


def add_gitinfo(config, workflow):
    config["__workflow_basedir__"] = workflow.basedir
    config["__workflow_workdir__"] = os.getcwd()
    config["__version__"] = version
    config["__workflow_commit_link__"] = "https://github.com/percyfal/wg-genealogy-smk"


class Data:
    _index = None
    _schema = None

    def __init__(self, *args):
        if len(args) == 1:
            args = args[0]
            if isinstance(args, Data):
                self._schemafile = args._schemafile
                self._data = args.data
            elif isinstance(args, str):
                self._read_tsv(args)
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
        self.schema.validate(self.data.reset_index())

    def _read_tsv(self, infile):
        sep = ","
        if infile.endswith(".tsv"):
            sep = "\t"
        self._data = pd.read_csv(infile, sep=sep).set_index(self._index, drop=False)
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
