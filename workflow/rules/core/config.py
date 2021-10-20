import os
import pandas as pd
import numpy as np
import contextlib
import subprocess as sp
from collections import OrderedDict
from snakemake.utils import logger, validate
from snakemake.io import _load_configfile


WORKFLOW_DIR = str(workflow.current_basedir)
SCHEMA_DIR = os.path.realpath(
    os.path.join(WORKFLOW_DIR, os.pardir, os.pardir, "schemas")
)


def wildcards_or(items, empty=False):
    items = list(set(items))
    if empty:
        items = [""] + items
    return f'({"|".join(items)})'


def x_runtime(wildcards, attempt):
    return attempt * workflow.default_resources.parsed.get("runtime", 100)


def x_mem_mb(wildcards, attempt):
    default = workflow.default_resources.defaults["mem_mb"]
    return attempt * workflow.default_resources.parsed.get("mem_mb", default)


def xx_mem_mb(wildcards, attempt):
    default = workflow.default_resources.defaults["mem_mb"]
    return (2 ** attempt) * workflow.default_resources.parsed.get("mem_mb", default)


# context manager for cd
@contextlib.contextmanager
def cd(path, logger):
    CWD = os.getcwd()
    logger.info("Changing directory from {} to {}".format(CWD, path))

    os.chdir(path)
    try:
        yield
    except:
        logger.warning("Exception caught: ".format(sys.exc_info()[0]))
    finally:
        logger.info("Changing directory back to {}".format(CWD))
        os.chdir(CWD)


def add_gitinfo(config):
    config["__workflow_basedir__"] = workflow.basedir
    config["__workflow_workdir__"] = os.getcwd()
    config["__worfklow_commit__"] = None
    config["__workflow_commit_link__"] = None

    try:
        with cd(workflow.basedir, logger):
            commit = sp.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
            commit_short = (
                sp.check_output(["git", "rev-parse", "--short", "HEAD"])
                .decode()
                .strip()
            )
            rc = sp.run(["git", "diff", "--quiet"])
            dirty = "-dirty" if rc.returncode == 1 else ""
            config["__workflow_commit__"] = commit_short + dirty
            config[
                "__workflow_commit_link__"
            ] = f"https://github.com/NBISweden/manticore-smk/commit/{commit}"
    except Exception as e:
        print(e)
        raise


class PropertyDict(OrderedDict):
    """Simple class that allows for property access"""

    def __init__(self, data=dict()):
        super().__init__(data)
        for k, v in data.items():
            if isinstance(v, dict):
                v = PropertyDict(v)
            elif isinstance(v, list):
                val = []
                for x in v:
                    if isinstance(x, PropertyDict):
                        val.append(x)
                    elif isinstance(x, dict):
                        val.append(PropertyDict(x))
                    else:
                        val.append(x)
                v = val
            else:
                pass
            self[k] = v

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        if key not in dir(dict()):
            try:
                setattr(self, key, value)
            except Exception as e:
                print(e)
                print(key, value)
                raise


class Schema(PropertyDict):
    def __init__(self, schemafile):
        self._schemafile = schemafile
        data = _load_configfile(self.schemafile, filetype="Schema")
        super().__init__(data)

    @property
    def schemafile(self):
        return self._schemafile


sample_schema = Schema(os.path.join(SCHEMA_DIR, "samples.schema.yaml"))
population_schema = Schema(os.path.join(SCHEMA_DIR, "populations.schema.yaml"))


##############################
## Data related classes
##############################
class Data:
    _index = None
    _schemafile = None

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
        except:
            pass
        validate(self.data.reset_index(), schema=self.schemafile)

    def _read_tsv(self, infile):
        sep = ","
        if infile.endswith(".tsv"):
            sep = "\t"
        self._data = pd.read_csv(infile, sep=sep).set_index(self._index, drop=False)
        self._data = self._data.replace({np.nan: None})
        self._data.index.names = self._index

    @property
    def schemafile(self):
        return self._schemafile

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
    _schemafile = sample_schema.schemafile

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
    _schemafile = population_schema.schemafile

    def __init__(self, *args):
        super().__init__(*args)

    @property
    def populations(self):
        return self.data.index


##############################
# Config related classes
##############################
class Analysis(PropertyDict):
    def __init__(self, name, default, **kw):
        default.update(**kw)
        super().__init__(default)
        if "subset" in self.keys():
            if "include" in self["subset"]:
                self["samples"] = self.samples.subset(**self.subset.include)
            if "exclude" in self["subset"]:
                self["samples"] = self.samples.subset(
                    **self.subset.exclude, invert=True
                )
        self._name = name

    @property
    def section(self):
        return self._name.split("/")[0]

    @property
    def shortname(self):
        return self._name.lstrip(f"{self.section}/")

    @property
    def name(self):
        return self._name


class ConfigRule(PropertyDict):
    def __init__(self, name, attempt=None, *args, **kw):
        self._name = name
        self._attempt = attempt
        super().__init__(*args, **kw)

    @property
    def attempt(self):
        if self._attempt is None:
            return 1
        return self._attempt

    @property
    def name(self):
        return self._name

    def resources(self, resource):
        assert isinstance(
            self[resource], int
        ), f"{self}: resource '{resource}' is not an int"
        return self.attempt * self[resource]

    def params(self, attr):
        return self[attr]


class Config(PropertyDict):
    def __init__(self, conf, samples):
        super().__init__(conf)
        self["samples"] = samples
        self._init_analyses()

    def _init_analyses(self):
        for k in self.keys():
            if not re.match(r"^(tsinfer|relate)/", k):
                continue
            default = {
                "chromosomes": self.chromosomes,
                "samples": self.samples,
            }
            self[k] = Analysis(k, default, **self[k])

    def ruleconf(self, rulename, attempt=None, analysis=None):
        """Retrieve rule configuration"""
        ruleobj = ConfigRule(rulename, attempt, self["resources.default"])
        if rulename in self.rules:
            ruleobj.update(**self.rules[rulename])
        if analysis is not None:
            if rulename in self.get_analysis(analysis).rules:
                ruleobj.update(**self.get_analysis(analysis).rules[rulename])
        return ruleobj

    def get_analysis(self, name):
        return self[name]

    def analyses(self, label=None):
        if label is not None:
            return [self[k] for k in self.keys() if k.startswith(label)]
        return [self[k] for k in self.keys() if re.match(r"^(tsinfer|relate)/", k)]

    @property
    def analysisnames(self):
        return [x.name for x in self.analyses()]

    @property
    def genome(self):
        """Return short genome name"""
        return os.path.splitext(os.path.basename(self.ref))[0]

    def params(self, wildcards, attr, rulename=None):
        """Get parameters in analysis context"""
        analysis = self.get_analysis(wildcards.analysis)
        default = None
        if rulename is not None:
            default = self.ruleconf(rulename).params(attr)
        return analysis_item.get(attr, default)
