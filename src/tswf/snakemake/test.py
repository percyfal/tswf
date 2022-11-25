import os
import shutil
import tempfile

ROOT = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
TESTDIR = os.path.join(ROOT, "tests")


def copytree_testdir(dname):
    """Copy test directory"""
    if dname is None:
        dname = tempfile.TemporaryDirectory().name
    try:
        outdir = shutil.copytree(TESTDIR, dname, dirs_exist_ok=True)
    except Exception as e:
        print(e)
        raise
    return outdir
