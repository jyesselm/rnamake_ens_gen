import pandas as pd
from rnamake_ens_gen import score, wrapper


def _test_scorer():
    path = "/Users/jyesselm/projects/python_packages/rnamake_ens_gen/test"
    w_opts = wrapper.Opts(extra_pdbs=f"{path}/resources/GAAA_tetraloop.pdb")
    opts = score.Opts(
        wrapper_opts=w_opts, build_files_path=path + "/resources/build_files"
    )
    df_constructs = pd.read_csv(f"{path}/resources/double.csv")
    s = score.Scorer(df_constructs, opts)
    s.setup()
    pdbs = [path + "/resources/pool/1.pdb"]
    s_val = s.score(pdbs)

def test_threads_scorer():
    path = "/Users/jyesselm/projects/python_packages/rnamake_ens_gen/test"
    w_opts = wrapper.Opts(extra_pdbs=f"{path}/resources/GAAA_tetraloop.pdb")
    opts = score.Opts(
            wrapper_opts=w_opts, build_files_path=path + "/resources/build_files",
            threads=3
    )
    df_constructs = pd.read_csv(f"{path}/resources/subset.csv")
    s = score.Scorer(df_constructs, opts)
    s.setup()
    pdbs = [path + "/resources/pool/1.pdb"]
    s_val = s.score(pdbs)
