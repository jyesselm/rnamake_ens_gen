import pandas as pd
from rnamake_ens_gen import simulate, wrapper, optimize


def test_simulate():
    path = "/Users/jyesselm/projects/python_packages/rnamake_ens_gen/test"
    w_opts = wrapper.Opts(extra_pdbs=f"{path}/resources/GAAA_tetraloop.pdb")
    opts = simulate.Opts(
        wrapper_opts=w_opts, build_files_path=path + "/resources/build_files"
    )
    df_constructs = pd.read_csv(f"{path}/resources/double.csv")
    s = simulate.Simulater(df_constructs, opts)
    s.setup()
    pdbs = [
        optimize.EnsembleMember(
            path + "/resources/pool/1.pdb", "A2-B14", "A8-B13"
        )
    ]
    s_val = s.score(pdbs)
    assert s_val[0] > 0


def _test_threads_simulate():
    path = "/Users/jyesselm/projects/python_packages/rnamake_ens_gen/test"
    w_opts = wrapper.Opts(extra_pdbs=f"{path}/resources/GAAA_tetraloop.pdb")
    opts = simulate.Opts(
        wrapper_opts=w_opts,
        build_files_path=path + "/resources/build_files",
        threads=2,
    )
    df_constructs = pd.read_csv(f"{path}/resources/double.csv")
    s = simulate.Simulater(df_constructs, opts)
    s.setup()
    pdbs = [path + "/resources/pool/1.pdb"]
    s_val = s.score(pdbs)
    print(s_val)
