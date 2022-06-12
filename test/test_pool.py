import pandas as pd
import os

from rnamake_ens_gen.optimize import PDBPool

def test_pool():
    path = "/Users/jyesselm/projects/python_packages/rnamake_ens_gen/test"
    pool_path = f"{path}/resources/pool"
    df = pd.read_csv(f"{path}/resources/pool.csv")
    p = PDBPool(df, 1, pool_path)
    inds = p.get_init_ensemble()
    assert len(inds) == 1
    assert inds[0] < len(df)
    # test swap
    new_inds = p.swap_ensemble_member(inds)
    assert new_inds != inds
    # test get pdb
    pdbs = p.get_pdbs(inds)
    assert os.path.isfile(pdbs[0].path)

