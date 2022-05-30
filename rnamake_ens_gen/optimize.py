from simanneal import Annealer
import random
import numpy as np
import pandas as pd

from rnamake_ens_gen import score, wrapper


class ScoreType:
    LINEAR_LINEAR = 0
    LOG_LINEAR = 1
    LOG_LOG = 2


class PDBPool(object):
    def __init__(self, df: pd.DataFrame, size, path):
        self.df = df
        self.size = size
        self.path = path

    def get_init_ensemble(self):
        return [random.randint(0, len(self.df) - 1) for _ in range(self.size)]

    def swap_ensemble_member(self, ens):
        while 1:
            pos = random.randint(0, len(ens) - 1)
            new_ind = random.randint(0, len(self.df) - 1)
            if new_ind not in ens:
                ens[pos] = new_ind
                break
        return ens

    def get_pdbs(self, ens):
        pdbs = []
        for i in ens:
            pdbs.append(self.path + "/" + str(i) + ".pdb")
        return pdbs




class SampleandSelect(Annealer):
    def __init__(self, scorer: score.Scorer, pool: PDBPool):
        self.scorer = scorer
        self.pool = pool
        super(SampleandSelect, self).__init__(self.pool.get_init_ensemble())

    def move(self):
        pass

    def energy(self):
        pass


def main():
    path = "/Users/jyesselm/projects/python_packages/rnamake_ens_gen/test"
    w_opts = wrapper.Opts(extra_pdbs=f"{path}/resources/GAAA_tetraloop.pdb")
    opts = score.Opts(
        wrapper_opts=w_opts,
        build_files_path=path + "/resources/build_files",
        threads=3,
    )
    df_constructs = pd.read_csv(f"{path}/resources/subset.csv")
    s = score.Scorer(df_constructs, opts)
    s.setup()

    pass


if __name__ == "__main__":
    main()
