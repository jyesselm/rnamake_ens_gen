from simanneal import Annealer
import random
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from rnamake_ens_gen import score, wrapper

def normalize_data(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

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
        new_ens = ens[:]
        while 1:
            pos = random.randint(0, len(new_ens) - 1)
            new_ind = random.randint(0, len(self.df) - 1)
            if new_ind not in new_ens:
                new_ens[pos] = new_ind
                break
        return new_ens

    def get_pdbs(self, ens):
        pdbs = []
        for i in ens:
            pos = self.df.iloc[i]["pdb"]
            pdbs.append(self.path + "/" + str(int(pos)) + ".pdb")
        return pdbs


class SampleandSelect(Annealer):
    def __init__(self, scorer: score.Scorer, pool: PDBPool):
        self.scorer = scorer
        self.pool = pool
        self.i = 0
        self.pos = 0
        self.hash = {}
        self.df_results = pd.DataFrame(
            columns="members,data,score".split(",")
        )
        super(SampleandSelect, self).__init__(self.pool.get_init_ensemble())

    def update(self, step, T, E, acceptance, improvement):
        if step == 0:
            return
        print("updates:", T, E, 100.0 * acceptance, 100.0 * improvement)

    def move(self):
        self.i += 1
        initial_energy = self.energy()
        self.state = self.pool.swap_ensemble_member(self.state)
        return self.energy() - initial_energy

    def __record_score(self, members, scores, score):
        self.df_results.loc[self.pos] = [members, scores, score]
        self.pos += 1
        output_num = self.scorer.opts.output_num
        self.df_results.to_json(f"log.{output_num}.json", orient="records")


    def energy(self):
        key = ""
        for i in self.state:
            key += str(i) + "_"
        if key in self.hash:
            return self.hash[key]
        pdbs = self.pool.get_pdbs(self.state)
        scores = self.scorer.score(pdbs)
        sum = np.sum(scores)
        if sum == scores[0]:
            score = len(self.state)*2 - np.log(sum)
        else:
            df = self.scorer.construct_df
            x = normalize_data(np.log(df["exp_score"]))
            y = normalize_data(-np.log(scores))
            score = 0
            for x1, y1 in zip(x, y):
                score += abs(x1 - y1)
        print(self.i, score, scores, self.state)
        self.__record_score(pdbs, scores, score)
        self.hash[key] = score
        return score


def main():
    path = "/Users/jyesselm/projects/python_packages/rnamake_ens_gen/test"
    w_opts = wrapper.Opts(extra_pdbs=f"{path}/resources/GAAA_tetraloop.pdb")
    opts = score.Opts(
        wrapper_opts=w_opts,
        build_files_path=path + "/resources/build_files",
        threads=2,
    )
    df_constructs = pd.read_csv(f"{path}/resources/test.csv")
    s = score.Scorer(df_constructs, opts)
    s.setup()
    pool_path = f"{path}/resources/pool"
    df = pd.read_csv(f"{path}/resources/pool.csv")
    pool = PDBPool(df, 5, pool_path)
    tsp = SampleandSelect(s, pool)
    tsp.steps = 10
    tsp.Tmax = 0.1
    tsp.Tmin = 0.001
    state, e = tsp.anneal()


if __name__ == "__main__":
    main()
