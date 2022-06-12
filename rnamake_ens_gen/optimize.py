import logging
import random
import numpy as np
import pandas as pd
from typing import List
from dataclasses import dataclass
from scipy.stats import pearsonr

from rnamake_ens_gen import simulate, wrapper, logger
from rnamake_ens_gen.anneal import Annealer, AnnealerParams


def normalize_data(data):
    if np.min(data) == np.max(data):
        return data
    return (data - np.min(data)) / (np.max(data) - np.min(data))


class ScoreType:
    LINEAR_LINEAR = 0
    LOG_LINEAR = 1
    LOG_LOG = 2


@dataclass(frozen=True, order=True)
class EnsembleMember:
    path: str
    end_1: str
    end_2: str
    end_3: str = ""


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

    def get_pdbs(self, ens) -> List[EnsembleMember]:
        ens_members = []
        for i in ens:
            row = self.df.iloc[i]
            path = self.path + "/" + str(row["pdb"]) + ".pdb"
            ens_members.append(EnsembleMember(path, row["end1"], row["end2"]))
        return ens_members


class Scorer(object):
    def __init__(self, constructs: pd.DataFrame, output_num):
        self.constructs = constructs
        self.output_num = output_num
        self.df_results = pd.DataFrame(columns="members,data,score".split(","))
        self.pos = 0

    def __record_score(self, members, scores, score):
        self.df_results.loc[self.pos] = [
            [x.path for x in members],
            scores,
            score,
        ]
        self.pos += 1
        self.df_results.to_json(f"log.{self.output_num}.json", orient="records")

    def score(self, scores, ens_mems):
        sum = np.sum(scores)
        if sum == scores[0]:
            score = len(ens_mems) * 2 - np.log(sum)
        else:
            df = self.constructs
            x = normalize_data(df["exp_score"])
            y = normalize_data(-np.log(scores))
            score = 0
            for x1, y1 in zip(x, y):
                score += abs(x1 - y1)
        # print(self.i, score, scores, self.state)
        self.__record_score(ens_mems, scores, score)
        return score


class SampleandSelect(Annealer):
    def __init__(self, simulator: simulate.Simulater, pool: PDBPool):
        self.simulator = simulator
        self.pool = pool
        self.hash = {}
        self.scorer = Scorer(simulator.construct_df, simulator.opts.output_num)
        super(SampleandSelect, self).__init__(self.pool.get_init_ensemble())

    def move(self):
        initial_energy = self.energy()
        self.state = self.pool.swap_ensemble_member(self.state)
        return self.energy() - initial_energy

    def update(self, *args, **kwargs):
        pass

    def energy(self):
        key = ""
        for i in self.state:
            key += str(i) + "_"
        if key in self.hash:
            return self.hash[key]
        ens_mems = self.pool.get_pdbs(self.state)
        scores = self.simulator.score(ens_mems)
        score = self.scorer.score(scores, ens_mems)
        self.hash[key] = score
        return score


def main():
    log = logger.setup_applevel_logger()
    log.setLevel(logging.DEBUG)
    path = "/Users/jyesselm/projects/python_packages/rnamake_ens_gen/test"
    w_opts = wrapper.Opts(extra_pdbs=f"{path}/resources/GAAA_tetraloop.pdb")
    opts = simulate.Opts(
        wrapper_opts=w_opts,
        build_files_path=path + "/resources/build_files",
        threads=2,
    )
    df_constructs = pd.read_csv(f"{path}/resources/test.csv")
    s = simulate.Simulater(df_constructs, opts)
    s.setup()
    pool_path = f"{path}/resources/pool"
    df = pd.read_csv(f"{path}/resources/pool.csv")
    pool = PDBPool(df, 5, pool_path)
    tsp = SampleandSelect(s, pool)
    tsp.steps = 1
    tsp.Tmax = 10
    tsp.Tmin = 0.011
    tsp.updates = 10
    params = AnnealerParams(1, 10, 0.011)
    state, e = tsp.anneal(params)


if __name__ == "__main__":
    main()
