import subprocess
import shutil
import numpy as np
import os
from dataclasses import dataclass
from multiprocessing import Pool

from rnamake_ens_gen import logger, wrapper

log = logger.get_logger("score")



@dataclass(frozen=True, order=True)
class Opts:
    wrapper_opts: wrapper.Opts = wrapper.Opts()
    output_num: int = 0
    runs: int = 1
    threads : int = 1
    build_files_path : str = ""


def write_ensemble_file(pdbs, output_num):
    f = open(f"ires.{output_num}.csv", "w")
    f.write("path,end_0,end_1,end_2\n")
    for pdb in pdbs:
        f.write(f"{pdb},A2-B14,A8-B13,\n")
    f.close()


def score(df, wrapper, opts):
    scores = []
    for i, row in df.iterrows():
        seq = row['sequence'][8:-8]
        build_file = f"{opts.build_files_path}/{row['topology']}.csv"
        ens_file = os.path.abspath("ires.0.csv")
        score = wrapper.run(seq, build_file, ens_file, opts.wrapper_opts)
        scores.append(score)
    return scores



class Scorer(object):
    def __init__(self, construct_df, opts: Opts):
        self.opts = opts
        self.construct_df = construct_df
        for col in ["sequence", "topology", "exp_score"]:
            if col not in self.construct_df:
                log.error(f"{col} must be included in construct dataframe")
                exit()

    def setup(self):
        self.wrapper = wrapper.BuildMotifGraphWrapper()
        self.wrapper.setup()
        self.p = Pool(processes=self.opts.threads)

    def score(self, pdbs):
        write_ensemble_file(pdbs, self.opts.output_num)
        if self.opts.threads == 1:
            scores = score(self.construct_df, self.wrapper, self.opts)
        else:
            wrappers = [self.wrapper for _ in range(self.opts.threads)]
            opts = [self.opts for _ in range(self.opts.threads)]
            dfs = np.array_split(self.construct_df, self.opts.threads)
            score_arrays = self.p.starmap(score, zip(dfs, wrappers, opts))
            scores = []
            for score_a in score_arrays:
                scores.extend(score_a)
        return scores

        #print(scores)
