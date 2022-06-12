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
    threads: int = 1
    build_files_path: str = ""


def write_ensemble_file(ens_mems, output_num):
    f = open(f"ires.{output_num}.csv", "w")
    f.write("path,end_0,end_1,end_2\n")
    for ens_mem in ens_mems:
        f.write(f"{ens_mem.path},{ens_mem.end_1},{ens_mem.end_2},\n")
    f.close()


def simulate(df, wrapper, opts):
    scores = []
    for i, row in df.iterrows():
        seq = row["sequence"][8:-8]
        build_file = f"{opts.build_files_path}/{row['topology']}.csv"
        ens_file = os.path.abspath(f"ires.{opts.output_num}.csv")
        avg = 0
        for i in range(opts.runs):
            score = wrapper.run(seq, build_file, ens_file, opts.wrapper_opts)
            if score == -1:
                log.warn('wrapper returned an error!')
                return np.zeros(len(df))
            avg += score

        scores.append(avg / opts.runs + 1)
    return scores


def simulate_single(df, wrapper, opts):
    scores = np.zeros(len(df))
    row = df.iloc[0]
    seq = row["sequence"][8:-8]
    build_file = f"{opts.build_files_path}/{row['topology']}.csv"
    ens_file = os.path.abspath(f"ires.{opts.output_num}.csv")
    avg = 0
    for i in range(opts.runs):
        score = wrapper.run(seq, build_file, ens_file, opts.wrapper_opts)
        if score == -1:
            log.warn('wrapper returned an error!')
            return np.zeros(len(df))
        avg += score
    scores[0] = avg / opts.runs + 1
    return scores


class Simulater(object):
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

    def score(self, ens_members):
        write_ensemble_file(ens_members, self.opts.output_num)
        #init_scores = simulate_single(self.construct_df, self.wrapper, self.opts)
        #if init_scores[0] < 200:
        #    return init_scores
        if self.opts.threads == 1:
            scores = simulate(self.construct_df, self.wrapper, self.opts)
        else:
            wrappers = [self.wrapper for _ in range(self.opts.threads)]
            opts = [self.opts for _ in range(self.opts.threads)]
            dfs = np.array_split(self.construct_df, self.opts.threads)
            score_arrays = self.p.starmap(simulate, zip(dfs, wrappers, opts))
            scores = []
            for score_a in score_arrays:
                scores.extend(score_a)
        return scores

        # print(scores)
