import click
import pandas as pd
import logging
import itertools

from rnamake_ens_gen import logger, simulate, optimize, wrapper
from rnamake_ens_gen.anneal import AnnealerParams


@click.group()
def cli():
    pass


@cli.command()
@click.argument("constructs")
@click.argument("pool")
@click.option("-bf", "--build-files", required=True)
@click.option("-pdbs", "--pdbs", required=True)
@click.option("-pp", "--pool-path", required=True)
@click.option("-t", "--threads", default=1, type=int)
@click.option("-es", "--ensemble-size", default=2, type=int)
@click.option("--steps", default=10000, type=int)
@click.option("--tmax", default=1, type=float)
@click.option("--tmin", default=0.011, type=float)
def sas(constructs, pool, build_files, pdbs, pool_path, **kwargs):
    log = logger.setup_applevel_logger()
    log.setLevel(logging.DEBUG)
    print(kwargs)
    w_opts = wrapper.Opts(extra_pdbs=pdbs)
    opts = simulate.Opts(
        wrapper_opts=w_opts,
        build_files_path=build_files,
        threads=kwargs["threads"],
    )
    df_constructs = pd.read_csv(constructs)
    df_pool = pd.read_csv(pool)
    s = simulate.Simulater(df_constructs, opts)
    s.setup()
    pool = optimize.PDBPool(df_pool, kwargs["ensemble_size"], pool_path)
    tsp = optimize.SampleandSelect(s, pool)
    params = AnnealerParams(kwargs["steps"], kwargs["tmax"], kwargs["tmin"])
    state, e = tsp.anneal(params)


@cli.command()
@click.argument("constructs")
@click.argument("pool")
@click.option("-bf", "--build-files", required=True)
@click.option("-pdbs", "--pdbs", required=True)
@click.option("-pp", "--pool-path", required=True)
@click.option("-t", "--threads", default=1, type=int)
@click.option("-es", "--ensemble-size", default=2, type=int)
@click.option("-r", "--runs", default=1, type=int)
def calc(constructs, pool, build_files, pdbs, pool_path, **kwargs):
    log = logger.setup_applevel_logger()
    log.setLevel(logging.DEBUG)
    w_opts = wrapper.Opts(extra_pdbs=pdbs, steps=10000000)
    opts = simulate.Opts(
        wrapper_opts=w_opts,
        build_files_path=build_files,
        threads=kwargs["threads"],
        runs=kwargs['runs']
    )
    df_constructs = pd.read_csv(constructs)
    df_pool = pd.read_csv(pool)
    s = simulate.Simulater(df_constructs, opts)
    s.setup()
    pool = optimize.PDBPool(df_pool, kwargs["ensemble_size"], pool_path)
    scorer = optimize.Scorer(df_constructs, 0)
    inds = range(0, len(df_pool))
    best = 1000
    count = 0
    for combo in itertools.combinations(inds, 2):
        ens_mems = pool.get_pdbs(combo)
        scores = s.score(ens_mems)
        score = scorer.score(scores, ens_mems)
        if best > score:
            best = score
        print(combo, score, best)


@cli.command()
@click.argument("constructs")
@click.argument("pool")
@click.option("-mem", "--member", multiple=True, required=True)
@click.option("-bf", "--build-files", required=True)
@click.option("-pp", "--pool-path", required=True)
@click.option("-pdbs", "--pdbs", required=True)
@click.option("-t", "--threads", default=1, type=int)
@click.option("--steps", default=1000000, type=int)
@click.option("-on", "--output-num", default=0, type=int)
@click.option("-r", "--runs", default=1, type=int)
def calc_one(constructs, pool, member, pdbs, build_files, pool_path, **kwargs):
    log = logger.setup_applevel_logger()
    log.setLevel(logging.DEBUG)
    w_opts = wrapper.Opts(extra_pdbs=pdbs, steps=kwargs["steps"])
    opts = simulate.Opts(
        wrapper_opts=w_opts,
        build_files_path=build_files,
        threads=kwargs["threads"],
        runs=kwargs["runs"]
    )
    df_pool = pd.read_csv(pool)
    df_constructs = pd.read_csv(constructs)
    ens = []
    for mem in member:
        index = df_pool.index[df_pool['pdb'] == mem].tolist()
        ens.append(index[0])
    s = simulate.Simulater(df_constructs, opts)
    s.setup()
    pool = optimize.PDBPool(df_pool, len(member), pool_path)
    scorer = optimize.Scorer(df_constructs, kwargs["output_num"])
    ens_mems = pool.get_pdbs(ens)
    scores = s.score(ens_mems)
    print(scores)
    score = scorer.score(scores, ens_mems)


if __name__ == "__main__":
    cli()
