import click
import pandas as pd

from rnamake_ens_gen import logger, score, optimize, wrapper


@click.command()
@click.argument("constructs")
@click.argument("pool")
@click.option("-bf", "--build-files", required=True)
@click.option("-pdbs", "--pdbs", required=True)
@click.option("-pp", "--pool-path", required=True)
def main(constructs, pool, build_files, pdbs, pool_path,**kwargs):
    logger.setup_applevel_logger()
    print(kwargs)
    w_opts = wrapper.Opts(extra_pdbs=pdbs)
    opts = score.Opts(
            wrapper_opts=w_opts,
            build_files_path=build_files,
            threads=5,
    )
    df_constructs = pd.read_csv(constructs)
    df_pool = pd.read_csv(pool)
    s = score.Scorer(df_constructs, opts)
    s.setup()
    pool = optimize.PDBPool(df_pool, 5, pool_path)
    tsp = optimize.SampleandSelect(s, pool)
    tsp.steps = 100
    tsp.Tmax = 0.01
    tsp.Tmin = 0.0001
    state, e = tsp.anneal()


if __name__ == '__main__':
    main()
