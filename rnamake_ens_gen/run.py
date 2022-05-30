import click

from rnamake_ens_gen import logger

def main():
    logger.setup_applevel_logger()

if __name__ == '__main__':
    main()