import shutil
import subprocess
from dataclasses import dataclass


def does_program_exist(prog_name):
    if shutil.which(prog_name) is None:
        return False
    else:
        return True


@dataclass(frozen=True, order=True)
class Opts:
    steps: int = 1000000
    sterics: bool = False
    extra_pdbs: str = ""
    connect_str : str = "0,A222-A251"



class BuildMotifGraphWrapper(object):
    def __init__(self):
        pass

    def setup(self):
        if not does_program_exist("build_motif_graph"):
            raise ValueError("must have build_motif_graph in path!")

    def run(self, sequence, build_file, ens_file, opts : Opts):
        cmd = (
            f"build_motif_graph --build {build_file} --ensembles {ens_file } "
            f"--seq {sequence} "
        )
        if opts.extra_pdbs != "":
            cmd += f"--pdbs {opts.extra_pdbs} "
        if opts.connect_str != "":
            cmd += f"--connect \"{opts.connect_str}\" "
        if opts.steps != 1000000:
            cmd += f"--steps {opts.steps} "
        if opts.sterics:
            cmd += "--sterics "
        output = subprocess.check_output(cmd, shell=True).decode("utf8")
        lines = output.split("\n")
        return int(lines[-2])
