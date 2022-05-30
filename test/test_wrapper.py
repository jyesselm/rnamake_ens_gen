from rnamake_ens_gen import wrapper


def test_wrapper():
    path = "/Users/jyesselm/projects/python_packages/rnamake_ens_gen"
    opts = wrapper.Opts(extra_pdbs=f"{path}/test/resources/GAAA_tetraloop.pdb")
    seq = (
        "ACUGAUAAUAUAUGGGAUUGCGAACUACCUGCGAGUCUGGAAACAGACUCGAACUACAGGCGCAAUCCCU"
        "AAGUAUUAUUAGU"
    )[8:-8]
    build_file = path + "/test/resources/build_files/626.csv"
    ens_file = path + "/test/resources/ires.csv"
    w = wrapper.BuildMotifGraphWrapper()
    w.setup()
    count = w.run(seq, build_file, ens_file, opts)
    assert count < 10
