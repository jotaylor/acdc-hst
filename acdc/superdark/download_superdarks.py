import glob
import gzip
import argparse
import wget
import os
import shutil

SUPERDARKS = {
    "superdark_FUVB_175_57650_58050_quiescent.asdf.gz": "https://stsci.box.com/shared/static/bhpjp43etyvt0mj7xye1mcj6t986qn0b.gz",
    "superdark_FUVB_175_56850_57100_active.asdf.gz": "https://stsci.box.com/shared/static/30bxwuh6k636nm14h08nu9z2a0m6162n.gz",
    "superdark_FUVB_169_56550_56800_active.asdf.gz": "https://stsci.box.com/shared/static/d5mllkh1ar7kl34o73twlccv6tzq3s5u.gz",
    "superdark_FUVB_169_57400_57700_quiescent.asdf.gz": "https://stsci.box.com/shared/static/0nl61crf5bu1g30bti3yn3hr912mjs1k.gz",
    "superdark_FUVB_163_58000_58300_quiescent.asdf.gz": "https://stsci.box.com/shared/static/xe88v32ggktwbadman10z9g40sywmyst.gz",
    "superdark_FUVB_167_57000_57300_active.asdf.gz": "https://stsci.box.com/shared/static/ual7xbt4bkmpsbsyt4kspt46fnv1yfoj.gz",
    "superdark_FUVB_167_55100_55400_quiescent.asdf.gz": "https://stsci.box.com/shared/static/b12bht6th8u5lja5v87sv6huc2jr5mlu.gz",
    "superdark_FUVB_163_56250_56420_quiescent.asdf.gz": "https://stsci.box.com/shared/static/madq09nk2rv5ugb5u2hro95nmzb0yzkp.gz",
    "superdark_FUVB_163_57050_57250_active.asdf.gz": "https://stsci.box.com/shared/static/xulxtemd8onwr8t5nb8x53va3zvf5ri5.gz",
    "superdark_FUVA_167_56520_56820_active.asdf.gz": "https://stsci.box.com/shared/static/7ifuaap3nnmm1w1jnpssxlcbasamkhnm.gz",
    "superdark_FUVA_171_57050_57300_active.asdf.gz": "https://stsci.box.com/shared/static/qug99msvg14rdqs6gqdujr276u8cov6h.gz",
    "superdark_FUVA_173_56900_57100_active.asdf.gz": "https://stsci.box.com/shared/static/5w1jxexyjdpduw9zrk4o6lh97trtb4p4.gz",
    "superdark_FUVA_167_57500_57800_quiescent.asdf.gz": "https://stsci.box.com/shared/static/0wq6duaao2ds9q5ayead3qm1oxo0gxj9.gz",
    "superdark_FUVA_163_58300_58600_active.asdf.gz": "https://stsci.box.com/shared/static/d1rmaqmu2zlzza0nc5b737d1g3icw39y.gz",
    "superdark_FUVA_169_55100_55400_quiescent.asdf.gz": "https://stsci.box.com/shared/static/n1yoeheto2ml8fgg5d98882fo8bzq6n7.gz",
    "superdark_FUVA_163_58000_58220_quiescent.asdf.gz": "https://stsci.box.com/shared/static/n8enh4s6pf2zk9w4z670js9h0l0eebhs.gz",
    "superdark_FUVA_169_55700_56010_active.asdf.gz": "https://stsci.box.com/shared/static/kosdxgsambik29cn5p618gupakpv8dv1.gz"}


def wget_superdarks(destdir):
    if not os.path.exists(destdir):
        os.makedirs(destdir)
    for item in SUPERDARKS:
        dest = os.path.join(destdir, item)
        wget.download(SUPERDARKS[item], dest)
    print(f"\n\n!!!\nIMPORTANT\n!!!\nMake sure to set the ACDC_SUPERDARKS environment variable:")
    print(f'export ACDC_SUPERDARKS="{destdir}"')
    print("For convenience, you may want to add this to your .bashrc/.cshrc/.bash_profile file")


def unzip_superdarks(destdir):
    zipped = glob.glob(os.path.join(destdir, "*gz"))
    for z_item in zipped:
        uz_item = z_item.split(".gz")[0]
        with gzip.open(z_item, "rb") as f_in, open(uz_item, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
            print("Uncompressing {} -> {}".format(z_item, uz_item))
        if os.path.isfile(uz_item):
            os.remove(z_item)


def download_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="destdir",
                        help="Directory to download superdarks into")
    args = parser.parse_args()
    wget_superdarks(args.destdir)
    unzip_superdarks(args.destdir)

if __name__ == "__main__":
    download_parser()
