from ullyses.coadd import COSSegmentList

def coadd_cos_spectra(grating, indir, outdir=None, overwrite=True):
    """Coadd multiple 1-D COS spectra of the same grating.
    Args:
        grating (str): COS grating name.
        indir (str): Path that holds all data to be coadded. All data in this
            directory that match grating will be coadded.
    Returns:
        None
    """

    prod = COSSegmentList(grating, path=indir)
    if len(prod.members) > 0:
        prod.create_output_wavelength_grid()
        prod.coadd()
        prod.target = prod.get_targname()
        prod.targ_ra, prod.targ_dec = prod.get_coords()
        if outdir is None:
            outdir = indir
        outname = os.path.join(outdir, f"{prod.target.lower()}_{prod.grating.lower()}.fits")
        prod.write(outname, overwrite)
        print(f"Wrote {outname}")
