# Run with nosetests
import os
import platform
from astropy.io import fits
import shutil

path = os.path.dirname(__file__)
bin_path = os.path.join(path, "../../build.%s/hstio/test/" % platform.platform())
join = os.path.join

def test_hstio_resize_hdu():
    shutil.copyfile(join(path, "test_hstio_resize_hdu.fits"),
                    join(path, "test_hstio_resize_hdu_tmp.fits"))

    assert os.system("%s %s" % (join(bin_path, "test_hstio_resize_hdu"),
                                join(path, "test_hstio_resize_hdu_tmp.fits"))) == 0

    with fits.open(join(path, "test_hstio_resize_hdu_tmp.fits")) as hdulist:
        assert hdulist[1].data is not None
        print (hdulist[1].data[34, 12])
        assert hdulist[1].data[34,12] == 43

def test_hstio_resize_hdu2():
    shutil.copyfile(join(path, "test_hstio_resize_hdu.fits"),
                    join(path, "test_hstio_resize_hdu_tmp.fits"))

    assert os.system("%s %s" % (join(bin_path, "test_hstio_resize_hdu2"),
                                join(path, "test_hstio_resize_hdu_tmp.fits"))) == 0

    with fits.open(join(path, "test_hstio_resize_hdu_tmp.fits")) as hdulist:
        assert hdulist[1].data is not None
        assert hdulist[1].data[1000,50] == 50
        print (hdulist[1].data[42, 42])
        assert hdulist[1].data[42, 42] == 42

