# Run with nosetests
import os
import platform
import pyfits
import shutil

path = os.path.dirname(__file__)
bin_path = os.path.join(path, "../../bin.%s" % platform.platform())

def test_hstio_resize_hdu():
    join = os.path.join
    
    shutil.copyfile(join(path, "test_hstio_resize_hdu.fits"),
                    join(path, "test_hstio_resize_hdu_tmp.fits"))

    os.system("%s %s" % (join(bin_path, "test_hstio_resize_hdu"),
                         join(path, "test_hstio_resize_hdu_tmp.fits")))

    hdulist = pyfits.open(join(path, "test_hstio_resize_hdu_tmp.fits"))
    assert hdulist[1].data is not None
    assert hdulist[1].data[34,12] == 43
    
