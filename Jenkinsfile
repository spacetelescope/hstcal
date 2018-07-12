if (utils.scm_checkout()) return

ASTROCONDA = "http://ssb.stsci.edu/astroconda"
OS_DEPS = "gcc gfortran libgomp pkgconfig"
CONDA_DEPS = "cmake cfitiso"
PREFIX = "/tmp/hstcal"

smoke = new BuildConfig()
smoke.nodetype = "linux"
smoke.name = 'Simple'
smoke.build_cmds = [
    "sudo yum install -y ${OS_DEPS}",
    "conda config --add channels ${ASTROCONDA}",
    "conda install ${CONDA_DEPS}",
    "mkdir _build",
    "cd _build && cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX}",
    "cd _build && make",
    "cd _build && make install",
]
smoke.test_cmds = [
    "for ex in ${PREFIX}/bin/*.e ; do \$ex --help || true; done",
]

utils.run([smoke])
