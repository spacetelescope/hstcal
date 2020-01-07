// Obtain files from source control system.
if (utils.scm_checkout()) return

// Config data to share between builds.
def configure_cmd = "cmake -DCMAKE_INSTALL_PREFIX=./runtime"

// Define each build configuration, copying and overriding values as necessary.
bc0 = new BuildConfig()
bc0.nodetype = "python3.6"
bc0.name = "debug"
bc0.env_vars = ['PATH=./clone/_install/bin:runtime/bin:$PATH',
                'PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig',
                'LD_LIBRARY_PATH=runtime/lib',
                'LDFLAGS="-Wl,-rpath=runtime/lib:$CONDA_PREFIX/lib"',
                'OMP_NUM_THREADS=8']
bc0.conda_channels = ['http://ssb.stsci.edu/astroconda']
bc0.conda_packages = ['python=3.6',
                      'cfitsio',
                      'cmake',
                      'pkg-config']
bc0.build_cmds = ["${configure_cmd} -DCMAKE_BUILD_TYPE=Debug .",
                  "make install",
                  "calacs.e --version"]


bc1 = utils.copy(bc0)
bc1.name = "release"
// Would be nice if Jenkins can access /grp/hst/cdbs/xxxx directly.
bc1.conda_packages = ['python=3.6',
                     'ci-watson',
                     'cfitsio',
                     'cmake',
                     'pkg-config',
                     'pytest',
                     'requests',
                     'astropy']
bc1.build_cmds[0] = "${configure_cmd} -DCMAKE_BUILD_TYPE=RelWithDebInfo ."
bc1.test_cmds = ["pytest tests --basetemp=tests_output --junitxml results.xml --bigdata -v"]
bc1.failedUnstableThresh = 1
bc1.failedFailureThresh = 6


bc2 = utils.copy(bc0)
bc2.name = "optimized"
bc2.env_vars += ['CFLAGS="-O3"']
bc2.build_cmds[0] = "${configure_cmd} -DCMAKE_BUILD_TYPE=Release ."


// Iterate over configurations that define the (distibuted) build matrix.
// Spawn a host of the given nodetype for each combination and run in parallel.
utils.run([bc0, bc1, bc2])
