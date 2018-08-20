// Obtain files from source control system.
if (utils.scm_checkout()) return

// Config data to share between builds.
CFLAGS = ''
LDFLAGS = ''
DEFAULT_FLAGS = "${CFLAGS} ${LDFLAGS}"
// Some waf flags cause a prompt for input during configuration, hence the 'yes'.
configure_cmd = "yes '' | ./waf configure --prefix=./_install ${DEFAULT_FLAGS}"


// Define each build configuration, copying and overriding values as necessary.
bc0 = new BuildConfig()
bc0.nodetype = "linux-stable"
bc0.build_mode = "debug"
bc0.env_vars = ['PATH=./_install/bin:$PATH']
bc0.build_cmds = ["conda config --add channels http://ssb.stsci.edu/astroconda",
                  "conda install -q -y cfitsio pkg-config",
                  "${configure_cmd} --debug",
                  "./waf build",
                  "./waf install",
                  "calacs.e --version"]

bc1 = utils.copy(bc0)
bc1.build_mode = "release"
// Would be nice if Jenkins can access /grp/hst/cdbs/xxxx directly.
bc1.env_vars = ['PATH=./_install/bin:$PATH',
                'OMP_NUM_THREADS=8',
                'TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory/scsb-hstcal']
bc1.build_cmds = ["conda config --add channels http://ssb.stsci.edu/astroconda",
                  "conda install -q -y cfitsio pkg-config pytest requests astropy",
                  "pip install ci-watson",
                  "${configure_cmd} --release-with-symbols",
                  "./waf build",
                  "./waf install",
                  "calacs.e --version"]
bc1.test_cmds = ["pytest tests --basetemp=tests_output --junitxml results.xml --bigdata -v"]
bc1.failedUnstableThresh = 1
bc1.failedFailureThresh = 6


bc2 = utils.copy(bc0)
bc2.build_mode = "optimized"
bc2.build_cmds[2] = "${configure_cmd} --O3"


// Iterate over configurations that define the (distibuted) build matrix.
// Spawn a host of the given nodetype for each combination and run in parallel.
utils.run([bc0, bc1, bc2])
