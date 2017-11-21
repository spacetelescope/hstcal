def nodes = ['linux']
def release_modes = ['debug', 'optimized', 'release']
def release_args = ['--debug', '--O3', '--release-with-symbols']

def CFLAGS = "CFLAGS=\"-m64\""
def LDFLAGS = "LDFLAGS=\"-m64\""
def DEFAULT_FLAGS = "${CFLAGS} ${LDFLAGS}"
def tasks = [:]

for(int i = 0; i < nodes.size(); i++)
{
    def node_is = nodes[i]

    for(int j = 0; j < release_modes.size(); j++)
    {
        def name = release_modes[j]
        def option = release_args[j]

        // By default, builds here only test for successful compilation.
        def test_cmd = 'calacs.e --version'

        // One of them should run "remote_data" tests when in release mode.
        // "slow" tests should only run in nightly build and not addressed here.
        if (name == 'release') {
            test_cmd = 'pytest tests --basetemp=tests_output --junitxml results.xml --remote-data'
        }

        tasks["${node_is}/${name}"] = {
            node {
                stage("Checkout") {
                    checkout scm
                }

                def prefix = pwd() + '/_install'
                def runtime = ["PATH=${prefix}/bin:${env.PATH}"]

                stage("Generate (${name})") {
                    sh "yes '' | ./waf configure --prefix=${prefix} ${option} ${DEFAULT_FLAGS}"
                    sh './waf build'
                    sh './waf install'
                }
                try {
                    stage("Test (${name})") {
                        if (name == 'release') {
                            sh 'conda install -q -y pytest requests astropy'
                            sh 'pip install -q pytest-remotedata'
                        }
                        withEnv(runtime) {
                            sh "${test_cmd}"
                        }
                    }
                }
                finally {
                    step([$class: 'XUnitBuilder',
                        thresholds: [
                        [$class: 'SkippedThreshold', failureThreshold: '0'],
                        [$class: 'FailedThreshold', unstableThreshold: '1'],
                        [$class: 'FailedThreshold', failureThreshold: '6']],
                        tools: [[$class: 'JUnitType', pattern: '*.xml']]])
                }
            }
        }
    }
}

stage("Matrix") {
    parallel tasks
}
