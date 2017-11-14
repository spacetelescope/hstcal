def args = [
    ['debug', '--debug'],
    ['release', '--release-with-symbols'],
    ['optimized', '--O3']
]
def defaults = 'CFLAGS="-m64" LDFLAGS="-m64"'
def name = ''
def option = ''

for(int i = 0; i < args.size(); i++)
{
    name = args[i][0]
    option = args[i][1]

    node {
        stage("Checkout") {
            checkout scm
        }

        def prefix = pwd() + '/_install'
        def runtime = ["PATH=${prefix}/bin:${env.PATH}"]

        stage("System (${name})") {
            sh 'uname -a'
            sh 'lscpu'
            sh 'free -m'
            sh 'df -hT'
        }
        stage("Configure (${name})") {
            sh "yes '' | ./waf configure --prefix=${prefix} ${defaults} ${option}"
        }
        stage("Build (${name})") {
            sh './waf build'
        }
        stage("Install (${name})") {
            sh './waf install'
        }
        try {
            stage("Test (${name})") {
                sh 'conda install -q -y pytest astropy'
                withEnv(runtime) {
                    sh 'pytest -s --basetemp=tests_output --junitxml results.xml --remote-data tests'
                }
            }
        }
        finally {
            stage("Ingest (${name}") {
                junit '*.xml'
                step([$class: 'XUnitBuilder',
                    thresholds: [
                    [$class: 'SkippedThreshold', failureThreshold: '0'],
                    [$class: 'FailedThreshold', failureThreshold: '2']],
                    tools: [[$class: 'JUnitType', pattern: '*.xml']]])
            }
        }
    }
}
