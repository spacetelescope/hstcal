// Not working? 1
pipeline {
    agent any
    stages {
        stage('Configure') {
            steps {
                sh './waf configure --prefix=$(pwd)/_install CFLAGS="-m64" LDFLAGS="-m64"'
            }
        }
        stage('Build') {
            steps {
                sh './waf build'
            }
        }
        stage('Install') {
            steps {
                sh './waf install'
            }
        }
    }
}
