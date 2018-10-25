#!groovy
/*
  README
  ======

  Parameterised variables are used to define the build configuration.
  Properties used are string or choice, for a choice, the first entry is default.

*/
properties([[
  $class: "ParametersDefinitionProperty",
  parameterDefinitions: [
  string(name: "branch", defaultValue: "default", description: "Mercurial branch to use"),
  string(name: "label", defaultValue: "default", description: "Label for Slack fallback"),
  string(name: "repo_host", defaultValue: "bitbucket.org", description: "Domain for repository"),
  string(name: "repo_path", defaultValue: "gyrokinetics/gs2", description: "Relative path to repository"),
  choice(name: "repo_method", choices: "ssh\nhttps", description: "Protocol used to access repository"),
  string(name: "arch", defaultValue: "linux_x86_64_gfortran7.0", description: "Manually set build arch"),
  choice(name: "opt", choices: "intermediate\ndebug\nfast\nwarnings\ncoverage", description: "Build type"),
  choice(name: "comms", choices: "serial\nmpi", description: "Communications paradigm"),
  string(name: "modules", defaultValue: "gcc/7.2.0 fftw openmpi/1.10.0", description: "Module environments to load"),
  string(name: "mathLibs", defaultValue: "openblas", description: "Select BLAS/LAPACK to use"),
  string(name: "mathLibDir", defaultValue: "", description: "Specify maths library location"),
  string(name: "fft", defaultValue: "fftw3", description: "Select FFT to use"),
  string(name: "fftLibDir", defaultValue: "", description: "Specify fft library location"),
  string(name: "netcdff_path", defaultValue: "/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf-fortran/4.4.4", description: "Install path for netcdff")
]]])

// Mercurial repository to use
String repo = "${params.repo_method}://${params.repo_host}/${params.repo_path}" 

/*
 End README
*/

// Declare variables
String slackURL
String startMessage
String channel
String payload
String command
String colourCode
String result
String resultMessage
String fallback
String regressionDetails = ""
String netcdff_path
int numErrors = 0
int testFailures = 0
int testSkipped = 0
boolean buildSuccess = true

// Do all the work on a node
node('sl7') {

  /*

    Stage 1
    =======

    Checkout the source code.
    Here we use a clone for test purposes.

  */

  // Add a known host in advance so that the machine is recognized.
  sh "mkdir -p ~/.ssh; ssh-keyscan -t rsa ${params.repo_host} >> ~/.ssh/known_hosts"

  
    stage('Checkout') {

        currentBuild.description = "Checking out from repository"

        // Old directory is kept
        sh "rm -rf gs2"
        checkout([$class: 'GitSCM', branches: [[name: '*/jenkins']], userRemoteConfigs: [[url: 'https://bitbucket.org/gyrokinetics/gs2']]])
        sh """
           git submodule update --init --recursive --remote
           payload=${BITBUCKET_PAYLOAD}
           payloadbranch=$(echo $payload | jq -r .push.changes[].new.name)
           """
        
    }

  stage('Install netcdff') {

        currentBuild.description = "Installing netcdf/netcdff, if necessary"

        sh """

        export TEST_FILE_NETCDF=/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf/4.6.1/include/netcdf.h
        export TEST_FILE_NETCDFF=/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf-fortran/4.4.4/include/netcdf.mod

        # true if either file missing
        if [[ -f \$TEST_FILE_NETCDF && -f \$TEST_FILE_NETCDFF ]]
        then 
              echo "netcdf + nectdff found"
        else
              module load gcc/7.2.0 openmpi/1.10.0

              rm -rf /home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf/4.6.1
              rm -rf /home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf-fortran/4.4.4

              mkdir -p /home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf/4.6.1
              wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.6.1.tar.gz
              tar xvf netcdf-4.6.1.tar.gz
              cd netcdf-4.6.1
              ./configure --prefix=/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf/4.6.1 --disable-netcdf-4
              make
              make install
              ls /home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf/4.6.1/include
              export LD_LIBRARY_PATH=/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf/4.6.1/lib
              export CPPFLAGS=-I/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf/4.6.1/include
              export LDFLAGS=-L/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf/4.6.1/lib
              cd ..

              mkdir -p /home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf-fortran/4.4.4
              wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz
              tar xvf netcdf-fortran-4.4.4.tar.gz
              cd netcdf-fortran-4.4.4
              ./configure --prefix=${params.netcdff_path}
              make
              make install
	fi
        """
  }
  

    stage('Build') {

        currentBuild.description = "Building"
    
	sh """
                module load gcc/7.2.0 openmpi/1.10.0
                module load fftw
                module avail
                .jenkins/build-gcc-7.2.0.sh
           """
    }
        
    stage('Test') {

        currentBuild.description = "Running tests"
    
        sh """
            module load gcc/7.2.0 openmpi/1.10.0
            module load fftw
            .jenkins/run-tests.sh
        """
    }
}

