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
  string(name: "netcdff_path", defaultValue: "/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf-fortran/4.4.4", description: "Install path for netcdff"),
  string(name: "baseURL", defaultValue: "https://gyrokinetics.slack.com/services/hooks/", description: "Base URL for Slack webhook"),
  string(name: "TOKEN", defaultValue: "cvxSY36OQb3FswQVi0wcZcrg", description: "Token for Slack webhook"),
  string(name: "CODECOV_TOKEN", defaultValue: "a848e445-4d00-4aa9-8d00-e9950fcf29f2", description: "Token for Codecov upload")
]]])

/*
 End README
*/
  
// Do all the work on a node
node {

	stage('Checkout') {

        	currentBuild.description = "Checking out from repository"

        	// Old directory is kept
        	sh "rm -rf gs2"

		checkout([$class: 'GitSCM',  
     			branches: [[name: "next" ]],  
     			doGenerateSubmoduleConfigurations: false,  
             		submoduleCfg: [],  
			userRemoteConfigs: [[url: 'https://bitbucket.org/gyrokinetics/gs2']]]
		)
        	sh '''
           		git submodule update --init --recursive --remote
           		cd Makefiles
           		git checkout jenkins
			git pull
           		cd ..

		'''
    	}

	stage('Merge') {

        	currentBuild.description = "Merging PR"
		
		try{
			sh '''
				git fetch
				git checkout ${sourceBranch}
				git pull
				git checkout origin/${targetBranch}
				git checkout -b jobid-${BUILD_NUMBER}-${sourceBranch}-${targetBranch}
				git merge ${sourceBranch}
				rm -rf tests/unit_tests/collisions/
			'''
		} catch(gitError) {
        		currentBuild.description = "PR merge failed"
			currentBuild.result = 'FAILURE'
		}

    	}
	
  stage('Install netcdff') {

        currentBuild.description = "Installing netcdf/netcdff, if necessary"

        sh """

        export TEST_FILE_NETCDF=/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1/include/netcdf.h
        export TEST_FILE_NETCDFF=/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf-fortran/4.4.4/include/netcdf.mod

        # true if either file missing
        if [[ -f \$TEST_FILE_NETCDF && -f \$TEST_FILE_NETCDFF ]]
        then 
              echo "netcdf + nectdff found"
        else
              module load gcc/4.8.5 openmpi/1.10.0

              rm -rf /home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1
              rm -rf /home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf-fortran/4.4.4

              mkdir -p /home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1
              wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.6.1.tar.gz
              tar xvf netcdf-4.6.1.tar.gz
              cd netcdf-4.6.1
              ./configure --prefix=/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1 --disable-netcdf-4
              make
              make install
              ls /home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1/include
              export LD_LIBRARY_PATH=/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1/lib
              export CPPFLAGS=-I/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1/include
              export LDFLAGS=-L/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1/lib
              cd ..

              mkdir -p /home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf-fortran/4.4.4
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
		
		try{
			sh '''
                              	export GK_SYSTEM=jenkins-4.8.5
                              	export MAKEFLAGS=-IMakefiles DEBUG=on
                              	export TESTEXEC='mpirun -np 6'

            			module load gcc/4.8.5 openmpi/1.10.0
            			module load fftw
                              	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1/lib:/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf-fortran/4.4.4/lib:/opt/software/gcc-4.8.5/fftw/3.3.5/lib

				make depend
				make test_make
				make -j
			'''
		} catch(buildError) {
        		currentBuild.description = "Build failed"
			currentBuild.result = 'FAILURE'
		}

    	}

	stage('Test') {

        	currentBuild.description = "Testing"
		
		try{
			sh '''#!/bin/bash
                              	export GK_SYSTEM=jenkins-4.8.5
                              	export MAKEFLAGS=-IMakefiles DEBUG=on
                              	export TESTEXEC='mpirun -np 6'

            			module load gcc/4.8.5 openmpi/1.10.0
            			module load fftw
                              	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1/lib:/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf-fortran/4.4.4/lib:/opt/software/gcc-4.8.5/fftw/3.3.5/lib

				make tests
				cd ../
				rm -rf gs2
			'''
		} catch(testError) {
        		currentBuild.description = "Tests failed"
			currentBuild.result = 'FAILURE'
		}

    	}

}

