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
  string(name: "baseURL", defaultValue: "https://gyrokinetics.slack.com/services/hooks/", description: "Base URL for Slack webhook"),
  string(name: "TOKEN", defaultValue: "cvxSY36OQb3FswQVi0wcZcrg", description: "Token for Slack webhook"),
  string(name: "CODECOV_TOKEN", defaultValue: "a848e445-4d00-4aa9-8d00-e9950fcf29f2", description: "Token for Codecov upload")
]]])

/*
 End README
*/
  
// Do all the work on a node
node ('sl7') {

	stage('Checkout') {

        	currentBuild.description = "Checking out from repository"

        	// Old directory is kept
        	sh "rm -rf gs2"

		checkout([$class: 'GitSCM',  
     			branches: [[name: "${env.PULL_REQUEST_FROM_HASH}" ]],  
     			doGenerateSubmoduleConfigurations: false,  
             		submoduleCfg: [],  
			userRemoteConfigs: [[url: 'https://bitbucket.org/gyrokinetics/gs2']]]
		)
        	sh '''
			git config --global user.email "joseph.parker@stfc.ac.uk"
			git config --global user.name "Joseph Parker"
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

	export NETCDF_VERSION=4.1.3
	export NETCDFF_VERSION=4.4.3
	export GCC_VERSION=gcc-7.2.0
	export MY_NETCDF_DIR=${WORKSPACE}/software/\$GCC_VERSION/netcdf/\$NETCDF_VERSION
        export MY_NETCDFF_DIR=${WORKSPACE}/software/\$GCC_VERSION/netcdf-fortran/\$NETCDFF_VERSION
	export MY_NETCDF_INC=\$MY_NETCDF_DIR/include
	export MY_NETCDFF_INC=\$MY_NETCDFF_DIR/include
	export MY_NETCDF_LIB=\$MY_NETCDF_DIR/lib
	export MY_NETCDFF_LIB=\$MY_NETCDFF_DIR/lib

        export TEST_FILE_NETCDF=\$MY_NETCDF_DIR/include/netcdf.h
        export TEST_FILE_NETCDFF=\$MY_NETCDFF_DIR/include/netcdf.mod

        # true if either file missing
        if [[ -f \$TEST_FILE_NETCDF && -f \$TEST_FILE_NETCDFF ]]
        then 
              echo "netcdf + nectdff found"
        else
              module load gcc/7.2.0 openmpi/1.10.0
	      module avail


              rm -rf \$MY_NETCDF_DIR
              rm -rf \$MY_NETCDFF_DIR

              mkdir -p \$MY_NETCDF_DIR
              wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-\$NETCDF_VERSION.tar.gz
              tar xvf netcdf-\$NETCDF_VERSION.tar.gz
              cd netcdf-\$NETCDF_VERSION
              ./configure --prefix=\$MY_NETCDF_DIR --disable-netcdf-4
              make
              make install
              ls \$MY_NETCDF_INC
              export LD_LIBRARY_PATH=\$MY_NETCDF_LIB
              export CPPFLAGS=-I\$MY_NETCDF_INC
              export LDFLAGS=-L\$MY_NETCDF_LIB
              cd ..

              mkdir -p \$MY_NETCDFF_DIR
              wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-\$NETCDFF_VERSION.tar.gz
              tar xvf netcdf-fortran-\$NETCDFF_VERSION.tar.gz
              cd netcdf-fortran-\$NETCDFF_VERSION
              ./configure --prefix=\$MY_NETCDFF_DIR
              make
              make install

	      nc-config --all
	      nf-config --all
	fi
        """
  }

	stage('Build') {

        	currentBuild.description = "Building"
		
		try{
			sh '''
				export NETCDF_VERSION=4.1.3
				export NETCDFF_VERSION=4.4.3
				export GCC_VERSION=gcc-7.2.0
				export MY_NETCDF_DIR=${WORKSPACE}/software/\$GCC_VERSION/netcdf/\$NETCDF_VERSION
				export MY_NETCDFF_DIR=${WORKSPACE}/software/\$GCC_VERSION/netcdf-fortran/\$NETCDFF_VERSION
				export MY_NETCDF_INC=\$MY_NETCDF_DIR/include
				export MY_NETCDFF_INC=\$MY_NETCDFF_DIR/include
				export MY_NETCDF_LIB=\$MY_NETCDF_DIR/lib
				export MY_NETCDFF_LIB=\$MY_NETCDFF_DIR/lib

                              	export GK_SYSTEM=jenkins
                              	export MAKEFLAGS=-IMakefiles DEBUG=on
                              	export TESTEXEC='mpirun -np 6'

            			module load gcc/7.2.0 openmpi/1.10.0
            			module load fftw
                              	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:\$MY_NETCDF_LIB:\$MY_NETCDFF_LIB:/opt/software/gcc-7.2.0/fftw/3.3.5/lib

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
			sh '''
				export NETCDF_VERSION=4.1.3
				export NETCDFF_VERSION=4.4.3
				export GCC_VERSION=gcc-7.2.0
				export MY_NETCDF_DIR=${WORKSPACE}/software/\$GCC_VERSION/netcdf/\$NETCDF_VERSION
				export MY_NETCDFF_DIR=${WORKSPACE}/software/\$GCC_VERSION/netcdf-fortran/\$NETCDFF_VERSION
				export MY_NETCDF_INC=\$MY_NETCDF_DIR/include
				export MY_NETCDFF_INC=\$MY_NETCDFF_DIR/include
				export MY_NETCDF_LIB=\$MY_NETCDF_DIR/lib
				export MY_NETCDFF_LIB=\$MY_NETCDFF_DIR/lib

                              	export GK_SYSTEM=jenkins
                              	export MAKEFLAGS=-IMakefiles DEBUG=on
                              	export TESTEXEC='mpirun -np 6'
            			module load gcc/7.2.0 openmpi/1.10.0
            			module load fftw
                              	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:\$MY_NETCDF_LIB:\$MY_NETCDFF_LIB:/opt/software/gcc-7.2.0/fftw/3.3.5/lib

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

