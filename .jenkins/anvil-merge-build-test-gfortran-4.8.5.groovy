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
  string(name: "label", defaultValue: "default", description: "Label for Slack fallback"),
  string(name: "repo_host", defaultValue: "bitbucket.org", description: "Domain for repository"),
  string(name: "repo_path", defaultValue: "gyrokinetics/gs2", description: "Relative path to repository"),
  choice(name: "repo_method", choices: "ssh\nhttps", description: "Protocol used to access repository"),
  string(name: "arch", defaultValue: "linux_x86_64_gfortran7.0", description: "Manually set build arch"),
  choice(name: "opt", choices: "intermediate\ndebug\nfast\nwarnings\ncoverage", description: "Build type"),
  choice(name: "comms", choices: "serial\nmpi", description: "Communications paradigm"),
  string(name: "testexec", defaultValue: "mpirun -np 2", description: "Test command"),
  string(name: "modules", defaultValue: "gcc/4.8.5 fftw python/3.6.0 openmpi/1.10.0", description: "Module environments to load"),
  string(name: "compiler", defaultValue: "gcc-7.2.0", description: "compiler-version, for consitent path names"),
  string(name: "GK_SYSTEM", defaultValue: "jenkins-4.8.5", description: "Value of GK_SYSTEM, which determines the local Makefile used"),
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
String PRname

import groovy.json.JsonBuilder // For json construction.
  
// Do all the work on a node
node ('sl7') {

// Store nicely formatted orginal name (easier than working out how it is passed to Jenkins!)
PRname = currentBuild.description

// Set up for Slack API
slackURL = "${params.baseURL}jenkins-ci?token=${env.TOKEN}"

startMessage = "Anvil build #${env.BUILD_NUMBER} for ${env.JOB_NAME}. (<${env.JOB_URL}|View>)"

// This fallback message shows in a Slack notification

fallback = "gs2 ${params.label} \"${params.opt}\" ${params.comms} build started for branch ${params.branch}"

channel = "#jenkins"
colourCode = "#439FE0" // Jenkins blue

// Construct the json payload
def mrkdwn = [ "pretext", "text", "fields" ]  // Required for Slack API
def attachment = [pretext: fallback,
                  text: startMessage,
		  fallback: fallback,
		  color: colourCode,
		  mrkdwn_in: mrkdwn]

def attachments = [ attachment ]
def load = [ channel: channel, attachments: attachments ]

try {

  def json = new JsonBuilder(load) // This may fail if the script isn't whitelisted
  payload = json.toString()
  json = null // Need to nullify the Json object as it is not "serializable" by Java

}

catch (jsonError) {

  echo "Error constructing json for Slack message - likely a script whitelist issue"



  // Construct the json as a string - much easier to get wrong, much harder to debug

  payload = "{\"channel\":\"#jenkins\",\"attachments\":[{\"pretext\":\"${fallback}\",\"text\":\"${startMessage}\",\"fallback\":\"${fallback}\",\"color\":\"#439FE0\",\"mrkdwn_in\":[\"pretext\",\"text\",\"fields\"]}]}"

}

// Construct the command as a string so that variables are expanded properly
command = "set +x; curl -X POST -H 'Content-Type: application/json' -d '${payload}' ${slackURL}"

// Set up netcdf paths
String NETCDF_VERSION="4.6.1"
String NETCDFF_VERSION="4.4.4"
String MY_NETCDF_DIR="${WORKSPACE}/software/${params.compiler}/netcdf/${NETCDF_VERSION}"
String MY_NETCDFF_DIR="${WORKSPACE}/software/${params.compiler}/netcdf-fortran/${NETCDFF_VERSION}"
String MY_NETCDF_INC="${MY_NETCDF_DIR}/include"
String MY_NETCDFF_INC="${MY_NETCDFF_DIR}/include"
String MY_NETCDF_LIB="${MY_NETCDF_DIR}/lib"
String MY_NETCDFF_LIB="${MY_NETCDFF_DIR}/lib"
String libpath="${MY_NETCDF_LIB}:${MY_NETCDFF_LIB}:/opt/software/gcc-7.2.0/openmpi/1.10.0/lib:/opt/software/gcc-7.2.0/fftw/3.3.5/lib:/opt/software/gcc/7.2.0/lib64:/opt/software/gcc/7.2.0/lib"


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
           		git checkout jenkins
			git pull
           		cd Makefiles
           		git checkout jenkins
			git pull
           		cd ..

		'''
    	}

	stage('Merge') {

        	currentBuild.description = PRname + " - Merging PR"
		
		try{
			sh '.jenkins/merge.sh'
		} catch(gitError) {
        		currentBuild.description = PRname + " - PR merge failed"
			currentBuild.result = 'FAILURE'
	
			runDuration=currentBuild.durationString[0..-14]
        		endMessage = ""
			colourCode = '#FF0000' // failure red
    			resultMessage = "Anvil build <${env.JOB_URL}/${env.BUILD_NUMBER}|#${env.BUILD_NUMBER}> for <${env.JOB_URL}|${env.JOB_NAME}> failed in git merge step after ${runDuration}"
    			fallback = "Anvil build #${env.BUILD_NUMBER} for ${env.JOB_NAME} failed in git merge step after ${runDuration}"

  			attachment = [text: resultMessage,
				pretext: endMessage,
				fallback: fallback,
				color: colourCode,
				mrkdwn_in: mrkdwn]
  			attachments = [ attachment ]
  			load = [ channel: channel, attachments: attachments ]
  			try {
    				endJson = new JsonBuilder(load)
    				payload = endJson.toString()
 		   		endJson = null // Need to nullify the Json object as it is not "serializable" by Java
 		 	} catch (jsonError) {
    				echo "Error contructing final Slack message - likely a library permissions issue"
		    		payload = "{\"channel\":\"#jenkins\",\"attachments\":[{\"Job Ended\",\"text\":\"${resultMessage}\",\"fallback\":\"Job Ended\",\"color\":\"${colourCode}\",\"mrkdwn_in\":[\"pretext\",\"text\",\"fields\"]}]}"
  			}
  			command = "set +x; curl -X POST -H 'Content-Type: application/json' -d '${payload}' ${slackURL}"
			try {
				sh command
			} catch (Exception slackError) {
				echo "Error posting to Slack"
			}
			throw gitError
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
              ./configure --prefix=/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1
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
              ./configure --prefix=/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf-fortran/4.4.4
              make
              make install
	fi
        """
  }

	stage('Build') {

        	currentBuild.description = "Building"
		
		try{
			sh """
				export GK_SYSTEM=${params.GK_SYSTEM}
                              	export MAKEFLAGS=-IMakefiles
                              	export TESTEXEC='mpirun -np 2'

            			module load gcc/4.8.5 openmpi/1.10.0
            			module load fftw
                              	export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1/lib:/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf-fortran/4.4.4/lib:/opt/software/gcc-system/fftw/3.3.5/lib

				make depend
				make test_make
				make -j
			"""
		} catch(buildError) {
        		currentBuild.description = PRname + " - build failed"
			currentBuild.result = 'FAILURE'
			throw buildError
		}

    	}

	stage('Test') {

        	currentBuild.description = "Testing"
		try{
			sh """
				export TESTEXEC="${params.testexec}"
				module load ${params.modules}
				export GK_SYSTEM=${params.GK_SYSTEM}
                              	export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf/4.6.1/lib:/home/jenkins/workspace/Gyrokinetics/Merge+Build+Test-gfortran-4.8.5/software/gcc-4.8.5/netcdf-fortran/4.4.4/lib:/opt/software/gcc-system/fftw/3.3.5/lib
                              	export MAKEFLAGS=-IMakefiles
				pip install --user netCDF4
				make tests TEST_LEVEL=2
				cd ../
				rm -rf gs2
			"""
		} catch(testError) {
        		currentBuild.description = PRname + " - tests failed"
			currentBuild.result = 'FAILURE'
			throw testError
		}
    	}

    stage('Finishing') {

        currentBuild.description = "Parsing test results"

	// This doesn't work, need a better way to parse files
	//sh '''
	//	export params.testFailures=$(grep 'failures="1"' test-results/*.xml | wc -l)
	//'''
	testFailures = sh ( script: """ grep 'failures="1"' test-results/*.xml | wc -l """, returnStdout: true ) as Integer

        currentBuild.description = "Sending notification to slack"
	
	// "currentBuild.durationString" unhelpfully returns the time followed by "and counting" because we're still in the pipeline. Bah.
	runDuration=currentBuild.durationString[0..-14]
        endMessage = ""
	colourCode = '#FF0000' // failure red
    	resultMessage = "Anvil build <${env.JOB_URL}/${env.BUILD_NUMBER}|#${env.BUILD_NUMBER}> for <${env.JOB_URL}|${env.JOB_NAME}> finished in mysterious circumstances."
    	fallback = "Anvil build #${env.BUILD_NUMBER} for ${env.JOB_NAME} finished in mysterious circumstances."

	result = currentBuild.result ?: 'SUCCESS'

  	if (result == 'ABORTED') {
    		colourCode = '#808080' // anonymous grey
    		currentBuild.description = "Build aborted"
    		resultMessage = "Anvil build <${env.JOB_URL}/${env.BUILD_NUMBER}|#${env.BUILD_NUMBER}> for <${env.JOB_URL}|${env.JOB_NAME}> was aborted."
        	fallback = "Anvil build #${env.BUILD_NUMBER} for ${env.JOB_NAME} was aborted." // same message but without hyperlinks, for pop up notifications
		currentBuild.result = 'STOPPED'
  	}
  	else if (testFailures > 0) {
    		colourCode = '#FF0000'
    		currentBuild.description = "Build fails"
    		resultMessage = "Anvil build <${env.JOB_URL}/${env.BUILD_NUMBER}|#${env.BUILD_NUMBER}> for <${env.JOB_URL}|${env.JOB_NAME}> finished with ${testFailures} failures in ${runDuration}."
        	fallback = "Anvil build #${env.BUILD_NUMBER} for ${env.JOB_NAME} finished with ${testFailures} failures in ${runDuration}." // same message but without hyperlinks, for pop up notifications
        	endMessage = '' // No title
		currentBuild.result = 'FAILURE'
	}
  	else if (testFailures == 0) {
    		colourCode = '#32CD32' // lovely successful green
    		currentBuild.description = "Stable build"
    		resultMessage = "Anvil build <${env.JOB_URL}/${env.BUILD_NUMBER}|#${env.BUILD_NUMBER}> for <${env.JOB_URL}|${env.JOB_NAME}> finished successfully in ${runDuration}."
        	fallback = "Anvil build #${env.BUILD_NUMBER} for ${env.JOB_NAME} finished successfully in ${runDuration}." // same message but without hyperlinks, for pop up notifications
        	endMessage = '' // No title
	}
  	else {
		println "Something is wrong with our counting"
  	}

  	attachment = [text: resultMessage,
		pretext: endMessage,
		fallback: fallback,
		color: colourCode,
		mrkdwn_in: mrkdwn]
  	attachments = [ attachment ]
  	load = [ channel: channel, attachments: attachments ]
  	try {
    		endJson = new JsonBuilder(load)
    		payload = endJson.toString()
    		endJson = null // Need to nullify the Json object as it is not "serializable" by Java
  	} catch (jsonError) {
    		echo "Error contructing final Slack message - likely a library permissions issue"
    		payload = "{\"channel\":\"#jenkins\",\"attachments\":[{\"Job Ended\",\"text\":\"${resultMessage}\",\"fallback\":\"Job Ended\",\"color\":\"${colourCode}\",\"mrkdwn_in\":[\"pretext\",\"text\",\"fields\"]}]}"
  	}
  	command = "set +x; curl -X POST -H 'Content-Type: application/json' -d '${payload}' ${slackURL}"
	try {
		sh command
	} catch (slackError) {
		echo "Error posting to Slack"
	}

	// Now everything is finished, set the name back to a useful hyperlink:
	currentBuild.description = PRname

    }
}

