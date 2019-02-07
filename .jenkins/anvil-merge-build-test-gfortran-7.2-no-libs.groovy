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
  choice(name: "opt", choices: "intermediate\ndebug\nfast\nwarnings\ncoverage", description: "Build type"),
  choice(name: "comms", choices: "serial\nmpi", description: "Communications paradigm"),
  string(name: "testexec", defaultValue: "mpirun -np 2", description: "Test command"),
  string(name: "modules", defaultValue: "gcc/7.2.0 python/3.6.0 openmpi/1.10.0", description: "Module environments to load"),
  string(name: "compiler", defaultValue: "gcc-7.2.0", description: "compiler-version, for consitent path names"),
  string(name: "GK_SYSTEM", defaultValue: "jenkins", description: "Module environments to load"),
  string(name: "fft", defaultValue: "fftw3", description: "Select FFT to use"),
  string(name: "baseURL", defaultValue: "https://gyrokinetics.slack.com/services/hooks/", description: "Base URL for Slack webhook")
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
String PRnameLink
String PRnameNoLink
String jobDescription
String jobDescriptionNoLink
String outcome

import groovy.json.JsonBuilder // For json construction.
  
// Do all the work on a node
node ('sl7') {

// PR name for display in Jenkins
PRname = currentBuild.description
// PR name for display in Slack
PRnameNoLink = "#${env.pullRequestId} ${env.pullRequestTitle}"
PRnameLink = "<https://${env.repo_host}/${env.repo_path}/pull-requests/${env.pullRequestId}|"+PRnameNoLink+">"

// Description of build and PR used in Slack messages
jobDescription = "Anvil build <${env.JOB_URL}/${env.BUILD_NUMBER}|#${env.BUILD_NUMBER}> for suite <${env.JOB_URL}|${env.JOB_NAME}> and PR "+PRnameLink
// Text-only version of jobDescription
jobDescriptionNoLink = "Anvil build #${env.BUILD_NUMBER} for suite ${env.JOB_NAME} and PR "+PRnameNoLink

// Set up for Slack API
withCredentials([string(credentialsId: 'slack-id', variable: 'token')]) { 
	slackURL = "${params.baseURL}jenkins-ci?token=${token}"
}

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
		'''
    	}

	stage('Merge') {

        	currentBuild.description = PRname + " - Merging PR"
		
		try{
			sh '.jenkins/merge.sh'
		} catch(gitError) {
		        def externalMethod = load("merge-error.groovy")
		        externalMethod.mergeError()
			throw gitError
		}
    	}
	
	stage('Build') {

        	currentBuild.description = "Building"
		
		try{
			sh """
				module load ${params.modules}
				export GK_SYSTEM=${params.GK_SYSTEM}
                              	export MAKEFLAGS=-IMakefiles USE_NETCDF= USE_FFT= USE_NEW_DIAG=
                              	export TESTEXEC='mpirun -np 6'

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
                              	export MAKEFLAGS=-IMakefiles USE_NETCDF= USE_FFT= USE_NEW_DIAG=
				pip install --user netCDF4
				make tests
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
		currentBuild.result = 'ABORTED'
  	}
  	else if (testFailures > 0) {
    		colourCode = '#FF0000'
    		currentBuild.description = "Build fails"
		outcome = " finished with ${testFailures} failures in ${runDuration}."
    		resultMessage = jobDescription+outcome
        	fallback = jobDescriptionNoLink+outcome // same message but without hyperlinks, for pop up notifications
        	endMessage = '' // No title
		currentBuild.result = 'FAILURE'
	}
  	else if (testFailures == 0) {
    		colourCode = '#32CD32' // lovely successful green
    		currentBuild.description = "Stable build"
		outcome = " finished successfully in ${runDuration}."
    		resultMessage = jobDescription+outcome
        	fallback = jobDescriptionNoLink+outcome // same message but without hyperlinks, for pop up notifications
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

