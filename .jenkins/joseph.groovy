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
  string(name: "netcdff_path", defaultValue: "/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf-fortran/4.4.4", description: "Install path for netcdff"),
  string(name: "baseURL", defaultValue: "https://gyrokinetics.slack.com/services/hooks/", description: "Base URL for Slack webhook"),
  string(name: "TOKEN", defaultValue: "cvxSY36OQb3FswQVi0wcZcrg", description: "Token for Slack webhook")
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

import groovy.json.JsonBuilder // For json construction.
import groovy.json.JsonSlurper
  
// Do all the work on a node
node {

  /*

    Stage 1
    =======

    Checkout the source code.
    Here we use a clone for test purposes.

  */

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

/*
    stage('Slack start notification') {
        currentBuild.description = "Sending notification to slack"

	try {
		sh command
	} catch (slackError) {
		echo "Error posting to Slack"
	}

    }
*/

	stage('Checkout') {

        	currentBuild.description = "Checking out from repository"

        	// Old directory is kept
        	sh "rm -rf gs2"
        	//checkout([$class: 'GitSCM', branches: [[name: '*/jenkins']], userRemoteConfigs: [[url: 'https://bitbucket.org/gyrokinetics/gs2']]])

		echo "${sourceBranch}"
		echo "${targetBranch}"
		echo "${repositoryOwner}"
		echo "${repositoryName}"
		echo "${pullRequestId}"
		echo "${pullRequestTitle}"
		echo "${pullRequestAuthor}"


		checkout([$class: 'GitSCM',  
     			branches: [[name: "${env.PULL_REQUEST_FROM_HASH}" ]],  
     			doGenerateSubmoduleConfigurations: false,  
             		submoduleCfg: [],  
			userRemoteConfigs: [[url: 'https://bitbucket.org/gyrokinetics/gs2']]]
		)
        	sh '''
           		git submodule update --init --recursive --remote
           		cd Makefiles
           		git checkout jenkins
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
				git checkout ${targetBranch}
				git pull
				git submodule update --init
				git checkout -b tmp
				git merge ${sourceBranch}
			'''
		} catch(gitError) {
        		currentBuild.description = "PR merge failed"
			currentBuild.result = 'FAILURE'
		}

    	}
	

	stage('Build') {

        	currentBuild.description = "Building"
		
		try{
			sh '''
				export GK_SYSTEM=joseph-jenkins-gfortran
				make depend -I Makefiles
				make -I Makefiles
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
				export GK_SYSTEM=joseph-jenkins-gfortran
				make tests -I Makefiles NTESTPROCS=6
			'''
		} catch(testError) {
        		currentBuild.description = "Tests failed"
			currentBuild.result = 'FAILURE'
		}

    	}

}
