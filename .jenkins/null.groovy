#!groovy
/*
  README
  ======

  The Jenkins Bitbucket Pull Request Builder plugin requires that the Jenkins 
  job be manually triggered and run successfully before PRs trigger jobs.
  However, our main job fails if triggered manually, because it tries to read
  information from the PR that triggered it.  Therefore, we need to initialize
  the Jenkins jobs with a blank groovy file that will execute successfully.

*/

/*
 End README
*/
  
// Do all the work on a node
node {
	stage('Null') {

    	}
}

