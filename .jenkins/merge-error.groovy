#!groovy

def mergeError(){

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
}

return this
