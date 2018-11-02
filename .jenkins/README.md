Setting up Jenkins-Bitbucket integration
----------------------------------------

These notes are for Jenkins v2.138.2 and for Bitbucket cloud circa November 2018.

Triggering tests on PRs requires the Jenkins plugin [Bitbucket Pull Request Builder](https://plugins.jenkins.io/bitbucket-pullrequest-builder).
Note that triggering tests on _commits_ requires a different Jenkins plugin
[Bitbucket Plugin](https://wiki.jenkins.io/display/JENKINS/BitBucket+Plugin).


Configuring a Jenkins job
-------------------------

On the project screen, select "New Item" and make a new "Pipeline" job. In the 
"configure" screen for this job, tick the box "Bitbucket Pull Requests Builder"
(this only appears if the Pull Request plugin is installed). In the Cron box,
type `* * * * *` (i.e. five space-separated asterisks). Fill in other details
as appropriate for the project. Also tick the box "Poll SCM".

In the "Pipeline" section, in the "definition" box, select "Pipeline script
from SCM". Under SCM, select "Git". In the "Repository URL", put
"https://bitbucket.org/gyrokinetics/gs2". Leave "credentials" as "none", and
for the branch specifier, put `*/jenkins`. Finally in the script path, put
`.jenkins/anvil.groovy`

Run the test once cleanly, by selecting "Build" on the left-hand side of in the
Job page.  Without doing this, automatic builds don't trigger!



Configuring the Bitbucket Webhook
---------------------------------

In the repository's settings page, select "Webhooks" from the contents list.
Click "Add webhook". Enter a descriptive title, and set the URL to
`https://<jenkins server address>/bitbucket-hook/`
For example, for Anvil
`https://anvil.softeng-support.ac.uk/jenkins/bitbucket-hook/`
Remember the trailing backslash on hook!

Under triggers, click "Choose from a full list of triggers", and then select
"Created" and "Updated" under Pull Request.  Click save.  Now the Jenkins job
will run automatically when a PR is created, or if new code is pushed to a PR.
The status and results of the job will display in the PR screen.



Using a local Jenkins server
----------------------------

For testing Jenkins jobs, it can be convenient to run a local instance of Jenkins.
Installing Jenkins from a package manager is very straightforward. To use your local
instance, in a browser go to the address `http://localhost:8080/`. Follow the on-screen
instructions to set up the server. When prompted, install all the recommended plugins,
and in addition select "Bitbucket Pull Request Builder".  Since the local
host cannot be seen by Bitbucket cloud, use the utility [ngrok](https://ngrok.com/)
to make a tunnel to your local host.  In a terminal, run `./ngrok http 8080`, and 
use copy the `ngrok.io` address generated in the "Forwarding" field:
`https://<hash>.ngrok.io`.  Use this as the server address in the Bitbucket webhook,
e.g. `https://2dd93fe4.ngrok.io/bitbucket-hook/`.  Now Bitbucket PRs will trigger 
builds on your local instance of Jenkins.  

Note that ngrok generated addresses expire after a few hours.
