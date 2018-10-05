# Contributing to GS2

GS2 is Free and Open Source and welcomes contributions from any
users. This file contains some help and guidelines on how to
effectively contribute to the GS2 project. There is also some help for
developers using git for the first time, those who have previously
used GS2 when it was under subversion, or for those who need a
refresh.

## Quickstart

A brief overview of how to contribute:

- Clone the GS2 repository

        git clone git@bitbucket.org:gyrokinetics/gs2.git

- Do your development in a new branch starting at `next`

        git checkout next
        git pull origin next
        git checkout -b my_new_branch

- Make sure the submodules are up-to-date

        git submodule update --init --recursive

- Commit your changes

        git add <files>
        git commit -m "Add my new feature"

- Push your changes to bitbucket

        git push origin my_new_branch

- Open a pull request on bitbucket: https://bitbucket.org/gyrokinetics/gs2/pull-requests/


## First steps with git

Getting a local copy of a repository in git is called cloning. The
command is:

    git clone https://bitbucket.org/gyrokinetics/gs2.git

This will clone the GS2 repository into a new directory called
`gs2`. This is called the working copy (or sometimes working tree).

We can check the status of the working copy using:

    git status

which will produce some output like this:

    $ git status
    On branch master
    Your branch is up to date with 'origin/master'.

    nothing to commit, working tree clean

Like `svn status`, the `git status` command will report tracked files
with changes as well as untracked files. To commit changes to files,
we need to first add the changes to a staging area:

    git add <file>...

Where `<file>...` is one or more filenames. Whereas `svn add` is only
used to start tracking a file, `git add` is used every time a file is
changed. Staging changes like this does not actually make a
commit. You might think of `git add` as adding the changes to a box,
where you are free to add and remove those changes as you wish until
you are happy with the contents of the box, when you then commit the
whole box of changes at once.

In general, `git status` will report something like the following:

    On branch master
    Your branch is up to date with 'origin/master'.

    Changes to be committed:
     (use "git reset HEAD <file>..." to unstage)

        new file:   MyFile.txt
        modified:   OtherFile.txt

    Changes not staged for commit:
     (use "git add <file>..." to update what will be committed)
     (use "git checkout -- <file>..." to discard changes in working directory)

           modified:   SomeFile.txt

    Untracked files:
     (use "git add <file>..." to include in what will be committed)

           NewFile.txt

As you can see, `git status` gives helpful reminders of common
operations.

Once the "Changes to be committed" section looks how you want it,
commit the changes:

    git commit

This will pop up your default text editor to write a commit
message. After saving the message and closing the editor, your changes
will be committed.

A tip for writing commit messages: it’s a strong convention to write
your commit messages as if they would complete the sentence "this
commit will…", e.g.:

    Add new method to calculate energy
    Fix bug with calculating potential

See the Resources section at the bottom for some more tips.

Committing changes in git is not like svn -- the commit stays on your
local machine until you are ready to push the changes:

    git push origin <branch>

## New Development Workflow

GS2 now uses a development model similar to the "Gitflow" model:
https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow
This might look a bit complicated at first, but in practice is very
easy to work with.

We have two "special" branches:

- `master` – the stable branch
- `next` – the development branch

These two branches are locked and cannot be modified
directly. Instead, work is initially done in feature branch starting
at next. When the feature or work is complete, a pull request (PR) is
opened to request the branch be merged into next. When we want to
release a new version, next is merged into master. Bug-fixes are a
special case and PRs can be created to merge them straight into
master.

The reason for working like this is that it allows us to maintain a
"known good", stable version in master, while allowing development to
carry on in next. We require pull requests in order to give the
community a chance to review the new code. This has several benefits:
firstly, it allows the maintainers to ensure the quality of the
codebase is kept high, and give less experienced developers feedback
on their code. It also has an important benefit of making sure more
than one person knows how a given part of the code works. Lastly, we
can automatically run tests on submitted code when a PR is opened. All
this is to say that the pull request model allows us to increase trust
in the GS2 source code while maintaining an open and transparent
development process.

Once a PR has been reviewed and accepted, a maintainer will approve it
and merge the code. Any member of the community may review a PR and
offer constructive feedback, but final approval must be from a member
of the Governing Committee (GC).

PRs should ideally provide tests to prove that they work correctly,
and new features must include documentation. Adding tests can
sometimes be tricky, so members of the GC will be happy to help with
these.

If you have questions about any part of the development workflow,
again, please don’t hesitate to ask a member of the GC, raise an issue
on bitbucket (https://bitbucket.org/gyrokinetics/gs2), or ask in the
Slack channel (https://gyrokinetics.slack.com).

We encourage developers to use one branch per feature or bug
fix. Multiple features in a branch make PRs difficult to review, as it
can be hard to see how changes are related to one another.

We strongly encourage commits to be atomic, i.e., contain one,
irreducible piece of code development. Small, orthogonal commits with
good messages act as an incredibly vital piece of documentation.
Example of the new development workflow Here is a complete example of
how to actually work with this new workflow, and get changes into GS2:

1. Starting from the next branch, pull the latest changes:

        git pull origin next

2. Make and checkout a new branch

        git checkout -b my_branch

3. Start implementing your new feature or bug fix
4. Commit your changes with a good commit message:

        git add <file>...
        git commit -m "Add new feature"

5. Push your branch to bitbucket

        git push -u origin my_branch

6. Open a pull request (PR) on bitbucket

    - With the GS2 repository open, click the "+" icon the left-hand
      sidebar and select "Create a Pull Request" under "Get to work"
    - Select your branch from the left drop-down box, and which branch
      you want to merge it into from the right box. You will normally
      want next
    - The title of the PR will be pre-populated with the name of the
      branch, while the body will contain the commit messages from the
      branch (another good reason to write commit messages carefully!)
    - Optionally: request a reviewer(s). You may like to consider one
      person who knows the section of code well and a second who
      doesn’t, to spread the knowledge around

## Making changes to the submodules

GS2 currently has two submodules, `utils` and `Makefiles`. These are
separate repositories shared with the rest of the Gyrokinetics
project. A git submodule can be thought of as a bit like a pointer to
a particular commit in another repository. We store a specific version
of `utils` and `Makefiles` in the gs2 repo. In order to bring a change
to one of the submodules into gs2, it is necessary to update the
particular version we store.

To update `utils`/`Makefiles`, you need to first submit the change to
the relevant main repo. After that has been approved and merged, you
can then update the submodule in gs2:

    cd <utils|Makefiles>        # Go into relevant submodule
    git checkout master         # Make sure you are on the master branch
    git pull origin master      # Get the latest changes
    cd ..                       # Go back to gs2
    git add <utils|Makefiles>   # Add the updated submodule version

If someone else has updated one of the submodules, you may need to
bring in those changes:

    git submodule update --init --recursive

In recent versions of `git`, you can run the following command to make
`git` pull updated submodules automatically:

    git config --global submodule.recurse true

### A note on developing the submodules

The default protocol of the submodules is HTTPS (rather than
SSH). This is necessary to enable GS2 users (non-developers) to
download the code successfully. Therefore, developers who wish to
modify code belonging to `utils` or `Makefiles` will either have to:

- clone a separate copy of `utils` and `Makefiles` via the git
  protocol, make changes to those, push those changes and then pull
  them into their `gs2` working copy; or

- add `git config --global url.git@bitbucket.org:gyrokinetics/`
  instead of `https://bitbucket.org/gyrokinetics/` to their git
  configuration before cloning `gs2`, thus allowing them to make and
  push changes directly from the `utils` and `Makefiles`
  subdirectories within their `gs2` working copy. You will still need
  to add the new submodule commit to `gs2`!

See the official git documentaion e-book for further details on
working with sub-modules:
https://git-scm.com/book/en/v2/Git-Tools-Submodules


## Moving uncommitted changes from an subversion checkout to git

If you have uncommitted changes in your local working copy of the
subversion repository, don’t worry, you can follow these simple steps
below to get them under git. If your changes are based on a branch
that’s up to date, do the following:

1. Create a new branch in git starting at the branch you were working
   from:

        git checkout -b <new branch name> <old branch>

2. Copy the modified files from your working copy of the subversion
   repository on top of the ones in your git working copy
3. Carry on from step 4 above

However, if your changes are based on a branch that’s not up-to-date
with the svn one, you’ll have to make a new branch based on the last
commit you do have. This will make the rest of the process much
cleaner.

To find the corresponding commit in the git repository, do

    git log --all --grep "@<svn commit> "

replacing `<svn commit>` with just the number of your svn commit,
e.g. if your commit is r123, run git log --grep "@123 " Note that you
may need a space after the svn commit ID to find the exact one. Make a
note of the git commit hash and then use it instead of the old branch
name in the steps above.

If you want to be really nice to the maintainers, you could split your
changes into several smaller commits. You could do this either with a
git tool, or a diff tool to copy parts of your changes instead of all
of them at once.

Please don’t hesitate to raise an issue on Bitbucket or ask a question
on Slack if you get stuck in this process!

## Resources

- The best resource for learning git is the git book:
    - https://git-scm.com/book/en/v2

- A couple of guides for equivalent commands in git and svn:
    - https://www.perforce.com/blog/list-of-equivalent-commands-in-git-mercurial-and-svn
    - https://backlog.com/git-tutorial/reference/commands/

- A very useful blog post on writing good commit messages:
    - https://chris.beams.io/posts/git-commit/

- Some information on pull requests:
    - https://confluence.atlassian.com/bitbucket/pull-requests-and-code-review-223220593.html

- Some guides on the workflow model:
    - https://confluence.atlassian.com/bitbucket/workflow-for-git-feature-branching-814201830.html
    - https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow
