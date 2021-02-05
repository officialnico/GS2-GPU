title: Release instructions

1. **Take stock of the open issues and PRs** in the GS2, Makefiles and utils repos. Decide which will be resolved/merged in this release, and which will be deferred. Also use this as an opportunity to tidy up the issues/PRs. Mark finished issues as resolved, and give open issues a milestone.
2. **Make release candidate for GS2, utils and Makefiles** in branches called `<version_number>-RC`. For feature releases, these should branch from `GS2 and  utils' `next` branch, while bugfix releases should branch from their ``master` branches. It should always branch from Makefiles' `master` as there is no `next` on Makefiles.
	- It may be convenient to make a branch
      `<version_number>-RC+latest-externals` that contains all the updates to
      utils which did not have an accompanying GS2 PR. Therefore its only change
      should be the `utils` and `Makefiles` git hashes. Make a PR for this
      branch into the release candidate but do not merge yet, as it may need to
      be updated. See [this
      example](https://bitbucket.org/gyrokinetics/gs2/pull-requests/250/).
3. **Make an issue on the GS2 repo for the release** called `<vesion_number>-RC` (see [this example](https://bitbucket.org/gyrokinetics/gs2/issues/65/milestone-802)). List and link to the PRs to be included, and mark those that have been merged with a tick. The issue both advertises the new release and allows discussion of proposed PRs. The issue description can be updated to reflect changing plans.  
It probably isn't necessary to make an issue for the release in `utils` or `Makefiles`, but it might be appropriate if there are many changes.
4. **Edit the PRs for this release to target the RC branch**. This helps clarify the work that needs to be done just from looking at the PR page. It also helps check if proposed changes are possible/reasonable. Note: admins can edit other users' PRs - the option appears in the "..." menu.
5. **Reserve a release DOI on Zenodo.** Go to [Zenodo](https://zenodo.org/record/2551066),
   click "New Version", enter the name as `GS2 v<version_number>`, and update
   the version number field. Make sure "Reserve DOI" is selected, and copy the
   DOI.  Click save. **DO NOT CLICK PUBLISH**
6. **Update `externals/utils` and `Makefiles`.**:
    1. Merge the release's PRs into the release candidate.
    2. Update `CHANGELOG.md`. The changelog file is an automatically-generated
      list of PRs, but **contains manual edits that must be kept**. To generate
      the changelog, `git tag` the current head with the version number, then
      run the generator script:  
       ```
         $ <gs2_topdir>/scripts/change-log-builder.sh > CHANGELOG-new.md
       ```
       Then copy the new commits from `CHANGELOG-new.md` into `CHANGELOG.md`,
       so that manual tweaks to CHANGELOG.md are kept. Delete the tag - we will 
       need it later for the final merge commit. 
    3.  Make a placeholder PR called `[WIP] <version_number>` from the RC branch into `master` (for both utils and makefiles`). 
    4.  Assess open PRs: merge to the release candidate or defer to the next release.  
    5.  After merging the `utils` RC branch into `utils` master, make a `gs2` PR
       updating `externals/utils` and `Makefiles` to their latest master commits
7. **Repeat step 6 for the main GS2 repo**.
8. **Update the list of authors in `CITATION.cff` and on Zenodo.** GS2 has been
   developed for over 20 years, and there will have been significant
   contributions of which you are unaware. Do not remove authors. The best we
   can achieve is to reflect the current state of authorship. The rule used for
   listing authors (from 8.0.1 onwards) is to use the following order: the GS2
   Governing Committee alphabetically (Michael Barnes, David Dickinson, Bill
   Dorland, Peter Hill, Joseph Parker, Colin Roach); recent contributors
   alphabetically (here "recent" is given some leeway - it probably means
   "committed to the current major version" - but should be used to list other
   active developers); all other contributors. To obtain a list of commit
   authors, use the commands:  
   ```
     $ git log <previous version tag>..HEAD --format='%aN' | sort | uniq
   ```
9. **Update Zenodo DOI in README and CITATION.cff.** Update the DOI in the
   Zenodo badge at the top of the README. Also update the DOI, release date and
   version number in CITATION.cff.
10. **Make a PR for the release candidate into `master`.** List main fixes and
   features. Highlight changes in defaults, or manual changes that users might
   need to make. Give this releases' DOI. Invite users to try the branch and
   raise issues. See [this
   example](https://bitbucket.org/gyrokinetics/gs2/pull-requests/135).
11. **Wait.** Allow some time for users to try the release candidate and for its
    manifold failings to be exposed. Fix these with PRs to the release candidate
    branch.
12. **Merge.** After merge, make sure to tag the merge commit with the version
    number. It is **necessary to use an annotated tag** (use -a), otherwise
   `git describe` does not work.
13. **Finish Zenodo release.** Download the tarball from the [Bitbucket downloads page](https://bitbucket.org/gyrokinetics/gs2/downloads/?tab=tags),
    along with the corresponding tarballs for utils` and `Makefiles`. Upload all
    three to Zenodo, and then publish the Zenodo draft.
14. **Merge master into next.** To do this, we need to make a new branch
    containing the release's changes so that we can make a pull request into
    `next`:  
    ```
      $ git checkout next
    ```  
    ```
      $ git checkout -b admin/update_next_with_8.0.6
    ```  
    ```
      $ git merge master
    ```  
    Then fix conflicts and make a PR to `next`.
15. **Make a new RC branch for the next release.**
