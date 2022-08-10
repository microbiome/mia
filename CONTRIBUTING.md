# Contributing to this R package

## Etiquette and conduct

Contributors and maintainers are expected to abide by the
[Bioconductor Code of Conduct](https://bioconductor.org/about/code-of-conduct/)

This instructions have been adapted from 

- [(Bioconductor contributor guidelines)](https://github.com/Bioconductor/Contributions/blob/master/CONTRIBUTING.md)
- [(BioJulia contributor guidelines)](https://github.com/BioJulia/Contributing)


## How can I contribute?


### Reporting Bugs

Following these guidelines will help us to better understand your
report, replicate the issues, and track down the source problems.

#### Before creating a bug report:

Please do the following:

1. Check the GitHub issues for this package

2. If you find an open issue that describes the same problem, kindly add a comment to let
  others know that you are experiencing the same issue.

3. If no **open** issue already exists for your problem 
   then kindly create a new issue.

   > **Note:** If you find a **Closed** issue that seems like it is the same thing
   > that you're experiencing, open a new issue and include a link to the original
   > issue in the body of the new one.

#### How to create a (good) new bug report:

Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/).

When you are creating a bug report, please do the following:

1. **Explain the problem**

   - *Use a clear and descriptive title* for the issue to identify the problem.
   - *Describe the exact steps which reproduce the problem* in as many details as possible.
   - *Provide a specific example*. (Includes links to pastebin, gists and so on.)
       If you're providing snippets in the issue, use
       [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).

   - *Describe the behaviour you observed after following the steps*
     - Point out what exactly is the problem with that behaviour.
     - *Explain which behaviour you expected to see instead and why.*
     - *OPTIONALLY: Include screenshots and animated GIFs* which show you
       following the described steps.

2. **Provide additional context for the problem (some of these may not always apply)**

   - *Did the problem start happening recently* (e.g. after updating to a new version)?
     - If the problem started recently, *can you reproduce the problem in older versions?*
     - Do you know the most recent package version in which the problem doesn't happen?

   - *Can you reliably reproduce the issue?* If not...
     - Provide details about how often the problem happens.
     - Provide details about under which conditions it normally happens.

   - Is the problem is related to *working with files*? If so....
     - Does the problem happen for all files and projects or only some?
     - Does the problem happen only when working with local or remote files?
     - Does the problem happen for files of a specific type, size, or encoding?
     - Is there anything else special about the files you are using?

3. **Include details about your configuration and environment**

- Packages and versions you are using? (see R command sessionInfo())

- Name and version of the OS you're using?


### Suggesting enhancements

Enhancements include new features as well as minor improvements to
existing functionality. Following these suggestions will help
maintainers and the community understand your suggestion and find
related suggestions.

#### Before Submitting An Enhancement Proposal

1. **Perform a cursory issue search** to see if the enhancement has already been suggested.
2. If it has not, open a new issue as per the guidance below.
3. If it has...
  1. Add a comment to the existing issue instead of opening a new one.
  2. If it was closed, take the time to understand why this was so (it's ok to
     ask! :) ), and consider whether anything has changed that makes the reason
     outdated. If you can think of a convincing reason to reconsider the
     enhancement, feel free to open a new issue as per the guidance below.

#### How to submit a (good) new enhancement proposal

Enhancement proposals are tracked as
[GitHub issues](https://guides.github.com/features/issues/).

1. **Explain the enhancement**
   - *Use a clear and descriptive title* for the issue to identify the suggestion.
   - *Provide a step-by-step description of the suggested enhancement* in as many details as possible.
   - *Provide specific examples to demonstrate the steps*.
     Include copy/pasteable snippets which you use in those examples, as
     [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).

   - If you want to change current behaviour...
     - Describe the *current* behaviour.
     - *Explain which behaviour you expected* to see instead and *why*.
     - *Will the proposed change alter APIs or existing exposed methods/types?*
       If so, this may cause dependency issues and breakages, so the maintainer
       will need to consider this when versioning the next release.
   - *OPTIONALLY: Include screenshots and animated GIFs*.

2. **Provide additional context for the enhancement**

   - *Explain why this enhancement would be useful* to users and
     isn't something that can or should be implemented as a separate package.

   - *Do you know of other projects where this enhancement exists?*

3. **Include details about your configuration and environment**

   - Package versions (see R sessionInfo()).

   - Name and version of the OS

### Making Pull Requests

All packages can be developed locally, and contributions suggested via pull requests.

Before you start working on code, it is often a good idea to open an enhancement
[suggestion](#suggest-an-enhancement)

Hit the 'Fork' button on the repositories page to create a forked copy
of the target package for your own Github account. This will ensure
your work and experiments won't hinder other users of the released and
stable package.

From there you can clone your fork of the package and work on it on your
machine using git.
Here's an example of cloning, assuming you already forked `mia`:

```sh
git clone https://github.com/<YOUR_GITHUB_USERNAME_HERE>/mia.git
```

Git will download or "clone" your fork and put it in your local folder.

It is beyond the scope of this document to describe good git and github use in
more specific detail, as the folks at Git and GitHub have already done that wonderfully
on their own sites.

If you have additional questions,
see the contact details at [microbiome.github.io](microbiome.github.io).

#### How to make (good) code contributions and new Pull-Requests

1. **In your code changes**

   - **Branch properly!**
     - If you are making a bug-fix, then you need to checkout your bug-fix branch
       from the last release tag.
     - If you are making a feature addition or other enhancement, checkout your
       branch from master.
     - See [here](#a-suggested-branching-model) for more information (or ask a package maintainer :smile:).

   - Please comment liberally for complex pieces of internal code to facilitate comprehension.

2. **In your pull request**

   - *Describe* the changes in the pull request

   - Provide a *clear, simple, descriptive title*.

   - Do not include issue numbers in the PR title.

   - If you have implemented *new features* or behaviour
     - *Provide a description of the addition* in as many details as possible.
     - *Provide justification of the addition*.
     - *Provide a runnable example of use of your addition*. This lets reviewers
       and others try out the feature before it is merged or makes it's way to release.

   - If you have *changed current behaviour*...
     - *Describe the behaviour prior to you changes*
     - *Describe the behaviour after your changes* and justify why you have made the changes.
     - *Does your change alter APIs or existing exposed methods/types?*
       If so, this may cause dependency issues and breakages, so the maintainer
       will need to consider this when versioning the next release.
     - If you are implementing changes that are intended to increase performance, you
       should provide the results of a simple performance benchmark exercise
       demonstrating the improvement. Especially if the changes make code less legible.

#### Reviews and merging

You can open a pull request early on and push changes to it until it is ready,
or you can do all your editing locally and make a pull request only when it is
finished - it is up to you.

When your pull request is ready on Github,
mention one of the maintainers of the repo in a comment e.g. `@antagomir`
and ask them to review it.
You can also use Github's review feature.
They will review the code and documentation in the pull request,
and will assess it.

Your pull request will be accepted and merged if:

1. The dedicated package maintainers approve the pull request for merging.
2. The automated build system confirms that all unit tests pass without any issues.

It may also be that the reviewers or package maintainers will want to
you to make changes to your pull request before they will merge it.
Take the time to understand why any such request has been made, and
freely discuss it with the reviewers.

## Styleguides

### Git Commit messages

* Use the present tense ("Add feature" not "Added feature").
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...").
* Limit the first line to 72 characters or less.
* Reference issues and pull requests liberally after the first line.
* Consider starting the commit message with an applicable emoji:
    * :art: `:art:` when improving the format/structure of the code
    * :racehorse: `:racehorse:` when improving performance
    * :memo: `:memo:` when writing docs
    * :penguin: `:penguin:` when fixing something on Linux
    * :apple: `:apple:` when fixing something on macOS
    * :checkered_flag: `:checkered_flag:` when fixing something on Windows
    * :bug: `:bug:` when fixing a bug
    * :fire: `:fire:` when removing code or files
    * :green_heart: `:green_heart:` when fixing the CI build
    * :white_check_mark: `:white_check_mark:` when adding tests
    * :arrow_up: `:arrow_up:` when upgrading dependencies
    * :arrow_down: `:arrow_down:` when downgrading dependencies
    * :exclamation: `:exclamation:` when removing warnings or depreciations

### Additional style suggestions

- Indent with 4 spaces.

- Separate logical blocks of code with one blank line. Although it is common
  and acceptable for short single-line functions to be defined together on
  consecutive lines with no blank lines between them.

- Function names, apart from constructors, are all camelCase.

- Generally try to keep lines below 80-columns, unless splitting a long line
  onto multiple lines makes it harder to read.

