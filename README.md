
#  iSCAM 2.0 

<img src="https://raw.githubusercontent.com/smartell/iSCAM/IPHC/src/R/logo/iscamLogoSmall.png" alt="Drawing" style="width: 180px;"/>



## What is iSCAM ##
iSCAM is short for integrated statistical catch age model. The software was orginally developed by Steven Martell at the University of British Columbia, and the project was intiated in the fall of 2010.  The code was originally written for use as a stock assessment model for BC herring stocks.  Since this time the code has evolved substantially, and can now be used to conduct assessments for multiple sexs, one or more stocks, or sub-stocks, and can be spatially explicit.

## Table of Contents ##
- [What is iSCAM](#what-is-iscam)
- [Obtaining iSCAM](#obtaining-iscam)
    - [Prerequisites](#prerequisites)
- [List of developers](#list-of-developers)

## Obtaining iSCAM##
Obtaining the latest version of iSCAM should be done using the [gitHub](https://github.com/smartell/iSCAM) repository

### Prerequisites ###
* C++ compiler (preferably `clang++` that supports C++11 compliance standards)
* AD Model Builder (version 11.1 or later)
* R (version 2.15 or later)
  * PBSmodelling package (and dependencies)
  * Hmisc package (and dependencies)
  * Shiny package (and dependencies)



## List of developers##
The following people have contributed source code or ideas to the iSCAM-project:

| **Individual**   | **Organization** | _email_                     |
|:-----------------|:----------------:| ---------------------------:|
| Steven Martell   | IPHC             | stevem@iphc.int             |
| Vivian Haist     | Consultant       | haistv@shaw.ca              | 
| Jaclyn Cleary    | DFO              | Jaclyn.Cleary@dfo-mpo.gc.ca |
| Chris Grandin    | DFO              | Chris.Grandin@dfo-mpo.gc.ca |
| Robyn Forrest    | DFO              | Robyn.Forrest@dfo-mpo.gc.ca |
| James Ianelli    | NOAA/NMFS        | Jim.Ianelli@noaa.gov        |
| Dave Fournier    | Otter Research   | davef@otter-rsch.com        |
| Carl Walters     | UBC              | c.walters@fisheries.ubc.ca  |
| Rob Kronlund     | DFO              | Allen.Kronlund@dfo-mpo.gc.ca|
| Sean Cox         | SFU              | spcox@sfu.ca                |
| Nathan Taylor    | DFO              | Nathan.Taylor@dfo-mpo.gc.ca |
| Catarina Wor     | UBC/IPHC         | catarinawor@gmail.com       |
| Richard Methot   | NOAA/NMFS        | Richard.Methot@noaa.gov     |
| Matthew Supernaw | NOAA/NMFS        | matthew.supernaw@noaa.gov   |
| Chris Francis    | NIWA             | chris.francis@clear.net.nz  |
| Mark Maunder     | IATTC            | mmaunder@iattc.org          |
| Marie Etienne    | France           | mp.etienne@gmail.com        |

---





	
## Cloning the repository

### On Mac OS and Linux 
Open terminal and run the following

	cd ~
	git clone https://github.com/smartell/iSCAM.git
	cd iSCAM-project
	make

### Windows
Obtain a copy of git for windows [here](http://git-scm.com/download/win), then in a command window

	git clone https://github.com/smartell/iSCAM.git

You might useful to install [Cygwin](http://www.cygwin.com) to give you Windows machine that Linux like feel and allow you machine to use much of the functionality builtin iSCAM. Compiling iSCAM and runnig examples is carried out using GNU Makefiles.


----
# Creating new projects
## Scripts
There are a number of shell scripts for creating new projects within the examples directory and the fba directory.  The makeproject script sets up the following directory tree:


- DATA
     - Example.ctl
     - Example.dat
     - Example.mpc
     - Example.pfc
     - Example.scn
     - Makefile
     - RUN.dat
- FIGS
- MISC
- PRESENTATION
     - Example.tex
     - iScamLogo.pdf
- R
     - collectAll.R
     - saveMSEdataframe.R
- TABLES
- WRITEUP

On \*nix operating systems, at the terminal simply type:

      ./makeproject <ProjectName>

By default the makeproject script copies a set of templates that are required for iSCAM to run.  The \*.ctl, \*.dat, \*.pfc, \*.scn, \*.mpc are input files required by iscam.  The RUN.dat file is the primary file that is opened by iSCAM to determine which files are used as the contord, data, projection, scenario and manamement procedure, respectively.



# Version control

* To check out a copy of the project code, open terminal and
  go to the directory (folder) where you want to keep a clone
  (copy) of the repository on your macbook pro.  If you are 
  using Mac OSX 10.6 or later you will already have "git" 
  installed on your computer. At the command prompt type:

* git clone https://code.google.com/p/iscam-project/

* The above command will now download a copy (referred to as a 
  clone in git terminology) onto your computer in the current
  directory. You only need to run the above command once.  To 
  obtain updates of the repository, you can run the following:

* git fetch https://code.google.com/p/iscam-project/

* To gain write access to the repository contact martell.steve@gmail.com

  THE FOLLOWING IS A LIST OF USEFUL GIT COMMANDS TO USE AT THE COMMAND LINE$:
      git help                      <list git commands>
      git clone url ILSMR/          <clone repository>
      git add filename              <add new file>
      git commit -m"message"        <commit with message>
      git log                       <view commit log>
      git status                    <view changes & staging>
      git branch                    <list, create or delete branches>
  
  Making an alias for a git command, eg:
      git config --global alias.ci commit
      git config --global alias.lol "log --oneline --graph"

### Version control system resources:

Before using git, I would highly recommend spending some time 
learning how to use git.  There are many online resources and
most of them can be found at: http://git-scm.com/documentation

A video tutorial: http://www.youtube.com/watch?v=ZDR433b0HJY

Another good [cheatsheet](http://cheat.errtheblog.com/s/git/)


## Compiling code using Makefiles:
There are several [GNU Makefiles](http://www.gnu.org/software/make/manual/make.html) throughout the directories in the iSCAM project that will greatly simplify recompiling the iscam code and running multiple models in parallel using the -j option in the Makefiles.

The first makefile is in the iSCAM-project root directory (shown below).  At the command line in the iSCAM-project directory, if you type make, the makefile will create a _dist_ directory and subdirectories, compile optimized and safe versions of the iSCAM.tpl file, then copy these files into the _dist/release_ and _dist/debug_ directories.


      ## Makefile for building distribtion folder for iscam
      ## TODO add verify target to run example models.
      .PHONY: dist clean


      ifndef DISK
        DISK=dist
      endif


      dist:
        mkdir -p ${DISK}/debug
        mkdir -p ${DISK}/R
        mkdir -p ${DISK}/release
        make --directory=src/admb-code --file=linux.mak 
        make --directory=src/admb-code --file=linux.mak opt
        cp -r ./src/r-code/ ${DISK}/R/

      clean:
        make --directory=src/admb-code --file=linux.mak clean
        rm -r dist

---
## Setting up a new project in fba:

- fba is a directory (short for "full blown assessments" ) that 
  contains project folders with a specific directory structure.
  It is important to maintain this directory structure within
  the project folders because the R and shell scripts use the
  relative paths for maintaining inputs and outputs to automate
  much of the running of iSCAM and its outputs.

- To set up a new project, use the makeproject shell script (linux or MacOSX)
  when you are inside the fba directory only.
  For example:

        cd fba
    	./makeproject MyAssessment  

  It may be necessary to change the permissions if you are 
  running MacOSX, eg: type "chmod 755 makeproject" in terminal
  at the fba directory to allow execution permissions for 
  makeproject.

- The shell script will create a number of directories including _DATA_ and copy various makefiles and default control and data files that will allow you to run the default model.  If you then navigate to the _MyAssessment/DATA_ directory, you can now simply type:

        make

and this will copy the executable, and run the model inside the _DATA_ directory.





---
# Contributing

You can send pull requests via GitHub. Patches should:

1. Follow the style of the existing code.
2. One commit should do exactly one thing.
3. Commit messages should start with a summary line below 80 characters followed by a blank line, and then the reasoning/analysis for why the change was made (if appropriate).
4. Commits that fix a bug in a previous commit (which has already been merged) should start with `fixup!` and then the summary line of the commit it fixes. If you are writing your commit message in TextMate then type `fix⇥` to get the prefix and a menu allowing you to pick the summary line from one of the last 15 commits.
5. Rebase your branch against the upstream’s master. We don’t want to pull redundant merge commits.
6. **Be clear about what license applies to your patch:** The files within this repository are under the [GPL 3][] (or later) but (as the original creator) we are still allowed to create non-free derivatives. However, if patches are given to us under GPL then those cannot make it into any non-free derivatives we may later wish to create. So to make it easier for us (and avoid any legal issues) we prefer if patches are released as public domain.

There is both the [textmate-dev][] mailing list and [#textmate][] IRC channel at [freenode.net][] where this project can be discussed.

## GitHub Workflow

Developing patches should follow this workflow:

### Initial Setup

1.	Fork on GitHub (click Fork button)
2.	Clone to computer: `git clone git@github.com:«github account»/textmate.git`
3.	cd into your repo: `cd textmate`
4.	Set up remote upstream: `git remote add -f upstream git://github.com/textmate/textmate.git`

### Adding a Feature

1.	Create a branch for the new feature: `git checkout -b my_new_feature`
2.	Work on your feature, add and commit as usual

Creating a branch is not strictly necessary, but it makes it easy to delete your branch when the feature has been merged into upstream, diff your branch with the version that actually ended in upstream, and to submit pull requests for multiple features (branches).

### Pushing to GitHub

8.	Push branch to GitHub: `git push origin my_new_feature`
9.	Issue pull request: Click Pull Request button on GitHub

### Useful Commands

If a lot of changes have happened upstream you can replay your local changes on top of these, this is done with `rebase`, e.g.:

	git fetch upstream
	git rebase upstream/master

This will fetch changes and re-apply your commits on top of these.

This is generally better than merge, as it will give a clear picture of which commits are local to your branch. It will also “prune” any of your local commits if the same changes have been applied upstream.

You can use `-i` with `rebase` for an “interactive” rebase. This allows you to drop, re-arrange, merge, and reword commits, e.g.:

	git rebase -i upstream/master

_____________________________________________________________
_____________________________________________________________
         Integrated Statistical Catch Age Model (iSCAM)

                        VERSION 1.5
                Tue Jul 19 22:24:46 PDT 2011

            Created by Steven Martell on 2010-04-09 
            Copyright (c) 2010. All rights reserved.

            Last changed on:
            Source code: https://github.com/smartell/iSCAM
_____________________________________________________________
