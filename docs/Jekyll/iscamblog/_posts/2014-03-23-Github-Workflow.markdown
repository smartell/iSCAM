---
layout: post
title:  "Github and workflow"
date:   2014-03-23 06:50:11
categories: iscam overview
---

#Github

The source code and documentation for iSCAM is hosted on [Github](https://github.com/smartell/iSCAM).  The project is completely open source and we encourage others to collaborate and contribute to this project.  Unlike other stock assessment programs, iSCAM is not intended to be a completely comprehensive model that will digest all forms of fisheries data known to man kind and produce all forms of outputs required by every Fisheries Agency around the world.  Therefore, one could easily expect custom modifications to the core source code that is maintained by the core developers in the primary repository.


Git and [Github](http://github.com) are perfect tools for working in parallel and maintaining your own custom modifications to the original source code.  If  you intend to develop and maintain your own custom features, you should fork the original repository, then periodically update your forked repo with updates from the original.  

Unlike SVN, Git is a distributed version control system where each computer/server maintains a cloned version of the source code.  Git also encourages the use of branching where you always work on a *developer* branch, then release new versions in a *master* branch.

If you prefer to work with a GUI, there are several different Github Apps available out there for Windows, Mac OSX, and Linux.  I personally prefer the terminal, but would opt for the simplicity of the Github app available from Github.

#Workflow
When collaborating with a team of programmers, it is important that everyone is aware of the parts/branches and how individuals collaborate together.  Resolving code conflicts is a manual process and can be time consuming and has the potential to undo the work of someone else.

You should have access to good tools for resolving conflicts (e.g., Filemerge)
![Filemerge](iscamLogoSmall.png)

