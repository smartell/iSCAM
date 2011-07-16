#!/bin/bash

#
#  bash script for creating the iscam archive
#  TO RUN:  press command R in textmate.
#
clear
echo "Archiving iscam files"

echo "** Cleaning file system"
rm *.b*
rm *.r0*
rm *.p0*
rm *.ret*
rm *.zip
rm *.tgz.gz
rm iscam.eva iscam.mc2 iscam.hst iscam.log iscam.mcm fmin.log eigv.rpt
rm check.tmp
rm admodel.cov admodel.dep admodel.hes
echo "** Finished removing files."

# script for creating the Tape Archive file
tar -czvf iscamArchive.$(date +%Y.%m.%d).tgz.gz Examples baranov.cxx\
 iscam.tpl iscam.dat iscam.r iSCAMViewTracker.txt iScamLogo.gif\
 stats.cxx iSCAMwin.txt Riscam_1.0.tar.gz Riscam.zip read.admb.r

# script for creating a windoz zip file
zip -r iscamArchive.$(date +%Y.%m.%d).zip Examples baranov.cxx\
 iscam.tpl iscam.dat iscam.r iSCAMViewTracker.txt iScamLogo.gif\
 stats.cxx iSCAMwin.txt Riscam_1.0.tar.gz Riscam.zip read.admb.r

zip -d iscamArchive.$(date +%Y.%m.%d).zip \*.psv \*.rep \*.cor \
 \*.par \*.std \*.mcrt \*.mcst \*.mcmc
exit 0