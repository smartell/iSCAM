#!/bin/bash

#script for zip file for Fisheries Research Submission.
clear

rm *.pdf
zip -r CAPAM_MS.zip \
../SelexCoverLetter.pdf \
Selex2.tex \
Selex2.bbl \
*.eps
#Selex.pdf \
#./Introduction/Intro.tex \
#./Methods/Methods.tex \
#./Results/Results.tex \
#./Discussion/Discussion.tex \