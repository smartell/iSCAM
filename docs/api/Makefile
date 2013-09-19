.PHONY: default clean


default: manDoxyfile Doxyfile
	mkdir -vp dox/html
	doxygen manDoxyfile
	doxygen Doxyfile


clean:
	rm -rf dox

# REV=`svnversion -n ../..`
# VER=`cat ../../VERSION`

	#(cat Doxyfile ; echo "PROJECT_NUMBER = "$(VER)"."$(REV); echo "FILE_VERSION_FILTER  = ./getID.sh"; ) | doxygen -

# default: Doxyfile ./dox-share/mainpage.txt 
# 	mkdir -vp dox/html
# #	Insert some linux-specific deltails into configuration file
# #	See Doxygen FAQ "Can I configure doxygen from the command line?"
# 	(cat Doxyfile ; echo "PROJECT_NUMBER = "$(VER)"."$(REV); echo "FILE_VERSION_FILTER  = ./getID.sh"; ) | doxygen -

# clean:
# 	rm -rf dox
