# README
## This project contains the unit tests for testing Com S 229 projects

### Project rules
-This project should only contain test code
-Try not to write duplicate tests
-Tabs should be 2 spaces

### Setting up gtest on Linux
**This setup has only been tested on Ubuntu**
**It may need to be modified for other distibutions**
-download gtest using wget http://googletest.googlecode.com/files/gtest-1.7.0.zip
-unzip the file
'''
	$ unzip gtest-1.7.0.zip
'''
- configure gtest
'''
	$ cd gtest-1.7.0
	$ ./configure
	$ make
'''
-add gtest to your include path
'''
	$ sudo cp -a include/gtest /usr/include/
	$ sudo cp -a lib/.libs/*	/usr/lib/
'''

-update the linker and check to see if libs are added
'''
	$ sudo ldconfig -v | grep gtest
'''
if the out put looks like:

'''
libgtest.so -> libgtest.so.0.0.0
libgtest_main.so.0 ->libgtest_main.so.0.0.0
'''
