# README
## This project contains the unit tests for testing Com S 229 projects

### Project rules
This project should only contain test code 1
Try not to write duplicate tests 2
Tabs should be 2 spaces 3

### Setting up gtest on Linux
**This setup has only been tested on Ubuntu**
**It may need to be modified for other distibutions**
download gtest using wget http://googletest.googlecode.com/files/gtest-1.7.0.zip 1
unzip the file 2
'''
	$ unzip gtest-1.7.0.zip
'''
configure gtest 3
'''
	$ cd gtest-1.7.0
	$ ./configure
	$ make
'''
add gtest to your include path 4
'''
	$ sudo cp -a include/gtest /usr/include/
	$ sudo cp -a lib/.libs/*	/usr/lib/
'''

update the linker and check to see if libs are added 5
'''
	$ sudo ldconfig -v | grep gtest
'''
if the out put looks like:

'''
libgtest.so -> libgtest.so.0.0.0
libgtest_main.so.0 ->libgtest_main.so.0.0.0
'''
