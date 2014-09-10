# README
## This project contains the unit tests for testing Com S 229 projects

### Project rules
1. This project should only contain test code 
2. Try not to write duplicate tests 2
3. Tabs should be 2 spaces 3

### Setting up gtest on Linux
**This setup has only been tested on Ubuntu**
**It may need to be modified for other distibutions**

1. download gtest using wget http://googletest.googlecode.com/files/gtest-1.7.0.zip

2. unzip the file 

```
	$ unzip gtest-1.7.0.zip
```

3. configure gtest 

```
	$ cd gtest-1.7.0
	$ ./configure
	$ make
```

4. add gtest to your include path 

```
	$ sudo cp -a include/gtest /usr/include/
	$ sudo cp -a lib/.libs/*	/usr/lib/
```

5. update the linker and check to see if libs are added 

```
	$ sudo ldconfig -v | grep gtest
```

if the out put looks like:

```
libgtest.so -> libgtest.so.0.0.0
libgtest_main.so.0 ->libgtest_main.so.0.0.0
```
