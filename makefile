all:
	g++ -isystem ./gtest/include -pthread test.cc libgtest.a -o test
