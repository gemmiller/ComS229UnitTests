PROJ = sample
TESTDIR = ./tests
CODEDIR = ./code
CODE= $(CODEDIR)/*.c
GTESTFLAGS = -lgtest -lpthread

all:

test: test_runner.cc $(TESTFILES)
	@g++ -o $(PROJ)Tests test_runner.cc $(CODE) $(GTESTFLAGS)
	./$(PROJ)Tests
	@rm $(PROJ)Tests
