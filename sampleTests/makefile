CC = g++
PROJ = sample
TESTDIR = ./tests
CODEDIR = ./code
CODE= $(CODEDIR)/*.c
GTESTFLAGS = -lgtest -lpthread

all: test_runner.cc $(TESTFILES)
	$(CC) -o $(PROJ)Tests test_runner.cc $(CODE) $(GTESTFLAGS)

test: test_runner.cc $(TESTFILES)
	@$(CC) -o $(PROJ)Tests test_runner.cc $(CODE) $(GTESTFLAGS)
	./$(PROJ)Tests
	@rm $(PROJ)Tests

clean:
	rm $(PROJ)Tests
