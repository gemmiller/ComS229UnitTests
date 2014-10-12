# README
## This project contains the unit tests for testing Com S 229 projects

### Setting up to test
   Create a students directory in the students folder
   To set up to test just copy template.project1.part1.h, mastercode.h, Makefile and the tests folder to the students dir 
   Now make sure the students code is named template.project1.part1.c rename it if you need to
   Once you're done with manual tests comment out the main function and the function prototypes from the studnents code
   Also add an include for template.project1.part1.h to the students code
   Now type make and it should build the test project
   Just type ./student_unittest and it should run the test. 

   Test Scores:
   getseq: 
   name -> 5pts
   seq -> 7pts
   len -> 3pts

   comput matrix: 
    S:
     mn entry -> 2pts
     last row -> 4 -> 2 partial
     last col -> 4 -> 2 partial
     full check -> 10 -> partial 8 -> 5

   I/D:
     mn entry -> 2pts
     last row -> 4 -> 2 partial
     last col -> 4 -> 2 partial
     full check -> 10 -> partial 8 -> 5

  produce alignment:
  General case -> 25 pts
  Insetion gap -> 10 pts
  Deletion gap -> 10 pts
  Perfect match -> 5 pts
