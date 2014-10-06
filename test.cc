// Include gtest for the tests
 #include <gtest/gtest.h>

 //include the files to test


 /**
  * By having the TEST macro gtest will be able to find the tests
  * with out them having to be registered
  */
TEST(IsPositiveTest, PositiveNum){
     EXPECT_EQ (1,1);
}

TEST(IsPositiveTest, ZeroAndNegativeNum){
     EXPECT_EQ (0, -1);
     EXPECT_EQ (0, 0);
}

TEST(IsNegative, NegativeNum){
     EXPECT_EQ (1, -1);
}

TEST(IsNegative, ZeroAndPositiveNum){
     EXPECT_EQ (0, 1);
     EXPECT_EQ (0, 0);
}
