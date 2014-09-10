// Include gtest for the tests
#include "gtest/gtest.h"

//include the files to test

#include "../code/sample1.h"

/**
 * By having the TEST macro gtest will be able to find the tests
 * with out them having to be registered
 */
TEST(IsPositiveTest, PositiveNum){
	EXPECT_EQ (1,isPositive(10));
}

TEST(IsPositiveTest, ZeroAndNegativeNum){
	EXPECT_EQ (0, isPositive(-1));
	EXPECT_EQ (0, isPositive(0));
}

TEST(IsNegative, NegativeNum){
	EXPECT_EQ (1, isNegative(-1));
}

TEST(IsNegative, ZeroAndPositiveNum){
	EXPECT_EQ (0, isNegative(1));
	EXPECT_EQ (0, isNegative(0));
}
