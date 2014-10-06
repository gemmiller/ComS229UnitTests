#include "gtest/gtest.h"

#include "./tests/sample1_unitTest.cc"

int main(int argc, char **argv) {
      ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
}
