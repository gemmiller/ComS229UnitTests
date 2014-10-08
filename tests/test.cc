// Include gtest for the tests
 #include <gtest/gtest.h>
#include "../template.project1.part1.h"
#include <string.h>
 //include the files to test


 /**
  * By having the TEST macro gtest will be able to find the tests
  * with out them having to be registered
  */
TEST(GETSEQ, Check_Name_5PTS){
    // Arrange
   struct seqtp seq; 

   // Act
   getseq("tests/A1.txt", &seq);

   // Assert
   ASSERT_STREQ ("A",seq.name);
}

TEST(GETSEQ, Check_Seq_7PTS){
    // Arrange
    struct seqtp seq;
    char exp[12] = "AGTAACGACCT";

    // Act
    getseq("tests/A1.txt",&seq);
    if(strlen(seq.seq)==12)
        seq.seq = seq.seq+1;
    if(strlen(seq.seq)==13){
        seq.seq[12]=0;
        seq.seq = seq.seq+1;
    }

    // Assert
    ASSERT_STREQ("AGTAACGACCT",seq.seq);
}

TEST(GETSEQ, Check_Length_3PTS){
    // Arrange
    struct seqtp seq;

    // Act
    getseq("tests/A1.txt",&seq);

    // Assert
    ASSERT_EQ(11,seq.slen);
}

TEST(COMPUTE_MAT, Check_S_MN_2PT){
   // Arrange
   struct seqtp seq1= {.name = "A",.seq = "AGTAACGACCT",.slen=11};
   struct seqtp seq2 = {.name = "A",.seq = "AGTAACGACCT",.slen=11};
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats;

   // Act
   computemats(&seq1,&seq2,&scores,&mats);

   // Assert
   ASSERT_EQ(0,mats.Spt[11][11]);
}

