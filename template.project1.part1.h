#ifndef TEMPLATE_H
#define TEMPLATE_H
struct seqtp   // a structure for three features of a sequence
 { char *name; // pointer to a string for keeping the name of a sequence
   char *seq;  // pointer to a char array for storing the sequence
   int  slen;  // the length of the sequence = char array length - 1,
 };            // where seq[1], seq[2], ... seq[n] are used

struct mattp    // a structure for three matrices
 {  int  rind;  // max row index
    int  cind;  // max column index
    int  **Dpt; // an int pointer to 2-dimensional array D[rind + 1][cind + 1]
    int  **Ipt; // an int pointer to 2-dimensional array I[rind + 1][cind + 1]
    int  **Spt; // an int pointer to 2-dimensional array S[rind + 1][cind + 1]
 };

struct scoretp  // a structure for scoring parameters
 { int  match;  // a positive score for a pair of identical DNA letters
   int  mismat; // a negative score for a pair of different DNA letters
   int  gopen;  // a negative score for a gap
   int  gext;   // a negative score for each letter in a gap
 };

struct algtp    // a structure for an optimal alignment  
 { int  score;  // the score of the alignment
   int  alen;   // the length of the alignment
   char *top;   // pointer to a char array for the top row of the alignment
   char *mid;   // pointer to a char array for the mid row of the alignment
   char *bot;   // pointer to a char array for the bottom row of the alignment
 }; 
void fatal(const char *msg);
FILE *ckopen(const char *name, const char *mode);
void *check_malloc(size_t amount);
void *check_realloc(void *ptr, size_t amount);
void getseq(const char *fname, struct seqtp *seqrec);
void **get2dspace(int rind, int cind);
void computemats(struct seqtp *sonept, struct seqtp *stwopt, struct scoretp *papt, struct mattp *matspt);
void producealg(struct seqtp *sonept, struct seqtp *stwopt, struct mattp *matspt, struct algtp *algopt);
void outputmat(char kind, struct seqtp *sonept, struct seqtp *stwopt, int **Mpt);
void outputalg(struct seqtp *sonept, struct seqtp *stwopt, struct scoretp *papt, struct algtp *algopt);
void freeseq(struct seqtp *seqrec);
void freealg(struct algtp *algopt);
void free2dspace(int rind, int **arr);
void freemat(struct mattp *matspt);


#endif
