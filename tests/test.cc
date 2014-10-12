// Include gtest for the tests
 #include <gtest/gtest.h>
#include "../template.project1.part1.h"
#include <string.h>
#include "../masterCode.h"
 //include the files to test
 
void m_getseq(const char *fname, struct seqtp *seqrec)
// Input parameter fname is a pointer to a string for storing the name of a DNA sequence file.
// Output parameter seqrec is a pointer to a structure whose members are used to
// save the name, length and DNA letters of the sequence in the file.
// A proper amount of memory needs to be allocated for each pointer member.

// This function opens the file whose name is pointed to by fname,
// reads the name, length and DNA letters of the sequence in the file,
// dynamically allocates a proper amount of memory for each pointer member
// in the structure pointed to by seqrec, and saves the data into the member variables.
// The function prints out a proper error message if fname is NULL,
// no file exists with the given name, or the sequence/its name is empty.
{
  int  M;           // Sequence length
  char *A;          // for sequence
  int  symbol;      // The next character
  FILE *Ap;
  int  i;
  int  size;
  char *dhead1;

	/* determine the sequence length */
	if ( fname == NULL )
	  m_fatal("getseq: file name pointer is NULL\n");
	if ( seqrec == NULL )
	  m_fatal("getseq: structure seqtp pointer is NULL\n");
	Ap = m_ckopen(fname, "r");
	if ( (symbol = getc(Ap) ) != '>' )
	  m_fatal("Sequence one must be in FASTA format, starting with >\n");
	for (size = 2; ( symbol = getc(Ap)) != EOF ; size++ )
	  if ( symbol == '\n' )
	   break;
	for (M = 0; ( symbol = getc(Ap)) != EOF ; )
	   if ( symbol != '\n' )
	      ++M;

	/* allocate space for A */
	A = ( char * ) m_check_malloc( (M + 1) * sizeof(char));
	dhead1 = ( char * ) m_check_malloc( size * sizeof(char));

	/* read the sequence into A */
	rewind(Ap);
	symbol = getc(Ap);
	for (i = 0; ( symbol = getc(Ap)) != EOF && symbol != '\n' ; )
	   dhead1[i++] = symbol;
	dhead1[i] = '\0';
	for (M = 0; ( symbol = getc(Ap)) != EOF ; )
	   if ( symbol != '\n' )
		A[++M] = symbol;
        fclose(Ap);
	seqrec->name = dhead1;
	seqrec->seq = A;
	seqrec->slen = M;
}

void m_freeseq(struct seqtp *seqrec)
// Frees each dynamically allocated memory block in the seqtp structure.
{
  free(seqrec->name);
  free(seqrec->seq);
  seqrec->name = seqrec->seq = NULL;
}

void **m_get2dspace(int rind, int cind)
// Allocates a 2-dimensional int array of (rind + 1) by (cind + 1)
// and returns its memory address.
// Prints out an error message if rind or cind is less than 1.
{
  void  **a2pt;
  int  j;
  if ( rind < 0 || cind < 0 )
    m_fatal("get2dspace: negative parameter\n");
  a2pt = (void**)m_check_malloc( (rind+1) * sizeof(int * ) );
  for ( j = 0; j <= rind; j++ ) 
      a2pt[j] = (int *) m_check_malloc( (cind+1) * sizeof(int) );
  return a2pt;
}

void m_free2dspace(int rind, int **arr)
// Deallocates a 2-dimensional int array of first dimension 
// (rind + 1) by releasing the memory for arr[j] for each j
// and then releasing the memory for arr.
// Prints out an error message if rind is less than 0.
{
  int  j;
  if ( rind < 0 )
    m_fatal("free2dspace: negative parameter\n");
  for ( j = 0; j <= rind; j++ )
      free(arr[j]);
  free(arr);
}

void m_freemat(struct mattp *matspt)
// Frees each dynamically allocated memory block in the mattp structure.
{
  m_free2dspace(matspt->rind, matspt->Spt);
  m_free2dspace(matspt->rind, matspt->Dpt);
  m_free2dspace(matspt->rind, matspt->Ipt);
  matspt->Spt = matspt->Dpt = matspt->Ipt = NULL;
}

void m_computemats(struct seqtp *sonept, struct seqtp *stwopt, struct scoretp *papt, struct mattp *matspt)
// Input parameter sonept is a pointer to a structure with the features of one DNA sequence.
// Input parameter stwopt is for another DNA sequence.
// Input parameter papt is for the scoring parameters.
// Output parameter matspt is a pointer to a structure with three matrices.

// The function sets the rind member of the structure pointed to
// by matspt to the slen member of the structure pointed to by sonept,
// and cind to the slen member by stwopt. Then it allocates space for each
// matrix in the mattp structure by calling get2dspace(). Next it computes,
// in reverse order, the three matrices by the dynamic programming algorithm
// for aligning the two sequences.
{
  char *A = sonept->seq;  // pointer to sequence A
  int   M = sonept->slen; // the length of A
  char *B = stwopt->seq;  // pointer to sequence B
  int   N = stwopt->slen; // the length of B
  int   i, j;
  int   *spt;  // pointer to the current row of S
  int   *dpt;  // pointer to the current row of D
  int   *ipt;  // pointer to the current row of I
  int   *sprpt; // pointer to the previous row of S
  int   *dprpt; // pointer to the previous row of D
unsigned char  a, b;   // bases at the current row and column
  int   q, r;   // gap open and extension penalties;
  int   match;  // match score
  int   mismat; // mismatch score
  int  svalue, dvalue, ivalue; // scores in S, D and I.
  int  tmp;  // temporary
  int  spres; // S(i+1, j+1)
  int  sub; // sigma(a_{i+1}, b_{j+1})

  if ( M == 0 && N == 0 )
     m_fatal("computemats: both sequences are of length 0\n");

  if ( matspt == NULL )
     m_fatal("computemats: matspt is NULL.\n");

  if ( papt == NULL )
     m_fatal("computemats: spt is NULL.\n");
  q = papt->gopen;
  r = papt->gext;
  match = papt->match;
  mismat = papt->mismat;
  int  qr = q + r;

  matspt->rind = M;
  matspt->cind = N;
  // allocates space for each matrix.
  matspt->Spt = (int **) m_get2dspace(M, N);
  matspt->Dpt = (int **) m_get2dspace(M, N);
  matspt->Ipt = (int **) m_get2dspace(M, N);

  spt = matspt->Spt[M]; // row M is the current row
  dpt = matspt->Dpt[M];
  ipt = matspt->Ipt[M];
  spt[N] = svalue = 0; // S(m, n) = 0
  ipt[N] = dpt[N] = dvalue = ivalue = svalue - q; // D(m, n) = I(m, n) = S(m, n) - q
  for ( j = N-1; j >= 0; j-- )
   { ivalue -= r; // I(m, j) = I(m, j+1) - r
     spt[j] = ipt[j] = svalue = ivalue; // S(m, j) = I(m, j)
     dpt[j] = svalue - q; // D(m, j) = S(m, j) - q
   } // computing row M

   for ( i = M-1; i >= 0; i-- )
    {
      sprpt = spt; // current row becomes the previous row for S
      spt = matspt->Spt[i]; // a new current row for S
      dprpt = dpt; // previous row for D
      dpt = matspt->Dpt[i]; // current row for D
      ipt = matspt->Ipt[i]; // the current row for I
      spt[N] = dpt[N] = svalue = dprpt[N] - r; // S(i, n) = D(i, n) = D(i+1, n) - r
      ivalue = ipt[N] = svalue - q; // I(i, n) = S(i, n) - q
      a = A[i+1]; // a_{i+1}
      for ( j = N-1; j >= 0; j-- )
       { // ivalue for I(i, j+1); svalue for S(i, j+1)
         if ( (ivalue = ivalue - r) < (tmp = svalue - qr) )
	   ivalue = tmp; // I(i, j) = max{I(i, j+1) - r, S(i, j+1) - q - r}
         ipt[j] = ivalue; // saves it in I(i, j)

	 if ( (dvalue = dprpt[j] - r) < (tmp = sprpt[j] - qr) )
	   dvalue = tmp; // D(i, j) = max{D(i+1, j) - r, S(i+1, j) -q - r}
         dpt[j] = dvalue; // saves it in D(i, j)
      
	 spres = sprpt[j+1]; // S(i+1, j+1)
	 b = B[j+1];
	 // sigma(a_{i+1}, b_{j+1})
	 sub = (tolower(a) == tolower(b)) ? match : mismat;
	 svalue = sub + spres; // sigma(a_{i+1}, b_{j+1}) + S(i+1, j+1)
	 if ( svalue < dvalue ) // less than D(i, j)
	    svalue = dvalue;
	 if ( svalue < ivalue ) // less than I(i, j)
	    svalue = ivalue;
	 spt[j] = svalue; // saves it in S(i, j)
       } // for each column j
    } // for each row i
}

void m_outputmat(char kind, struct seqtp *sonept, struct seqtp *stwopt, int **Mpt)
// Input parameter kind denotes matrix type: 'D', 'I', or 'S'.
// Input parameter sonept is a pointer to a structure with a sequence of length M
// Input parameter sonept is a pointer to a structure with a sequence of length N
// Input parameter Mpt is a pointer to a 2-dimensional array (matrix)
// of (M + 1) by (N + 1).

// The function reports the matrix type and each value in the matrix
// on the stdout. The values are reported first in order of row index i 
// and then in order of column index j in the form of (i, j, Mpt[i][j])
// so that the row and column coordinates of a value are easily recognized.
// This function is for your own use to check on each matrix,
// so any format is OK.
{
  int i, j;
  int *row;
  int rind, cind;
  char *A, *B;

  A = sonept->seq;
  rind = sonept->slen;
  B = stwopt->seq;
  cind = stwopt->slen;
  printf("Matrix type: %c\n\n", kind);
  printf(" B         ");
  for ( j = 0; j <= cind; j++ )
    if ( j )
      printf("  %c ", B[j]);
    else
      printf("    ");
  printf("\n");
  printf("   Column  ");
  for ( j = 0; j <= cind; j++ )
      printf(" %2d ", j);
  printf("\n");
  printf("A\n");
  for ( i = 0; i <= rind; i++ )
   {  
      row = Mpt[i];
      if ( i )
         printf("%c row %2d  ",A[i], i);
      else
         printf("  row %2d  ", i);
      for ( j = 0; j <= cind; j++ )
        printf("%4d", row[j]);
      printf("\n\n");
   }
  printf("\n");
}

void m_producealg(struct seqtp *sonept, struct seqtp *stwopt, struct mattp *matspt, struct algtp *algopt)
// Input parameter sonept is a pointer to a structure with the features of one DNA sequence.
// Input parameter stwopt is for another DNA sequence.
// Input parameter matspt is a pointer to a structure with three matrices.
// Output parameter algopt is a pointer to a structure for an optimal alignment.

// The function allocates memory for each row in the algtp structure,
// produces an optimal alignment by tracing through the matrices in the mattp structure,
// and saves the alignment along with it score and length in the algtp structure.
// The alignment is represented by using three rows of characters (three char arrays).

// Note that the length of an optimal alignment cannot exceed the sum of
// the lengths of the two DNA sequences.
{
  char *A;  // pointer to sequence A
  int   M;  // the length of A
  char *B;  // pointer to sequence B
  int   N;  // the length of B
  int   i, j;
  int   **Spt;  // pointer to the S array
  int   **Dpt;  // pointer to the D array
  int   **Ipt;  // pointer to the I array
unsigned char  a, b;   // bases at the current row and column
  int   q, r;   // gap open and extension penalties;
  char *top, *mid, *bot; // pointers to char arrays for alignment
  int  pos;	// the current index for the arrays
  char mat;	// matrix type
  int  qr;      // q + r
  int  alen; // max length of the three char arrays in a structure pointed to by algopt.

  if ( sonept == NULL )
     m_fatal("producealg: sonept is NULL.\n");
  if ( stwopt == NULL )
     m_fatal("producealg: stwopt is NULL.\n");
  if ( matspt == NULL )
     m_fatal("producealg: matspt is NULL.\n");
  if ( algopt == NULL )
     m_fatal("producealg: algopt is NULL.\n");

  A = sonept->seq;
  M = sonept->slen;
  B = stwopt->seq;
  N = stwopt->slen;
  if ( M == 0 && N == 0 )
     m_fatal("producealg: both sequences are of length 0\n");

  alen = M + N;
  algopt->top = top = ( char * ) m_check_malloc( (alen + 2) * sizeof(char));
  algopt->mid = mid = ( char * ) m_check_malloc( (alen + 2) * sizeof(char));
  algopt->bot = bot = ( char * ) m_check_malloc( (alen + 2) * sizeof(char));

  pos = i = j = 0;
  mat = 'S';
  Spt = matspt->Spt;
  Dpt = matspt->Dpt;
  Ipt = matspt->Ipt;
  qr = M > 0? (- Dpt[M-1][N]) : ( - Ipt[M][N-1]);
  if ( qr < 0 )
     m_fatal("producealg: wrong q + r value\n");
  while ( i <= M && j <= N )
   {
     if ( mat == 'S' )
      {
        if ( i == M && j == N )
	   break;
        if ( j == N || Spt[i][j] == Dpt[i][j] )
	 { mat = 'D'; continue; }
        if ( i == M || Spt[i][j] == Ipt[i][j] )
	 { mat = 'I'; continue; }
        top[pos] = a = A[i+1];
        bot[pos] = b = B[j+1];
        mid[pos] = (tolower(a) == tolower(b)) ? '|' : ' ';
	if ( ++pos > alen )
          m_fatal("producealg: The alignment is too long.\n");
        i++; j++;
	continue;
      }

     if ( mat == 'D' )
      {
        top[pos] = A[i+1];
        bot[pos] = ' ';
        mid[pos] = '-';
	if ( ++pos > alen )
          m_fatal("producealg: The alignment is too long.\n");
        if ( i == M-1 || Dpt[i][j] == Spt[i+1][j] - qr )
	     mat = 'S';
        i++;
	continue;
      }

     if ( mat == 'I' )
      {
        top[pos] = ' ';
        bot[pos] = B[j+1];
        mid[pos] = '-';
	if ( ++pos > alen )
          m_fatal("producealg: The alignment is too long.\n");
        if ( j == N-1 || Ipt[i][j] == Spt[i][j+1] - qr )
             mat = 'S';
        j++;
	continue;
      }
   }
  
  if ( i < M || j < N )
      m_fatal("producealg: An error in traceback.\n");
  algopt->alen = pos;
  algopt->score = Spt[0][0];
}

void m_freealg(struct algtp *algopt)
// Frees each dynamically allocated memory block in the algtp structure.
{
  free(algopt->top);
  free(algopt->mid);
  free(algopt->bot);
  algopt->top = algopt->mid = algopt->bot = NULL;
}

void m_outputalg(struct seqtp *sonept, struct seqtp *stwopt, struct scoretp *papt, struct algtp *algopt)
// Input parameter sonept is a pointer to a structure with the features of one DNA sequence.
// Input parameter stwopt is for another DNA sequence.
// Input parameter matspt is a pointer to a structure with three matrices.
// Input parameter algopt is a pointer to a structure for an optimal alignment.

// The function reports, on the stdout, a summary of sequence and alignment information
// and the alignment. The summary includes the name and length of each sequence,
// the values of the four scoring parameters, and the score and length of the alignment.
// The alignment is reported in sections of 70 characters, with each section consisting
// of three rows. The sequence positions of the first DNA letters in each section are
// reported in the left margin of 10 spaces.

// Students are welcome to use the functions below.
// Remember to include each function prototype before main().
{
   int pos;      // alignment position
   int i, j, k;  // sequence positions
   char *top, *mid, *bot; // pointers to char arrays for alignment
   char tline[100], mline[100], bline[100]; // a section of alignment
   int  linelen = 70; // length of the section
   int  limit;   // limit on alignment position
   char *tpt, *mpt, *bpt; // pointers to the three rows in an alignment section
   char a, b;    // the current sequence letters
   int  ast, bst; // the positions of the leftmost letters in the sequences.
   char sline[100]; // line position marks
   char *spt;     // pointer to a sline position
   int  s;       // a sline position

  if ( sonept == NULL )
     m_fatal("producealg: sonept is NULL.\n");
  if ( stwopt == NULL )
     m_fatal("producealg: stwopt is NULL.\n");
  if ( papt == NULL )
     m_fatal("producealg: papt is NULL.\n");
  if ( algopt == NULL )
     m_fatal("producealg: algopt is NULL.\n");

  printf("Match Score  Mismatch Score  Gap-Open Penalty Gap-Extension Penalty\n");
  printf("    %d          %d           %d               %d\n\n",
        papt->match, papt->mismat, papt->gopen, papt->gext);

  printf("Sequence A: %s\n", sonept->name);
  printf("    Length: %d\n", sonept->slen);
  printf("Sequence B: %s\n", stwopt->name);
  printf("    Length: %d\n\n", stwopt->slen);

  printf("Alignment Score: %d\n", algopt->score);
  printf("         Length: %d\n\n", algopt->alen);

  i = j = 1;
  pos = 0;
  top = algopt->top;
  bot = algopt->bot;
  mid = algopt->mid;
  while ( pos < algopt->alen )
   {
     limit = pos + linelen;
     if ( limit > algopt->alen )
        limit = algopt->alen;
     tpt = tline;
     mpt = mline;
     bpt = bline;
     spt = sline;
     ast = i;
     bst = j;
     s = 0;

     // makes an alignment section.
     for ( k = pos; k < limit; k++ )
      {
        *tpt++ = a = top[k];
        *mpt++ = mid[k];
        *bpt++ = b = bot[k];
	if ( a != ' ' ) i++;
	if ( b != ' ' ) j++;
	if ( ++s % 10 == 0 )
	  *spt++ = ':';
        else
	  if ( s % 5 == 0 )
	     *spt++ = '.';
          else
	     *spt++ = ' ';
      }
     *tpt = '\0';
     *mpt = '\0';
     *bpt = '\0';
     *spt = '\0';

     printf("%9d %s\n", pos+1, sline);
     printf("%9d %s\n", ast, tline); // reports it
     printf("          %s\n", mline);
     printf("%9d %s\n", bst, bline);
     printf("\n");
     pos = limit; // updates the current alignment position
   }
}

// Checks whether realloc() succeeds if amount is not 0.  
 void *m_check_realloc(void *ptr, size_t amount)
{
  void *tpt;
  // Returns the old memory block and allocates a new one in amount bytes.
  tpt = realloc(ptr, amount);

  /* Checks if it was successful. */
  if ( amount != 0 && tpt == NULL )
   { // prints a message to standard error device (console).
     fprintf(stderr, "No memory of %lu bytes\n", amount);
//   printf("No memory of %lu bytes\n", amount);
     exit(1); // stops unexpectedly.
   }

  return tpt;
}

// Checks if malloc() succeeds.
void *m_check_malloc(size_t amount)
{
  void *tpt;
  /* Allocates a memory block in amount bytes. */
  tpt = malloc( amount );

  /* Checks if it was successful. */
  if ( tpt == NULL )
   { // prints a message to standard error device (console).
     fprintf(stderr, "No memory of %lu bytes\n", amount);
//   printf("No memory of %lu bytes\n", amount);
     exit(1); // stops unexpectedly.
   }

  return tpt;
}

// Prints out an error message and stops.
void m_fatal(const char *msg)
{
     fprintf(stderr, "Error message: %s\n", msg);
//   printf("Error message: %s\n", msg);
     exit(1); // stops unexpectedly.
}

/* ckopen - open file; check for success */
FILE *m_ckopen(const char *name, const char *mode)
{
	FILE *fp;

	if ((fp = fopen(name, mode)) == NULL)
	{ fprintf(stderr, "Cannot open %s\n", name);
	  exit(1);
	}
	return(fp);
}


 /**
  * By having the TEST macro gtest will be able to find the tests
  * with out them having to be registered
  */
TEST(GETSEQ, Check_Name_5PTS){
    // Arrange
   struct seqtp seq_act; 
   struct seqtp seq_exp;

   // Act
   m_getseq("tests/A1.txt",&seq_exp);
   getseq("tests/A1.txt", &seq_act);

   // Assert
   ASSERT_STREQ (seq_exp.name,seq_act.name);
}

TEST(GETSEQ, Check_Seq_7PTS){
    // Arrange
    struct seqtp seq_exp;
    struct seqtp seq_act;

    // Act
    m_getseq("tests/A1.txt",&seq_exp);
    seq_exp.seq[12]=0;
    seq_exp.seq = seq_exp.seq+1;

    getseq("tests/A1.txt",&seq_act);

    if(strlen(seq_act.seq)==12)
        seq_act.seq = seq_act.seq+1;
    if(strlen(seq_act.seq)==13){
        seq_act.seq[12]=0;
        seq_act.seq = seq_act.seq+1;
    }

    // Assert
    ASSERT_STREQ(seq_exp.seq,seq_act.seq);
}

TEST(GETSEQ, Check_Length_3PTS){
    // Arrange
    struct seqtp seq_act;
    struct seqtp seq_exp;

    // Act
    m_getseq("tests/A1.txt",&seq_exp);
    getseq("tests/A1.txt",&seq_act);

    // Assert
    ASSERT_EQ(seq_exp.slen,seq_act.slen);
}

TEST(COMPUTE_MAT, Check_S_MN_2PT){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

    m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);

   // Assert
   ASSERT_EQ(mats_exp.Spt[seq1.slen][seq2.slen],mats_act.Spt[seq1.slen][seq2.slen]);
}

TEST(COMPUTE_MAT, Check_S_LastRow_PARTIAL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq2.slen; i++){
       if(mats_act.Spt[seq1.slen][i] == mats_exp.Spt[seq1.slen][i]){
           count++;
       }
   }

   ASSERT_EQ(true, count >= seq2.slen/2);
}

TEST(COMPUTE_MAT, Check_S_LastRow_FULL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq2.slen; i++){
       if(mats_act.Spt[seq1.slen][i] == mats_exp.Spt[seq1.slen][i]){
           count++;
       }
   }

   ASSERT_EQ(true, count == seq2.slen);
}

TEST(COMPUTE_MAT, Check_S_LastCol_PARTIAL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       if(mats_act.Spt[i][seq2.slen] == mats_exp.Spt[i][seq2.slen]){
           count++;
       }
   }

   ASSERT_EQ(true, count >= seq2.slen/2);
}

TEST(COMPUTE_MAT, Check_S_LastCol_FULL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       if(mats_act.Spt[i][seq2.slen] == mats_exp.Spt[i][seq2.slen]){
           count++;
       }
   }

   ASSERT_EQ(true, count == seq2.slen);
}

TEST(COMPUTE_MAT, Check_S_Full_PARTIAL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int j = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       for(j =0; j< seq2.slen; j++){
            if(mats_exp.Spt[i][j] == mats_exp.Spt[i][j]){
                count++;
            }
       }
   }

   ASSERT_EQ(true, count >= seq2.slen);
}

TEST(COMPUTE_MAT, Check_S_Full_FULL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int j = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       for(j =0; j< seq2.slen; j++){
            if(mats_exp.Spt[i][j] == mats_exp.Spt[i][j]){
                count++;
            }
       }
   }

   ASSERT_EQ(true, count == seq2.slen*seq1.slen);
}

TEST(COMPUTE_MAT, Check_I_MN){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

    m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);

   // Assert
   ASSERT_EQ(mats_exp.Ipt[seq1.slen][seq2.slen],mats_act.Ipt[seq1.slen][seq2.slen]);
}

TEST(COMPUTE_MAT, Check_I_LastRow_PARTIAL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq2.slen; i++){
       if(mats_act.Ipt[seq1.slen][i] == mats_exp.Ipt[seq1.slen][i]){
           count++;
       }
   }

   ASSERT_EQ(true, count >= seq2.slen/2);
}

TEST(COMPUTE_MAT, Check_I_LastRow_FULL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq2.slen; i++){
       if(mats_act.Ipt[seq1.slen][i] == mats_exp.Ipt[seq1.slen][i]){
           count++;
       }
   }

   ASSERT_EQ(true, count == seq2.slen);
}

TEST(COMPUTE_MAT, Check_I_LastCol_PARTIAL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       if(mats_act.Ipt[i][seq2.slen] == mats_exp.Ipt[i][seq2.slen]){
           count++;
       }
   }

   ASSERT_EQ(true, count >= seq2.slen/2);
}

TEST(COMPUTE_MAT, Check_I_LastCol_FULL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       if(mats_act.Ipt[i][seq2.slen] == mats_exp.Ipt[i][seq2.slen]){
           count++;
       }
   }

   ASSERT_EQ(true, count == seq2.slen);
}

TEST(COMPUTE_MAT, Check_I_Full_PARTIAL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int j = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       for(j =0; j< seq2.slen; j++){
            if(mats_exp.Ipt[i][j] == mats_exp.Ipt[i][j]){
                count++;
            }
       }
   }

   ASSERT_EQ(true, count >= seq2.slen);
}

TEST(COMPUTE_MAT, Check_I_Full_FULL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int j = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       for(j =0; j< seq2.slen; j++){
            if(mats_exp.Ipt[i][j] == mats_exp.Ipt[i][j]){
                count++;
            }
       }
   }

   ASSERT_EQ(true, count == seq2.slen*seq1.slen);
}

TEST(COMPUTE_MAT, Check_D_MN){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

    m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);

   // Assert
   ASSERT_EQ(mats_exp.Dpt[seq1.slen][seq2.slen],mats_act.Dpt[seq1.slen][seq2.slen]);
}

TEST(COMPUTE_MAT, Check_D_LastRow_PARTIAL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq2.slen; i++){
       if(mats_act.Dpt[seq1.slen][i] == mats_exp.Dpt[seq1.slen][i]){
           count++;
       }
   }

   ASSERT_EQ(true, count >= seq2.slen/2);
}

TEST(COMPUTE_MAT, Check_D_LastRow_FULL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq2.slen; i++){
       if(mats_act.Dpt[seq1.slen][i] == mats_exp.Dpt[seq1.slen][i]){
           count++;
       }
   }

   ASSERT_EQ(true, count == seq2.slen);
}

TEST(COMPUTE_MAT, Check_D_LastCol_PARTIAL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       if(mats_act.Dpt[i][seq2.slen] == mats_exp.Dpt[i][seq2.slen]){
           count++;
       }
   }

   ASSERT_EQ(true, count >= seq2.slen/2);
}

TEST(COMPUTE_MAT, Check_D_LastCol_FULL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       if(mats_act.Dpt[i][seq2.slen] == mats_exp.Dpt[i][seq2.slen]){
           count++;
       }
   }

   ASSERT_EQ(true, count == seq2.slen);
}

TEST(COMPUTE_MAT, Check_D_Full_PARTIAL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int j = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       for(j =0; j< seq2.slen; j++){
            if(mats_exp.Dpt[i][j] == mats_exp.Dpt[i][j]){
                count++;
            }
       }
   }

   ASSERT_EQ(true, count >= seq2.slen);
}

TEST(COMPUTE_MAT, Check_D_Full_FULL){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats_exp;
   struct mattp mats_act;
   int i = 0;
   int j = 0;
   int count =0;

   // Act
   m_getseq("tests/A1.txt",&seq1);
   m_getseq("tests/A1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats_exp);
   computemats(&seq1,&seq2,&scores,&mats_act);
   // count matches from expected and actual
   for(i =0; i < seq1.slen; i++){
       for(j =0; j< seq2.slen; j++){
            if(mats_exp.Dpt[i][j] == mats_exp.Dpt[i][j]){
                count++;
            }
       }
   }

   ASSERT_EQ(true, count == seq2.slen*seq1.slen);
}

TEST(PRODUCE_ALG, Check_Alignment_Del_Gap){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats;
   int i = 0;
   int j = 0;
   int count =0;
   struct algtp alg_exp;
   struct algtp alg_act;

   // Act
   m_getseq("tests/gap2.txt",&seq1);
   m_getseq("tests/gap1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats);
   
   m_producealg(&seq1,&seq2,&mats,&alg_exp);
   producealg(&seq1,&seq2,&mats,&alg_act);

   // Assert
   ASSERT_STREQ(alg_exp.top,alg_act.top);
   ASSERT_STREQ(alg_exp.mid,alg_act.mid);
   ASSERT_STREQ(alg_exp.bot,alg_act.bot);

}

TEST(PRODUCE_ALG, Check_Alignment_In_Gap){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats;
   int i = 0;
   int j = 0;
   int count =0;
   struct algtp alg_exp;
   struct algtp alg_act;

   // Act
   m_getseq("tests/gap1.txt",&seq1);
   m_getseq("tests/gap2.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats);
   
   m_producealg(&seq1,&seq2,&mats,&alg_exp);
   producealg(&seq1,&seq2,&mats,&alg_act);

   // Assert
   ASSERT_STREQ(alg_exp.top,alg_act.top);
   ASSERT_STREQ(alg_exp.mid,alg_act.mid);
   ASSERT_STREQ(alg_exp.bot,alg_act.bot);
}

TEST(PRODUCE_ALG, CHECK_Alignment_Perfect_Match){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats;
   int i = 0;
   int j = 0;
   int count =0;
   struct algtp alg_exp;
   struct algtp alg_act;

   // Act
   m_getseq("tests/gap1.txt",&seq1);
   m_getseq("tests/gap1.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats);
   
   m_producealg(&seq1,&seq2,&mats,&alg_exp);
   producealg(&seq1,&seq2,&mats,&alg_act);

   // Assert
   ASSERT_STREQ(alg_exp.top,alg_act.top);
   ASSERT_STREQ(alg_exp.mid,alg_act.mid);
   ASSERT_STREQ(alg_exp.bot,alg_act.bot);
}

TEST(PRODUCE_ALG, CHECK_Alignment_Genral_Case){
   // Arrange
   struct seqtp seq1;
   struct seqtp seq2;
   struct scoretp scores = {10,-20,40,2};
   struct mattp mats;
   int i = 0;
   int j = 0;
   int count =0;
   struct algtp alg_exp;
   struct algtp alg_act;

   // Act
   m_getseq("tests/test1.txt",&seq1);
   m_getseq("tests/test2.txt",&seq2);

   m_computemats(&seq1,&seq2,&scores,&mats);
   
   m_producealg(&seq1,&seq2,&mats,&alg_exp);
   producealg(&seq1,&seq2,&mats,&alg_act);

   // Assert
   ASSERT_STREQ(alg_exp.top,alg_act.top);
   ASSERT_STREQ(alg_exp.mid,alg_act.mid);
   ASSERT_STREQ(alg_exp.bot,alg_act.bot);
}




