#ifndef MASTER_CODE_H
#define MASTER_CODE_H
#include "template.project1.part1.h"

void m_fatal(const char *msg);
FILE *m_ckopen(const char *name, const char *mode);
void *m_check_malloc(size_t amount);
void *m_check_realloc(void *ptr, size_t amount);
void m_getseq(const char *fname, struct seqtp *seqrec);
void **m_get2dspace(int rind, int cind);
void m_computemats(struct seqtp *sonept, struct seqtp *stwopt, struct scoretp *papt, struct mattp *matspt);
void m_producealg(struct seqtp *sonept, struct seqtp *stwopt, struct mattp *matspt, struct algtp *algopt);
void m_outputmat(char kind, struct seqtp *sonept, struct seqtp *stwopt, int **Mpt);
void outputalg(struct seqtp *sonept, struct seqtp *stwopt, struct scoretp *papt, struct algtp *algopt);
void m_freeseq(struct seqtp *seqrec);
void m_freealg(struct algtp *algopt);
void m_free2dspace(int rind, int **arr);
void m_freemat(struct mattp *matspt);

#endif
