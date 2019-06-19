#ifndef __HASH_H__
#define __HASH_H__

typedef struct _hashelem {
    char *name;
    struct _hashelem *next;
    struct _column *col;
    struct _bound *bnd;
    int must_be_int;
    int index; /* for row and column name hash tables */
} hashelem;

typedef struct _hashstruct {
    hashelem **table;
    int size;
} hashstruct;

hashstruct *create_hash_table(int size);
void free_hash_table(hashstruct *ht);
hashelem *findhash(const char *name, hashstruct *ht);
hashelem *puthash(const char *name, hashstruct *ht);
hashstruct *copy_hash_table(hashstruct *ht);

#endif
