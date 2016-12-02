#include "espresso.h"

/*
    The cofactor of a cover against a cube "c" is a cover formed by the
    cofactor of each cube in the cover against c.  The cofactor of two
    cubes is null if they are distance 1 or more apart.  If they are
    distance zero apart, the cofactor is the restriction of the cube
    to the minterms of c.

    The cube list contains the following information:

	T[0] = pointer to a cube identifying the variables that have
		been cofactored against
	T[1] = pointer to just beyond the sentinel (i.e., T[n] in this case)
	T[2]
	  .
	  .  = pointers to cubes
	  .
	T[n-2]
	T[n-1] = NULL pointer (sentinel)


    Cofactoring involves repeated application of "cdist0" to check if a
    cube of the cover intersects the cofactored cube.  This can be
    slow, especially for the recursive descent of the espresso
    routines.  Therefore, a special cofactor routine "scofactor" is
    provided which assumes the cofactor is only in a single variable.
*/


/* cofactor -- compute the cofactor of a cover with respect to a cube */
pcube *cofactor(pcube* T, pcube c)
{
    pcube temp = cube.temp[0], *Tc_save, *Tc, *T1;
    pcube p;
    int listlen;

    listlen = static_cast<int>(CUBELISTSIZE(T) + 5);

    /* Allocate a new list of cube pointers (max size is previous size) */
    Tc_save = Tc = new pcube[listlen];

    /* pass on which variables have been cofactored against */
    *Tc++ = set_or(new_cube(), T[0], set_diff(temp, cube.fullset, c));
    Tc++;

    /* Loop for each cube in the list, determine suitability, and save */
    for(T1 = T+2; (p = *T1++) != NULL; ) {
	if (p != c) {

#ifdef NO_INLINE
	if (! cdist0(p, c)) goto false;
#else
    {int w,last;unsigned int x;if((last=cube.inword)!=-1)
    {x=p[last]&c[last];if(~(x|x>>1)&cube.inmask)goto false0;for(w=1;w<last;w++)
    {x=p[w]&c[w];if(~(x|x>>1)&DISJOINT)goto false0;}}}{int w,var,last;
    pcube mask;for(var=cube.num_binary_vars;var<cube.num_vars;var++){
    mask=cube.var_mask[var];last=cube.last_word[var];for(w=cube.first_word[var
    ];w<=last;w++)if(p[w]&c[w]&mask[w])goto nextvar;goto false0;nextvar:;}}
#endif

	    *Tc++ = p;
	false0: ;
	}
    }

    *Tc++ = (pcube) NULL;                       /* sentinel */
    Tc_save[1] = (pcube) Tc;                    /* save pointer to last */
    return Tc_save;
}

/*
    scofactor -- compute the cofactor of a cover with respect to a cube,
    where the cube is "active" in only a single variable.

    This routine has been optimized for speed.
*/

pcube *scofactor(pcube * T, pcube c, int var)
{
    pcube *Tc, *Tc_save;
    pcube p, mask = cube.temp[1], *T1;
    int first = cube.first_word[var], last = cube.last_word[var];
    int listlen;

    listlen = static_cast<int>(CUBELISTSIZE(T) + 5);

    /* Allocate a new list of cube pointers (max size is previous size) */
    Tc_save = Tc = new pcube[listlen];

    /* pass on which variables have been cofactored against */
    *Tc++ = set_or(new_cube(), T[0], set_diff(mask, cube.fullset, c));
    Tc++;

    /* Setup for the quick distance check */
    (void) set_and(mask, cube.var_mask[var], c);

    /* Loop for each cube in the list, determine suitability, and save */
    for(T1 = T+2; (p = *T1++) != NULL; )
	if (p != c) {
	    int i = first;
	    do
		if (p[i] & mask[i]) {
		    *Tc++ = p;
		    break;
		}
	    while (++i <= last);
	}

    *Tc++ = (pcube) NULL;                       /* sentinel */
    Tc_save[1] = (pcube) Tc;                    /* save pointer to last */
    return Tc_save;
}

void massive_count(pcube * T)
{
    int *count = cdata.part_zeros;
    pcube *T1;

    /* Clear the column counts (count of # zeros in each column) */
 {  int i;
    for(i = cube.size - 1; i >= 0; i--)
	count[i] = 0;
 }

    /* Count the number of zeros in each column */
 {  int i, *cnt;
    unsigned int val;
    pcube p, cof = T[0], full = cube.fullset;
    for(T1 = T+2; (p = *T1++) != NULL; )
	for(i = LOOP(p); i > 0; i--)
	    if ((val = full[i] & ~ (p[i] | cof[i]))) {
		cnt = count + ((i-1) << LOGBPI);
#if BPI == 32
	    if (val & 0xFF000000) {
		if (val & 0x80000000) cnt[31]++;
		if (val & 0x40000000) cnt[30]++;
		if (val & 0x20000000) cnt[29]++;
		if (val & 0x10000000) cnt[28]++;
		if (val & 0x08000000) cnt[27]++;
		if (val & 0x04000000) cnt[26]++;
		if (val & 0x02000000) cnt[25]++;
		if (val & 0x01000000) cnt[24]++;
	    }
	    if (val & 0x00FF0000) {
		if (val & 0x00800000) cnt[23]++;
		if (val & 0x00400000) cnt[22]++;
		if (val & 0x00200000) cnt[21]++;
		if (val & 0x00100000) cnt[20]++;
		if (val & 0x00080000) cnt[19]++;
		if (val & 0x00040000) cnt[18]++;
		if (val & 0x00020000) cnt[17]++;
		if (val & 0x00010000) cnt[16]++;
	    }
#endif
	    if (val & 0xFF00) {
		if (val & 0x8000) cnt[15]++;
		if (val & 0x4000) cnt[14]++;
		if (val & 0x2000) cnt[13]++;
		if (val & 0x1000) cnt[12]++;
		if (val & 0x0800) cnt[11]++;
		if (val & 0x0400) cnt[10]++;
		if (val & 0x0200) cnt[ 9]++;
		if (val & 0x0100) cnt[ 8]++;
	    }
	    if (val & 0x00FF) {
		if (val & 0x0080) cnt[ 7]++;
		if (val & 0x0040) cnt[ 6]++;
		if (val & 0x0020) cnt[ 5]++;
		if (val & 0x0010) cnt[ 4]++;
		if (val & 0x0008) cnt[ 3]++;
		if (val & 0x0004) cnt[ 2]++;
		if (val & 0x0002) cnt[ 1]++;
		if (val & 0x0001) cnt[ 0]++;
	    }
	}
 }

    /*
     * Perform counts for each variable:
     *    cdata.var_zeros[var] = number of zeros in the variable
     *    cdata.parts_active[var] = number of active parts for each variable
     *    cdata.vars_active = number of variables which are active
     *    cdata.vars_unate = number of variables which are active and unate
     *
     *    best -- the variable which is best for splitting based on:
     *    mostactive -- most # active parts in any variable
     *    mostzero -- most # zeros in any variable
     *    mostbalanced -- minimum over the maximum # zeros / part / variable
     */

 {  int var, i, lastbit, active, maxactive;
    int best = -1, mostactive = 0, mostzero = 0, mostbalanced = 32000;
    cdata.vars_unate = cdata.vars_active = 0;

    for(var = 0; var < cube.num_vars; var++) {
        if (var < cube.num_binary_vars) { /* special hack for binary vars */
            i = count[var*2];
            lastbit = count[var*2 + 1];
            active = (i > 0) + (lastbit > 0);
            cdata.var_zeros[var] = i + lastbit;
            maxactive = std::max(i, lastbit);
        } else {
            maxactive = active = cdata.var_zeros[var] = 0;
            lastbit = cube.last_part[var];
            for(i = cube.first_part[var]; i <= lastbit; i++) {
                cdata.var_zeros[var] += count[i];
                active += (count[i] > 0);
                if (active > maxactive) maxactive = active;
            }
        }

        /* first priority is to maximize the number of active parts */
        /* for binary case, this will usually select the output first */
        if (active > mostactive) {
            best = var, mostactive = active, mostzero = cdata.var_zeros[best],
            mostbalanced = maxactive;
        } else if (active == mostactive) {
            /* secondary condition is to maximize the number zeros */
            /* for binary variables, this is the same as minimum # of 2's */
            if (cdata.var_zeros[var] > mostzero) {
                best = var, mostzero = cdata.var_zeros[best],
                mostbalanced = maxactive;
            } else if (cdata.var_zeros[var] == mostzero) {
                /* third condition is to pick a balanced variable */
                /* for binary vars, this means roughly equal # 0's and 1's */
                if (maxactive < mostbalanced) {
                    best = var, mostbalanced = maxactive;
                }
            }
        }

        cdata.parts_active[var] = active;
        cdata.is_unate[var] = (active == 1);
        cdata.vars_active += (active > 0);
        cdata.vars_unate += (active == 1);
    }
    cdata.best = best;
 }
}

int binate_split_select(pcube * T, pcube cleft, pcube cright, int debug_flag)
{
    int best = cdata.best;
    int i, lastbit = cube.last_part[best], halfbit = 0;
    pcube cof=T[0];

    /* Create the cubes to cofactor against */
    set_diff(cleft, cube.fullset, cube.var_mask[best]);
    set_diff(cright, cube.fullset, cube.var_mask[best]);
    for(i = cube.first_part[best]; i <= lastbit; i++)
	if (! is_in_set(cof,i))
	    halfbit++;
    for(i = cube.first_part[best], halfbit = halfbit/2; halfbit > 0; i++)
	if (! is_in_set(cof,i))
	    halfbit--, set_insert(cleft, i);
    for(; i <= lastbit; i++)
	if (! is_in_set(cof,i))
	    set_insert(cright, i);

    if (debug & debug_flag) {
	printf("BINATE_SPLIT_SELECT: split against %d\n", best);
	if (verbose_debug)
	    printf("cl=%s\ncr=%s\n", pc1(cleft), pc2(cright));
    }
    return best;
}


pcube *cube1list(pcover A)
{
    pcube last, p, *plist, *list;

    list = plist = new pcube[A->count + 3];
    *plist++ = new_cube();
    plist++;
    foreach_set(A, last, p) {
	*plist++ = p;
    }
    *plist++ = NULL;                    /* sentinel */
    list[1] = (pcube) plist;
    return list;
}


pcube *cube2list(pcover A, pcover B)
{
    pcube last, p, *plist, *list;

    list = plist = new pcube[A->count + B->count + 3];
    *plist++ = new_cube();
    plist++;
    foreach_set(A, last, p) {
	*plist++ = p;
    }
    foreach_set(B, last, p) {
	*plist++ = p;
    }
    *plist++ = NULL;
    list[1] = (pcube) plist;
    return list;
}


pcube *cube3list(pcover A, pcover B, pcover C)
{
    pcube last, p, *plist, *list;

    plist = new pcube[A->count + B->count + C->count + 3];
    list = plist;
    *plist++ = new_cube();
    plist++;
    foreach_set(A, last, p) {
	*plist++ = p;
    }
    foreach_set(B, last, p) {
	*plist++ = p;
    }
    foreach_set(C, last, p) {
	*plist++ = p;
    }
    *plist++ = NULL;
    list[1] = (pcube) plist;
    return list;
}


pcover cubeunlist(pcube * A1)
{
    int i;
    pcube p, pdest, cof = A1[0];
    pcover A;

    A = new_cover(static_cast<int>(CUBELISTSIZE(A1)));
    for(i = 2; (p = A1[i]) != NULL; i++) {
	pdest = GETSET(A, i-2);
	INLINEset_or(pdest, p, cof);
    }
    A->count = static_cast<int>(CUBELISTSIZE(A1));
    return A;
}

void simplify_cubelist(pcube * T)
{
    pcube *Tdest;
    int i, ncubes;

    set_copy(cube.temp[0], T[0]);		/* retrieve cofactor */

    ncubes = static_cast<int>(CUBELISTSIZE(T));
    qsort((char *) (T+2), ncubes, sizeof(pset), d1_order);

    Tdest = T+2;
    /*   *Tdest++ = T[2];   */
    for(i = 3; i < ncubes; i++) {
	if (d1_order(&T[i-1], &T[i]) != 0) {
	    *Tdest++ = T[i];
	}
    }

    *Tdest++ = NULL;				/* sentinel */
    Tdest[1] = (pcube) Tdest;			/* save pointer to last */
}