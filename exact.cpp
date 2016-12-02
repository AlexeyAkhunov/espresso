#include "espresso.h"
#include <memory>
#include <fstream>

static void dump_irredundant(pcover E, pcover Rt, pcover Rp, sm_matrix * table);
static pcover do_minimize(pcover F, pcover D, pcover R, int exact_cover, int weighted);


/*
 *  minimize_exact -- main entry point for exact minimization
 *
 *  Global flags which affect this routine are:
 *
 *      debug
 *      skip_make_sparse
 */

pcover
minimize_exact(pcover F, pcover D, pcover R, int exact_cover)
{
    return do_minimize(F, D, R, exact_cover, /*weighted*/ 0);
}


pcover
minimize_exact_literals(pcover F, pcover D, pcover R, int exact_cover)
{
    return do_minimize(F, D, R, exact_cover, /*weighted*/ 1);
}



static pcover
do_minimize(pcover F, pcover D, pcover R, int exact_cover, int weighted)
{
    pcover newF, E, Rt, Rp;
    pset p, last;
    int heur, level, *weights;
    sm_matrix *table;
    sm_row *cover;
    sm_element *pe;
    int debug_save = debug;

    if (debug & EXACT) {
	debug |= (IRRED | MINCOV);
    }
#if defined(sun) || defined(bsd4_2)			/* hack ... */
    if (debug & MINCOV) {
	setlinebuf(stdout);
    }
#endif
    level = (debug & MINCOV) ? 4 : 0;
    heur = ! exact_cover;

    /* Generate all prime implicants */
    EXEC(F = primes_consensus(cube2list(F, D)), "PRIMES     ", F);

    /* Setup the prime implicant table */
    EXEC(irred_split_cover(F, D, &E, &Rt, &Rp), "ESSENTIALS ", E);
    EXEC(table = irred_derive_table(D, E, Rp),  "PI-TABLE   ", Rp);

    /* Solve either a weighted or nonweighted covering problem */
    if (weighted) {
	/* correct only for all 2-valued variables */
	weights = new int[F->count];
	foreach_set(Rp, last, p) {
	    weights[SIZE(p)] = cube.size - set_ord(p);
	}
    } else {
	weights = nullptr;
    }
    EXEC(cover=sm_minimum_cover(table,weights,heur,level), "MINCOV     ", F);
    if (weights != 0) {
	delete weights;
    }

    if (debug & EXACT) {
	dump_irredundant(E, Rt, Rp, table);
    }

    /* Form the result cover */
    newF = new_cover(100);
    foreach_set(E, last, p) {
	newF = sf_addset(newF, p);
    }
    sm_foreach_row_element(cover, pe) {
	newF = sf_addset(newF, GETSET(F, pe->col_num));
    }

    free_cover(E);
    free_cover(Rt);
    free_cover(Rp);
    sm_free(table);
    sm_row_free(cover);
    free_cover(F);

    /* Attempt to make the results more sparse */
    debug &= ~ (IRRED | SHARP | MINCOV);
    if (! skip_make_sparse && R != 0) {
	newF = make_sparse(newF, D, R);
    }

    debug = debug_save;
    return newF;
}

static void
dump_irredundant(pcover E, pcover Rt, pcover Rp, sm_matrix * table)
{
    std::unique_ptr<std::ofstream> fp_pi_table_ptr;
    std::unique_ptr<std::ofstream> fp_primes_ptr;
    std::ostream *fp_pi_table, *fp_primes;
    pPLA PLA;
    pset last, p;
    char *file;
    
    if (filename == 0 || strcmp(filename, "(stdin)") == 0) {
        fp_pi_table = fp_primes = &std::cout;
    } else {
	file = new char[strlen(filename)+20];
	(void) sprintf(file, "%s.primes", filename);
        fp_primes_ptr = std::unique_ptr<std::ofstream>(new std::ofstream(file));
        if (fp_primes_ptr -> bad()) {
            fprintf(stderr, "espresso: Unable to open %s\n", file);
            fp_primes = &std::cout;
        } else {
            fp_primes = fp_primes_ptr.get();
        }
	(void) sprintf(file, "%s.pi", filename);
        fp_pi_table_ptr = std::unique_ptr<std::ofstream>(new std::ofstream(file));
	if (fp_pi_table_ptr -> bad()) {
	    fprintf(stderr, "espresso: Unable to open %s\n", file);
        fp_pi_table = &std::cout;
    } else {
        fp_pi_table = fp_pi_table_ptr.get();
    }
	delete file;
    }

    PLA = new_PLA();
    PLA_labels(PLA);

    fpr_header(*fp_primes, PLA, F_type);
    free_PLA(PLA);

    *fp_primes << "# Essential primes are\n";
    foreach_set(E, last, p) {
	*fp_primes << pc1(p) << "\n";
    }
    *fp_primes << "# Totally redundant primes are\n";
    foreach_set(Rt, last, p) {
	*fp_primes << pc1(p) << "\n";
    }
    *fp_primes << "# Partially redundant primes are\n";
    foreach_set(Rp, last, p) {
	*fp_primes << pc1(p) << "\n";
    }
	
    sm_write(*fp_pi_table, table);
}
