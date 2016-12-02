/*
    module: cvrout.c
    purpose: cube and cover output routines
*/

#include <ostream>
#include "espresso.h"

void fprint_pla(std::ostream& fp, pPLA PLA, int output_type)
{
    int num;
    pcube last, p;

    if ((output_type & CONSTRAINTS_type) != 0) {
	output_symbolic_constraints(fp, PLA, 0);
	output_type &= ~ CONSTRAINTS_type;
	if (output_type == 0) {
	    return;
	}
    }

    if ((output_type & SYMBOLIC_CONSTRAINTS_type) != 0) {
	output_symbolic_constraints(fp, PLA, 1);
	output_type &= ~ SYMBOLIC_CONSTRAINTS_type;
	if (output_type == 0) {
	    return;
	}
    }

    if (output_type == PLEASURE_type) {
	pls_output(PLA);
    } else if (output_type == EQNTOTT_type) {
	eqn_output(PLA);
    } else if (output_type == KISS_type) {
	kiss_output(fp, PLA);
    } else {
	fpr_header(fp, PLA, output_type);

	num = 0;
	if (output_type & F_type) num += (PLA->F)->count;
	if (output_type & D_type) num += (PLA->D)->count;
	if (output_type & R_type) num += (PLA->R)->count;
    fp << ".p " << num << "\n";

	/* quick patch 01/17/85 to support TPLA ! */
	if (output_type == F_type) {
	    foreach_set(PLA->F, last, p) {
		print_cube(fp, p, "01");
	    }
        fp << ".p " << num << "\n";
	    fp << ".e\n";
	} else {
	    if (output_type & F_type) {
		foreach_set(PLA->F, last, p) {
		    print_cube(fp, p, "~1");
		}
	    }
	    if (output_type & D_type) {
		foreach_set(PLA->D, last, p) {
		    print_cube(fp, p, "~2");
		}
	    }
	    if (output_type & R_type) {
		foreach_set(PLA->R, last, p) {
		    print_cube(fp, p, "~0");
		}
	    }
	    fp << ".end\n";
	}
    }
}

void fpr_header(std::ostream& fp, pPLA PLA, int output_type)
{
    int i, var;
    int first, last;

    /* .type keyword gives logical type */
    if (output_type != F_type) {
	fp << ".type ";
	if (output_type & F_type) fp.put('f');
	if (output_type & D_type) fp.put('d');
	if (output_type & R_type) fp.put('r');
	fp.put('\n');
    }

    /* Check for binary or multiple-valued labels */
    if (cube.num_mv_vars <= 1) {
	fp << ".i " << cube.num_binary_vars << "\n";
	if (cube.output != -1)
	    fp << ".o " << cube.part_size[cube.output] << "\n";
    } else {
	fp << ".mv " << cube.num_vars << " " << cube.num_binary_vars;
	for(var = cube.num_binary_vars; var < cube.num_vars; var++)
	    fp << " " << cube.part_size[var];
	fp << "\n";
    }

    /* binary valued labels */
    if (PLA->label != nullptr && PLA->label[1] != nullptr
	    && cube.num_binary_vars > 0) {
	fp << ".ilb";
	for(var = 0; var < cube.num_binary_vars; var++)
	    fp << " " << INLABEL(var);
	fp.put('\n');
    }

    /* output-part (last multiple-valued variable) labels */
    if (PLA->label != nullptr &&
	    PLA->label[cube.first_part[cube.output]] != nullptr
		&& cube.output != -1) {
	fp << ".ob";
	for(i = 0; i < cube.part_size[cube.output]; i++)
	    fp << " " << OUTLABEL(i);
	fp.put('\n');
    }

    /* multiple-valued labels */
    for(var = cube.num_binary_vars; var < cube.num_vars-1; var++) {
	first = cube.first_part[var];
	last = cube.last_part[var];
	if (PLA->label != NULL && PLA->label[first] != NULL) {
	    fp << ".label var=" << var;
	    for(i = first; i <= last; i++) {
		fp << " " << PLA->label[i];
	    }
	    fp.put('\n');
	}
    }

    if (PLA->phase != (pcube) NULL) {
	first = cube.first_part[cube.output];
	last = cube.last_part[cube.output];
	fp << "#.phase ";
	for(i = first; i <= last; i++)
	    fp.put(is_in_set(PLA->phase,i) ? '1' : '0');
	fp << "\n";
    }
}

void pls_output(pPLA PLA)
{
    pcube last, p;

    printf(".option unmerged\n");
    makeup_labels(PLA);
    pls_label(PLA, std::cout);
    pls_group(PLA, std::cout);
    printf(".p %d\n", PLA->F->count);
    foreach_set(PLA->F, last, p) {
        print_expanded_cube(std::cout, p, PLA->phase);
    }
    printf(".end\n");
}


void pls_group(pPLA PLA, std::ostream& fp)
{
    int var, i;
    std::size_t col, len;

    fp << "\n.group";
    col = 6;
    for(var = 0; var < cube.num_vars-1; var++) {
	fp << " (", col += 2;
	for(i = cube.first_part[var]; i <= cube.last_part[var]; i++) {
	    len = strlen(PLA->label[i]);
	    if (col + len > 75)
		fp << " \\\n", col = 0;
	    else if (i != 0)
		fp.put(' '), col += 1;
	    fp << PLA->label[i], col += len;
	}
	fp << ")", col += 1;
    }
    fp << "\n";
}


void pls_label(pPLA PLA, std::ostream& fp)
{
    int var, i;
    std::size_t col, len;

    fp << ".label";
    col = 6;
    for(var = 0; var < cube.num_vars; var++)
	for(i = cube.first_part[var]; i <= cube.last_part[var]; i++) {
	    len = strlen(PLA->label[i]);
	    if (col + len > 75)
		fp << " \\\n", col = 0;
	    else
		fp.put(' '), col += 1;
	    fp << PLA->label[i], col += len;
	}
}



/*
    eqntott output mode -- output algebraic equations
*/
void eqn_output(pPLA PLA)
{
    pcube p, last;
    int i, var;
    std::size_t col, len;
    int x;
    bool firstand, firstor;

    if (cube.output == -1)
	fatal("Cannot have no-output function for EQNTOTT output mode");
    if (cube.num_mv_vars != 1)
	fatal("Must have binary-valued function for EQNTOTT output mode");
    makeup_labels(PLA);

    /* Write a single equation for each output */
    for(i = 0; i < cube.part_size[cube.output]; i++) {
	printf("%s = ", OUTLABEL(i));
	col = strlen(OUTLABEL(i)) + 3;
	firstor = TRUE;

	/* Write product terms for each cube in this output */
	foreach_set(PLA->F, last, p)
	    if (is_in_set(p, i + cube.first_part[cube.output])) {
		if (firstor)
		    printf("("), col += 1;
		else
		    printf(" | ("), col += 4;
		firstor = FALSE;
		firstand = TRUE;

		/* print out a product term */
		for(var = 0; var < cube.num_binary_vars; var++)
		    if ((x=GETINPUT(p, var)) != DASH) {
			len = strlen(INLABEL(var));
			if (col+len > 72)
			    printf("\n    "), col = 4;
			if (! firstand)
			    printf("&"), col += 1;
			firstand = FALSE;
			if (x == ZERO)
			    printf("!"), col += 1;
			printf("%s", INLABEL(var)), col += len;
		    }
		printf(")"), col += 1;
	    }
	printf(";\n\n");
    }
}


char *fmt_cube(pcube c, const char * out_map, char * s)
{
    int i, var, last, len = 0;

    for(var = 0; var < cube.num_binary_vars; var++) {
	s[len++] = "?01-" [GETINPUT(c, var)];
    }
    for(var = cube.num_binary_vars; var < cube.num_vars - 1; var++) {
	s[len++] = ' ';
	for(i = cube.first_part[var]; i <= cube.last_part[var]; i++) {
	    s[len++] = "01" [is_in_set(c, i) != 0];
	}
    }
    if (cube.output != -1) {
	last = cube.last_part[cube.output];
	s[len++] = ' ';
	for(i = cube.first_part[cube.output]; i <= last; i++) {
	    s[len++] = out_map [is_in_set(c, i) != 0];
	}
    }
    s[len] = '\0';
    return s;
}


void print_cube(std::ostream& fp, const pcube& c, const char * out_map)
{
    int i, var, ch;
    int last;

    for(var = 0; var < cube.num_binary_vars; var++) {
	ch = "?01-" [GETINPUT(c, var)];
	fp.put(ch);
    }
    for(var = cube.num_binary_vars; var < cube.num_vars - 1; var++) {
	fp.put(' ');
	for(i = cube.first_part[var]; i <= cube.last_part[var]; i++) {
	    ch = "01" [is_in_set(c, i) != 0];
	    fp.put(ch);
	}
    }
    if (cube.output != -1) {
	last = cube.last_part[cube.output];
	fp.put(' ');
	for(i = cube.first_part[cube.output]; i <= last; i++) {
	    ch = out_map [is_in_set(c, i) != 0];
	    fp.put(ch);
	}
    }
    fp.put('\n');
}


void print_expanded_cube(std::ostream& fp, pcube c, pcube phase)
{
    int i, var, ch;
    const char *out_map;

    for(var = 0; var < cube.num_binary_vars; var++) {
	for(i = cube.first_part[var]; i <= cube.last_part[var]; i++) {
	    ch = "~1" [is_in_set(c, i) != 0];
	    fp.put(ch);
	}
    }
    for(var = cube.num_binary_vars; var < cube.num_vars - 1; var++) {
	for(i = cube.first_part[var]; i <= cube.last_part[var]; i++) {
	    ch = "1~" [is_in_set(c, i) != 0];
	    fp.put(ch);
	}
    }
    if (cube.output != -1) {
	var = cube.num_vars - 1;
	fp.put(' ');
	for(i = cube.first_part[var]; i <= cube.last_part[var]; i++) {
	    if (phase == (pcube) NULL || is_in_set(phase, i)) {
		out_map = "~1";
	    } else {
		out_map = "~0";
	    }
	    ch = out_map[is_in_set(c, i) != 0];
	    fp.put(ch);
	}
    }
    fp.put('\n');
}


char *pc1(pcube c)
{static char s1[256];return fmt_cube(c, "01", s1);}
char *pc2(pcube c)
{static char s2[256];return fmt_cube(c, "01", s2);}


void debug_print(pcube * T, const char * name, int level)
{
    pcube *T1, p, temp;
    long cnt;

    cnt = CUBELISTSIZE(T);
    temp = new_cube();
    if (verbose_debug && level == 0)
	printf("\n");
    printf("%s[%d]: ord(T)=%ld\n", name, level, cnt);
    if (verbose_debug) {
	printf("cofactor=%s\n", pc1(T[0]));
	for(T1 = T+2, cnt = 1; (p = *T1++) != (pcube) NULL; cnt++)
	    printf("%4ld. %s\n", cnt, pc1(set_or(temp, p, T[0])));
    }
    free_cube(temp);
}


void debug1_print(pcover T, const char * name, int num)
{
    int cnt = 1;
    pcube p, last;

    if (verbose_debug && num == 0)
	printf("\n");
    printf("%s[%d]: ord(T)=%d\n", name, num, T->count);
    if (verbose_debug)
	foreach_set(T, last, p)
	    printf("%4d. %s\n", cnt++, pc1(p));
}


void cprint(pcover T)
{
    pcube p, last;

    foreach_set(T, last, p)
	printf("%s\n", pc1(p));
}


void makeup_labels(pPLA PLA)
{
    int var, i, ind;

    if (PLA->label == (char **) NULL)
	PLA_labels(PLA);

    for(var = 0; var < cube.num_vars; var++)
	for(i = 0; i < cube.part_size[var]; i++) {
	    ind = cube.first_part[var] + i;
	    if (PLA->label[ind] == (char *) NULL) {
		PLA->label[ind] = new char[15];
		if (var < cube.num_binary_vars)
		    if ((i % 2) == 0)
			(void) sprintf(PLA->label[ind], "v%d.bar", var);
		    else
			(void) sprintf(PLA->label[ind], "v%d", var);
		else
		    (void) sprintf(PLA->label[ind], "v%d.%d", var, i);
	    }
	}
}

void kiss_output(std::ostream& fp, pPLA PLA)
{
    pset last, p;

    foreach_set(PLA->F, last, p) {
	kiss_print_cube(fp, PLA, p, "~1");
    }
    foreach_set(PLA->D, last, p) {
	kiss_print_cube(fp, PLA, p, "~2");
    }
}


void kiss_print_cube(std::ostream& fp, pPLA PLA, pcube p, const char* out_string)
{
    int i, var;
    int part, x;

    for(var = 0; var < cube.num_binary_vars; var++) {
	x = "?01-" [GETINPUT(p, var)];
	fp.put(x);
    }

    for(var = cube.num_binary_vars; var < cube.num_vars - 1; var++) {
	fp.put(' ');
	if (setp_implies(cube.var_mask[var], p)) {
	    fp.put('-');
	} else {
	    part = -1;
	    for(i = cube.first_part[var]; i <= cube.last_part[var]; i++) {
		if (is_in_set(p, i)) {
		    if (part != -1) {
			fatal("more than 1 part in a symbolic variable\n");
		    }
		    part = i;
		}
	    }
	    if (part == -1) {
		fp.put('~');	/* no parts, hope its an output ... */
	    } else {
		fp << PLA->label[part];
	    }
	}
    }

    if ((var = cube.output) != -1) {
	fp.put(' ');
	for(i = cube.first_part[var]; i <= cube.last_part[var]; i++) {
	    x = out_string [is_in_set(p, i) != 0];
	    fp.put(x);
	}
    }

    fp.put('\n');
}

void output_symbolic_constraints(std::ostream& fp, pPLA PLA, int output_symbolic)
{
    pset_family A;
    int i, j;
    int size, var, npermute, *permute, *weight, noweight;

    if ((cube.num_vars - cube.num_binary_vars) <= 1) {
	return;
    }
    makeup_labels(PLA);

    for(var=cube.num_binary_vars; var < cube.num_vars-1; var++) {

	/* pull out the columns for variable "var" */
	npermute = cube.part_size[var];
	permute = new int[npermute];
	for(i=0; i < npermute; i++) {
	    permute[i] = cube.first_part[var] + i;
	}
	A = sf_permute(sf_save(PLA->F), permute, npermute);
	delete permute;


	/* Delete the singletons and the full sets */
	noweight = 0;
	for(i = 0; i < A->count; i++) {
	    size = set_ord(GETSET(A,i));
	    if (size == 1 || size == A->sf_size) {
		sf_delset(A, i--);
		noweight++;
	    }
	}


	/* Count how many times each is duplicated */
	weight = new int[A->count];
	for(i = 0; i < A->count; i++) {
	    RESET(GETSET(A, i), COVERED);
	}
	for(i = 0; i < A->count; i++) {
	    weight[i] = 0;
	    if (! TESTP(GETSET(A,i), COVERED)) {
		weight[i] = 1;
		for(j = i+1; j < A->count; j++) {
		    if (setp_equal(GETSET(A,i), GETSET(A,j))) {
			weight[i]++;
			SET(GETSET(A,j), COVERED);
		    }
		}
	    }
	}


	/* Print out the contraints */
	if (! output_symbolic) {
	    fp <<
	    "# Symbolic constraints for variable " << var << "(Numeric form)\n";
	    fp << "# unconstrained weight = " << noweight << "\n";
	    fp << "num_codes=" << cube.part_size[var] << "\n";
	    for(i = 0; i < A->count; i++) {
		if (weight[i] > 0) {
		    fp << "weight=" << weight[i] << ": ";
		    for(j = 0; j < A->sf_size; j++) {
			if (is_in_set(GETSET(A,i), j)) {
			    fp << " " << j;
			}
		    }
		    fp << "\n";
		}
	    }
	} else {
	    fp <<
	    "# Symbolic constraints for variable " << var << " (Symbolic form)\n";
	    for(i = 0; i < A->count; i++) {
		if (weight[i] > 0) {
		    fp << "#   w=" << weight[i] << ": (";
		    for(j = 0; j < A->sf_size; j++) {
			if (is_in_set(GETSET(A,i), j)) {
			    fp << " " <<
				PLA->label[cube.first_part[var]+j];
			}
		    }
		    fp << " )\n";
		}
	    }
	    delete weight;
	}
    }
}
