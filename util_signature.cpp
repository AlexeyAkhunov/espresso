/* Module:util_signature.c
 * Purpose:
 *	contains miscellaneous utility routines 
 * Routines:
 * void set_time_limit(seconds)
 *	sets the cpu time limit
 * int sf_equal():
 *	 checks equlaity of two set families.
 * int print_cover():
 *	prints cover.
 * int mem_usage():
 *	current Memory usage.
 *	Initialized on the first call.
 * int time_usage():
 *	current time usage.
 *	Initialized on the first call.
 */

#include <stdio.h>
#include <math.h>
#include "espresso.h"
#include "signature.h"
#include <sys/time.h>
#include <sys/resource.h>

void
set_time_limit(int seconds)
{
	struct rlimit rlp_st, *rlp = &rlp_st;
	rlp->rlim_cur = seconds;
	setrlimit(RLIMIT_CPU, rlp);
}

void
print_cover(pcover F, char * name)
{
    pcube last, p;
	printf("%s:\t %d\n",name,F->count);
    foreach_set(F,last,p){
        print_cube(std::cout,p,"~0");
	}
	printf("\n\n");
}

/* sf_equal: Check equality of two set families */
int
sf_equal(pcover F1, pcover F2)
{
    int i;
    int count = F1->count;
    pcube *list1,*list2;
    
    if(F1->count != F2->count){
        return(FALSE);
    }
	
	list1 = sf_sort(F1,descend);
	list2 = sf_sort(F2,descend);

	for(i = 0; i < count; i++){
		if(!setp_equal(list1[i],list2[i])){
			return false;
		}
	}

	return TRUE;
}

/* mem_usage:
 * 	Initialized on first call. Prints current memory usage.
 */
int
mem_usage(char * name)
{
	static int memory_init;
	int memory_current;
	static int flag = 1;

	if(flag){
		memory_init = 0;
		flag = 0;
	}

	memory_current = 0;

	printf("Memory %s\t %d\n", name, memory_current - memory_init);

	return memory_current;

}

/* time_usage:
 * 	Initialized on first call. Prints current time usage.
 */
int
time_usage(char * name)
{
    static std::chrono::time_point<std::chrono::high_resolution_clock> time_init;
	static int flag = 1;

	if(flag){
		time_init = ptime();
		flag = 0;
		return 0;
	}

	const auto time_current = ptime();

	printf("%s\t %ld\n", name, static_cast<long>((time_current - time_init).count()/1000.0));

	return static_cast<int>((time_current - time_init).count()/1000.0);

}

/* s_totals : add time spent in the function and update call count */
void
s_totals(const std::chrono::time_point<std::chrono::high_resolution_clock>& time, int i)
{
    total_time[i] += (ptime() - time);
    total_calls[i]++;
}


