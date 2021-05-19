/*

  mmt.cpp : expects call from c file mmt_main.c and runs the MMT algorithm

*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include "../color.h"
#include "../graph.h"
#include "../color_defs.h"
#include "../color_parms.h"
#include "../color_private.h"
#include "../bleile/utils/mmt_read.h"

// static functions copied from color_main.c
static int parseargs (int ac, char **av, COLORparms* parms);
static int get_problem_name(char* pname,const char* efname);
// static functions copied from color_jobkiller.c
static void usage (char *f);

int read_graph(int ac, char **av, int *pncount, int *pecount, int **pelist)
{
      int rval = 0;

      COLORproblem colorproblem;
      COLORparms *parms = &(colorproblem.parms);
      colordata *cd = &(colorproblem.root_cd);

      int ncolors = 0;
      COLORset *colorclasses = (COLORset *)NULL;

      rval = COLORprogram_header(ac, av);
      COLORcheck_rval(rval, "Failed in COLORprogram_header");

      COLORproblem_init(&colorproblem);
      cd->id = 0;
      colorproblem.ncolordata = 1;

      rval = parseargs(ac, av, parms);
      if (rval)
          goto CLEANUP;

      //setzt durch Kommandozeilenparameter übergebene upper_bound oder Default
      cd->upper_bound = parms->initial_upper_bound;
      get_problem_name(cd->pname, parms->edgefile);

      printf("Rounding strategy: ");
      switch (parms->rounding_strategy)
      {
      case COLOR_neighbor_rounding:
          printf("neighbor\n");
          break;
      case COLOR_uniform_rounding:
          printf("uniform\n");
          break;
      case COLOR_no_rounding:
          printf("none\n");
          break;
      }

      if (COLORdbg_lvl() > 1)
          printf("Debugging turned on\n");
      fflush(stdout);


      rval = COLORread_dimacs(parms->edgefile, pncount, pecount,
                              pelist, (int **)NULL);
      // rval = COLORread_dimacs(parms->edgefile, pncount, pecount,
      //                         pelist, (int **)NULL);
      COLORcheck_rval(rval, "COLORread_dimacs failed");

      // for (int i = 0; i < 2*(cd->ecount); i++) {
      //   printf("%i\n", (cd->elist)[i]);
      // }

      // pelist = &(cd->elist);
      // pncount = &(cd->ncount);
      // pecount = &(cd->ecount);

      // if (cd->upper_bound > cd->ncount)
      // {
      //     cd->upper_bound = cd->ncount;
      // }
      //
      // if (colorproblem.parms.backupdir)
      // {
      //     recover_colordata(cd, &colorproblem);
      // }

      /*

        Added snippet to make use of already coded functionalities

      */
      // *pncount = cd->ncount;
      // *pecount = cd->ecount;
      // *pelist = cd->elist;


  CLEANUP:
      COLORproblem_free(&colorproblem);
      COLORfree_sets(&colorclasses, &ncolors);
      COLORlp_free_env();

      return rval;
}

static int parseargs (int ac, char **av, COLORparms* parms)
{
    int c;
    int rval = 0;
    int debug = COLORdbg_lvl();

    while ((c = getopt (ac, av, "admpo:r:w:c:u:b:l:s:R:")) != EOF) {
        switch (c) {
        case 'd':
           /* each -d increases the verbosity by one.*/
           ++(debug);
           COLORset_dbg_lvl(debug);
           break;
        case 'o':
           rval = COLORparms_set_outfile(parms,optarg);
           COLORcheck_rval(rval,"Failed in COLORparms_set_outfile");
           break;
        case 'm':
           COLORparms_set_write_mwis(parms,1);
           break;
        case 'r':
           rval = COLORparms_set_cclasses_infile(parms,optarg);
           COLORcheck_rval(rval,"Failed in COLORparms_set_cclasses_infile");
           break;
        case 'w':
           rval = COLORparms_set_cclasses_outfile(parms,optarg);
           COLORcheck_rval(rval,"Failed in COLORparms_set_cclasses_outfile");
           break;
        case 'c':
           rval = COLORparms_set_color_infile(parms,optarg);
           COLORcheck_rval(rval,"Failed in COLORparms_set_color_infile");
           break;
        case 'u':
           rval = COLORparms_set_initial_upper_bound(parms,atoi(optarg));
           COLORcheck_rval(rval,"Failed in COLORparms_set_initial_upper_bound");
           break;
        case 'a':
           parms->upper_bounds_only  = 1;
           parms->branching_strategy = COLOR_hybrid_strategy;
           break;
        case 'p':
           rval = COLORparms_set_parallel(parms,1);
           COLORcheck_rval(rval,"Failed in COLORparms_set_initial_upper_bound");
           break;
        case 'b':
           rval = COLORparms_set_backupdir(parms,optarg);
           COLORcheck_rval(rval,"Failed in COLORparms_set_backupdir");
           break;
        case 'l':
           rval = COLORparms_set_branching_cpu_limit(parms,atof(optarg));
           COLORcheck_rval(rval,"Failed in COLORparms_set_branching_cpu_limit");
           break;
        case 's':
           rval = COLORparms_set_branching_strategy(parms,atoi(optarg));
           COLORcheck_rval(rval,"Failed in COLORparms_set_branching_strategy");
           break;
        case 'R':
           rval = COLORparms_set_rounding_strategy(parms,atoi(optarg));
           COLORcheck_rval(rval,"Failed in COLORparms_set_rounding_strategy");
           break;
        default:
          rval = 1;
           goto CLEANUP;
        }
    }

    if (ac <= optind) {
        rval = 1; goto CLEANUP;
    } else {
       rval =  COLORparms_set_edgefile(parms,av[optind++]);
       COLORcheck_rval(rval,"Failed in COLORparms_set_edgefile");
    }

CLEANUP:

    if (rval) {usage (av[0]);}
    return  (rval);
}

static int get_problem_name(char* pname,const char* efname)
{
  int    rval = 0;
  int    len = 0;
  const char * fname = strrchr(efname,'/');
  const char * lastdot = strrchr(efname,'.');
  if(!fname) {
     /* no slashes in efname.*/
     fname = efname;
  } else {
     fname++;
  }

  if (lastdot) {
     len = lastdot - fname + 1;
  } else {
     len = strlen(fname);
  }

  if (snprintf(pname,len,"%s",fname) < 0) {
     rval = 1;
  }
  printf("Extracted problem name %s\n",pname);

  return rval;
}

static void usage (char *f)
{
   fprintf (stderr, "Usage %s: <coloring job id>\n", f);
}
