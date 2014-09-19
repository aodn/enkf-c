/******************************************************************************
 *
 * File:        grid.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <values.h>
#include <string.h>
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#if !defined(NO_GRIDUTILS)
#include "gridnodes.h"
#endif
#include "gridmap.h"
#include "nan.h"

struct grid {
    char* name;
    int htype;

    grid_xy2fij_fn xy2fij_fn;
    grid_z2fk_fn z2fk_fn;
    grid_fij2xy_fn fij2xy_fn;
    grid_tocartesian_fn tocartesian_fn;
    void* gridnodes;
    int lontype;

    int** numlevels;
    float** depth;
};

typedef struct {
    int nx;
    int ny;
    int nz;
    int periodic_x;
    int periodic_y;
    int regular;

    double* x;
    double* y;
    double* z;
    double* xb;
    double* yb;
    double* zb;
} gnxy;                         /* for "Grid Nodes X Y" */

#if !defined(NO_GRIDUTILS)
typedef struct {
    gridnodes* gn;
    gridmap* gm;

    int nz;
    double* z;
    double* zb;
} gnc;                          /* for "Grid Nodes Curvilinear" */
#endif

/**
 */
static gnxy* gnxy_create(int nx, int ny, int nz, double* x, double* y, double* z, int periodic_x, int periodic_y, int regular)
{
    gnxy* nodes = malloc(sizeof(gnxy));
    int i;

    nodes->nx = nx;
    nodes->ny = ny;
    nodes->nz = nz;
    nodes->x = x;
    nodes->y = y;
    nodes->z = z;
    nodes->periodic_x = periodic_x;
    nodes->periodic_y = periodic_y;
    nodes->regular = regular;

    nodes->zb = malloc((nz + 1) * sizeof(double));
    assert(fabs(z[nz - 1]) > fabs(z[0]));
    nodes->zb[0] = 0;
    for (i = 1; i < nz + 1; ++i)
        nodes->zb[i] = 2 * z[i - 1] - nodes->zb[i - 1];

    if (regular) {
        nodes->xb = NULL;
        nodes->yb = NULL;
    } else {
        nodes->xb = malloc((nx + 1) * sizeof(double));
        nodes->xb[0] = x[0] * 1.5 - x[1] * 0.5;
        for (i = 1; i < nx + 1; ++i)
            nodes->xb[i] = 2 * x[i - 1] - nodes->xb[i - 1];

        nodes->yb = malloc((ny + 1) * sizeof(double));
        nodes->yb[0] = y[0] * 1.5 - y[1] * 0.5;
        for (i = 1; i < ny + 1; ++i)
            nodes->yb[i] = 2 * y[i - 1] - nodes->yb[i - 1];
    }

    return nodes;
}

/**
 */
void gnxy_destroy(gnxy * nodes)
{
    free(nodes->x);
    free(nodes->y);
    free(nodes->z);
    free(nodes->zb);
    if (!nodes->regular) {
        free(nodes->xb);
        free(nodes->yb);
    }
    free(nodes);
}

#if !defined(NO_GRIDUTILS)
/**
 */
static gnc* gnc_create(int nx, int ny, int nz, double** x, double** y, double* z)
{
    gnc* nodes = malloc(sizeof(gnc));
    gridnodes* gn_new;
    int i;

    nodes->nz = nz;
    nodes->zb = malloc((nz + 1) * sizeof(double));
    assert(fabs(z[nz - 1]) > fabs(z[0]));
    nodes->zb[0] = 0;
    for (i = 1; i < nz + 1; ++i)
        nodes->zb[i] = 2 * z[i - 1] - nodes->zb[i - 1];

    nodes->gn = gridnodes_create2(nx, ny, NT_CEN, x, y);
    gn_new = gridnodes_transform(nodes->gn, NT_COR);
    gridnodes_destroy(nodes->gn);
    nodes->gn = gn_new;
    gridnodes_validate(nodes->gn);
    nodes->gm = gridmap_build(gridnodes_getnx(nodes->gn), gridnodes_getny(nodes->gn), gridnodes_getx(nodes->gn), gridnodes_gety(nodes->gn));

    return nodes;
}

/**
 */
void gnc_destroy(gnc * nodes)
{
    gridnodes_destroy(nodes->gn);
    gridmap_destroy(nodes->gm);
    free(nodes->z);
    free(nodes->zb);
    free(nodes);
}
#endif

/**
 */
grid* grid_create(char name[])
{
    grid* g = malloc(sizeof(grid));

    g->name = strdup(name);
    g->htype = GRIDHTYPE_NONE;
    g->xy2fij_fn = NULL;
    g->gridnodes = NULL;
    g->lontype = LONTYPE_NONE;
    g->numlevels = NULL;
    g->depth = NULL;

    return g;
}

/**
 */
void grid_destroy(grid* g)
{
    free(g->name);
    if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR)
        gnxy_destroy(g->gridnodes);
#if !defined(NO_GRIDUTILS)
    else if (g->htype == GRIDHTYPE_CURVILINEAR)
        gnc_destroy(g->gridnodes);
#endif
    else
        enkf_quit("programming_error");
    if (g->numlevels != NULL)
        free2d(g->numlevels);
    if (g->depth != NULL)
        free2d(g->depth);
    free(g);
}

/**
 */
void grid_print(grid* g, char offset[])
{
    int nx, ny, nz;

    enkf_printf("%sgrid info:\n", offset);
    switch (g->htype) {
    case GRIDHTYPE_LATLON_REGULAR:
        enkf_printf("%s  type = LATLON_REGULAR\n", offset);
        break;
    case GRIDHTYPE_LATLON_IRREGULAR:
        enkf_printf("%s  type = LATLON_IRREGULAR\n", offset);
        break;
#if !defined(NO_GRIDUTILS)
    case GRIDHTYPE_CURVILINEAR:
        enkf_printf("%s  type = CURVILINEAR\n", offset);
        break;
#endif
    default:
        enkf_printf("%s  type = NONE\n", offset);
    }
    enkf_printf("%s  periodic by X = %s\n", offset, grid_isperiodic_x(g) ? "yes" : "no");
    enkf_printf("%s  periodic by Y = %s\n", offset, grid_isperiodic_y(g) ? "yes" : "no");
    grid_getdims(g, &nx, &ny, &nz);
    enkf_printf("%s  dims = %d x %d x %d\n", offset, nx, ny, nz);
    if (g->lontype == LONTYPE_180)
        enkf_printf("%s  longitude range = [-180, 180]\n", offset);
    else if (g->lontype == LONTYPE_360)
        enkf_printf("%s  longitude range = [0, 360]\n", offset);
    else if (g->lontype == LONTYPE_NONE)
        enkf_printf("%s  longitude range = any\n", offset);
}

/**
 */
void grid_describeprm(void)
{
    enkf_printf("\n");
    enkf_printf("  Grid parameter file format for z-model:\n");
    enkf_printf("\n");
    enkf_printf("    NAME             = <name>\n");
    enkf_printf("    DATA             = <data file name>\n");
    enkf_printf("    XDIMNAME         = <x dimension name>\n");
    enkf_printf("    YDIMNAME         = <y dimension name>\n");
    enkf_printf("    ZDIMNAME         = <z dimension name>\n");
    enkf_printf("    XVARNAME         = <x variable name>\n");
    enkf_printf("    YVARNAME         = <y variable name>\n");
    enkf_printf("    ZVARNAME         = <z variable name>\n");
    enkf_printf("    DEPTHVARNAME     = <depth variable name>\n");
    enkf_printf("    NUMLEVELSVARNAME = <# of levels variable name>\n");
    enkf_printf("\n");
    enkf_printf("  Notes:\n");
    enkf_printf("    1. < ... > denotes a description of an entry\n");
    enkf_printf("\n");
}

/**
 */
void grid_getdims(grid* g, int* ni, int* nj, int* nk)
{
    if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR) {
        gnxy* nodes = (gnxy *) g->gridnodes;

        *ni = nodes->nx;
        *nj = nodes->ny;
        if (nk != NULL)
            *nk = nodes->nz;
#if !defined(NO_GRIDUTILS)
    } else if (g->htype == GRIDHTYPE_CURVILINEAR) {
        gnc* nodes = (gnc *) g->gridnodes;

        *ni = gridnodes_getnx(nodes->gn);
        *nj = gridnodes_getny(nodes->gn);
        if (nk != NULL)
            *nk = nodes->nz;
#endif
    } else
        enkf_quit("programming error");
}

/**
 */
static double z2fk(int n, double* zt, double* zb, double z)
{
    int ascending, i1, i2, imid;

    ascending = (zt[n - 1] > zt[0]) ? 1 : 0;

    if ((ascending && z < zt[0]) || (!ascending && z > zt[0]))
        return 0.0;
    if ((ascending && z > zt[n - 1]) || (!ascending && z < zt[n - 1]))
        return (double) (n - 1);

    i1 = 0;
    i2 = n - 1;
    if (ascending) {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (z > zb[imid])
                i1 = imid;
            else
                i2 = imid;
        }
    } else {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (z > zb[imid])
                i2 = imid;
            else
                i1 = imid;
        }
    }

    if (z < zb[i1 + 1])
        return (double) i1 + (z - zt[i1]) / (zb[i1 + 1] - zb[i1]);
    else
        return (double) i1 + 0.5 + (z - zb[i1 + 1]) / (zb[i1 + 2] - zb[i1]);
}

/**
 */
static double x2fi_irreg(int n, double* v, double* vb, double x, int periodic)
{
    int ascending, i1, i2, imid;

    if (n < 2)
        return NaN;

    if (periodic) {
        if (x < 0.0)
            x = x + 360.0;
        else if (x >= 360.0)
            x = x - 360.0;
    }

    ascending = (v[n - 1] > v[0]) ? 1 : 0;

    if ((ascending && x < v[0]) || (!ascending && x > v[0])) {
        double fi = (x - v[0]) / (v[1] - v[0]);

        if (periodic) {
            assert(fi >= -0.5);
            return fi;
        }
        return (fi >= -0.5) ? 0.0 : NaN;
    }
    if ((ascending && x > v[n - 1]) || (!ascending && x < v[n - 1])) {
        double fi = (double) (n - 1) + (x - v[n - 1]) / (v[n - 1] - v[n - 2]);

        if (periodic) {
            assert(fi < (double) n - 0.5);
            return fi;
        }
        return (fi <= (double) n - 0.5) ? (double) (n - 1) : NaN;
    }

    i1 = 0;
    i2 = n - 1;
    if (ascending) {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (x > vb[imid])
                i1 = imid;
            else
                i2 = imid;
        }
        if (x < vb[i1 + 1])
            return (double) i1 + (x - v[i1]) / (vb[i1 + 1] - vb[i1]);
        else
            return (double) i1 + 0.5 + (x - vb[i1 + 1]) / (vb[i1 + 2] - vb[i1]);
    } else {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (x > vb[imid])
                i2 = imid;
            else
                i1 = imid;
        }
        if (x > vb[i1 + 1])
            return (double) i1 + (x - v[i1]) / (vb[i1 + 1] - vb[i1]);
        else
            return (double) i1 + 0.5 + (x - vb[i1 + 1]) / (vb[i1 + 2] - vb[i1]);
    }
}

/**
 */
static double x2fi_reg(int n, double* v, double x, int periodic)
{
    double fi;

    if (n < 2)
        return NaN;

    if (periodic) {
        if (x < 0.0)
            x = x + 360.0;
        else if (x >= 360.0)
            x = x - 360.0;
    }

    fi = (x - v[0]) / (v[n - 1] - v[0]) * (double) (n - 1);

    if (fi < -0.5 || fi > (double) n - 0.5)
        return NaN;

    if (!periodic) {
        if (fi < 0.0)
            fi = 0.0;
        else if (fi > (double) (n - 1))
            fi = (double) (n - 1);
    }

    return fi;
}

/**
 */
static double fi2x(int n, double* v, double fi, int periodic)
{
    double ifrac;
    int i;

    if (n < 2)
        return NaN;

    if (fi < -1.0 || fi > (double) n)
        return NaN;

    ifrac = (int) (fi - floor(fi));

    if (fi < 0.0 || fi >= (double) (n - 1)) {
        double v1, v2, x;

        if (!periodic)
            return NaN;

        v1 = v[n - 1];
        v2 = v[0];
        if (v1 - v2 > 180.0)
            v2 += 360.0;
        x = v1 + ifrac * (v2 - v1);

        return (x < 360.0) ? x : x - 360.0;
    }

    i = (int) fi;
    assert(i < n - 1);

    return v[i] + ifrac * (v[i + 1] - v[i]);
}

/**
 */
static void g1_xy2fij(void* p, double x, double y, double* fi, double* fj)
{
    gnxy* nodes = (gnxy *) ((grid*) p)->gridnodes;

    *fi = x2fi_reg(nodes->nx, nodes->x, x, nodes->periodic_x);
    *fj = x2fi_reg(nodes->ny, nodes->y, y, nodes->periodic_y);
}

/**
 */
static void g12_fij2xy(void* p, double fi, double fj, double* x, double* y)
{
    gnxy* nodes = (gnxy *) ((grid*) p)->gridnodes;

    *x = fi2x(nodes->nx, nodes->x, fi, nodes->periodic_x);
    *y = fi2x(nodes->ny, nodes->y, fj, nodes->periodic_y);
}

/**
 */
static void g2_xy2fij(void* p, double x, double y, double* fi, double* fj)
{
    gnxy* nodes = (gnxy *) ((grid*) p)->gridnodes;

    *fi = x2fi_irreg(nodes->nx, nodes->x, nodes->xb, x, nodes->periodic_x);
    *fj = x2fi_irreg(nodes->ny, nodes->y, nodes->yb, y, nodes->periodic_y);
}

#if !defined(NO_GRIDUTILS)
/**
 */
static void gc_xy2fij(void* p, double x, double y, double* fi, double* fj)
{
    gnc* nodes = (gnc *) ((grid*) p)->gridnodes;

    gridmap_xy2fij(nodes->gm, x, y, fi, fj);
}

/**
 */
static void gc_fij2xy(void* p, double fi, double fj, double* x, double* y)
{
    gnc* nodes = (gnc *) ((grid*) p)->gridnodes;

    gridmap_fij2xy(nodes->gm, fi, fj, x, y);
}
#endif

/**
 */
static void z2fk_1d(void* p, double fi, double fj, double z, double* fk)
{
    gnxy* nodes = (gnxy *) ((grid*) p)->gridnodes;

    *fk = (z == 0.0) ? 0.0 : z2fk(nodes->nz, nodes->z, nodes->zb, z);
}

/**
 */
static void grid_setlontype(grid* g)
{
    double xmin = DBL_MAX;
    double xmax = -DBL_MAX;

    if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR) {
        double* x = ((gnxy *) g->gridnodes)->x;
        int nx = ((gnxy *) g->gridnodes)->nx;

        if (xmin < x[0])
            xmin = x[0];
        if (xmin < x[nx - 1])
            xmin = x[nx - 1];
        if (xmax > x[0])
            xmax = x[0];
        if (xmax > x[nx - 1])
            xmax = x[nx - 1];
#if !defined(NO_GRIDUTILS)
    } else if (g->htype == GRIDHTYPE_CURVILINEAR) {
        double** x = gridnodes_getx(((gnc *) g->gridnodes)->gn);
        int nx = gridnodes_getnce1(((gnc *) g->gridnodes)->gn);
        int ny = gridnodes_getnce2(((gnc *) g->gridnodes)->gn);
        int i, j;

        for (j = 0; j < ny; ++j) {
            for (i = 0; i < nx; ++i) {
                if (xmin < x[j][i])
                    xmin = x[j][i];
                if (xmax > x[j][i])
                    xmax = x[j][i];
            }
        }
#endif
    }
    if (xmin < 0.0 && xmax <= 180.0)
        g->lontype = LONTYPE_180;
    else if (xmin >= 0 && xmax <= 360.0)
        g->lontype = LONTYPE_360;
}

/**
 */
void grid_setcoords(grid* g, int htype, int periodic_x, int periodic_y, int nx, int ny, int nz, void* x, void* y, double* z)
{
    g->htype = htype;
    if (htype == GRIDHTYPE_LATLON_REGULAR) {
        g->xy2fij_fn = g1_xy2fij;
        g->fij2xy_fn = g12_fij2xy;
        g->gridnodes = gnxy_create(nx, ny, nz, x, y, z, periodic_x, periodic_y, 1);
    } else if (htype == GRIDHTYPE_LATLON_IRREGULAR) {
        g->xy2fij_fn = g2_xy2fij;
        g->fij2xy_fn = g12_fij2xy;
        g->gridnodes = gnxy_create(nx, ny, nz, x, y, z, periodic_x, periodic_y, 0);
#if !defined(NO_GRIDUTILS)
    } else if (htype == GRIDHTYPE_CURVILINEAR) {
        g->xy2fij_fn = gc_xy2fij;
        g->fij2xy_fn = gc_fij2xy;
        g->gridnodes = gnc_create(nx, ny, nz, x, y, z);
#endif
    } else
        enkf_quit("programming error");

    grid_setlontype(g);
    g->z2fk_fn = z2fk_1d;
}

/**
 */
void grid_setdepth(grid* g, float** depth)
{
    g->depth = depth;
}

/**
 */
void grid_setnumlevels(grid* g, int** numlevels)
{
    g->numlevels = numlevels;
}

/**
 */
char* grid_getname(grid* g)
{
    return g->name;
}

/**
 */
int grid_gethtype(grid* g)
{
    return g->htype;
}

/**
 */
float** grid_getdepth(grid* g)
{
    return g->depth;
}

/**
 */
int** grid_getnumlevels(grid* g)
{
    return g->numlevels;
}

/**
 */
int grid_getlontype(grid* g)
{
    return g->lontype;
}

/**
 */
grid_xy2fij_fn grid_getxy2fijfn(grid* g)
{
    return g->xy2fij_fn;
}

/**
 */
grid_z2fk_fn grid_getz2fkfn(grid* g)
{
    return g->z2fk_fn;
}

/**
 */
grid_fij2xy_fn grid_getfij2xyfn(grid* g)
{
    return g->fij2xy_fn;
}

/**
 */
int grid_isperiodic_x(grid* g)
{
    if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR) {
        gnxy* nodes = (gnxy *) ((grid*) g)->gridnodes;

        return nodes->periodic_x;
    }

    return 0;
}

/**
 */
int grid_isperiodic_y(grid* g)
{
    if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR) {
        gnxy* nodes = (gnxy *) ((grid*) g)->gridnodes;

        return nodes->periodic_y;
    }

    return 0;
}

/**
 */
void grid_settocartesian_fn(grid* g, grid_tocartesian_fn fn)
{
    g->tocartesian_fn = fn;
}

/**
 */
void grid_tocartesian(grid* g, double* in, double* out)
{
    g->tocartesian_fn(in, out);
}
