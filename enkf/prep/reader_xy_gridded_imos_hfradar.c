/******************************************************************************
 *
 * File:        reader_xy_gridded_hfradar.c
 *
 * Created:     08/03/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *              Hugo Oliveira
 *              AODN
 *
 * Description: A Generic reader for gridded surface velocity observations from HF-Radar.
 *                It currently assumes the following:
 *              - There is a rotation mode set by the user.
 *              - there is only one data record (2D field);
 *              - longitude is the inner ("fast") coordinate of the variable.
 *                There are a number of parameters that must (++) or can be
 *              specified if they differ from the default value (+). Some
 *              parameters are optional (-):
 *              - VARNAME (++)
 *              - TIMENAME ("time") (+)
 *              - ROTATION (True 1 | False 0)  (+)
 *                  a boolean to rotate the vectors to the current grid when reading.
 *                If rotation is allowed, several [U_,V_][VARNAME] variables can be provided (see below). A variable angle in the grid.prm file is also required [ANGLE]. Otherwise, standard obs.prm parameters are used.
 *              - U_VARNAME ("UCUR") (-)
 *                  Name of the variable for the zonal velocity. This will only be used if rotation is True. The default follows IMOS conventions.
 *              - V_VARNAME ("VCUR") (-)
 *                  Name of the variable for the meridional velocity. same as above.
 *              - U_VARSHIFT (-)
 *                  u data offset to be added to u values only. This will only be used if rotation is True. The default is 0 and it's performed together with the standard "VARSHIFT" parameter.
 *              - V_VARSHIFT (-)
 *                 v data offset to be added to v values. same as above.
 *              -----------------------------------------------------------------
 *              - NPOINTSNAME ("npoints") (-)
 *                  number of collated points for each datum; used basically as
 *                  a data mask n = 0
 *              - LONNAME ("lon" | "longitude") (+)
 *              - LATNAME ("lat" | "latitude") (+)
 *              - STDNAME ("std") (-)
 *              - ESTDNAME ("error_std") (-)
 *                  error STD; if absent then needs to be specified externally
 *                  in the observation data parameter file
 *                  If both STDNAME or ESTDNAME are absent, the code will try to
 *                  read a VARNAME.'error_std' attribute value.
 *              - VARSHIFT (-)
 *                  data offset to be added
 *              - MINDEPTH (-)
 *                  minimal allowed depth
 *              - INSTRUMENT (-)
 *                  instrument string that will be used for calculating
 *                  instrument stats
 *              - QCFLAGNAME (-)
 *                  name of the QC flag variable, 0 <= qcflag <= 31
 *              - QCFLAGVALS (-)
 *                  the list of allowed values of QC flag variable
 *              Note: it is possible to have multiple entries of QCFLAGNAME and
 *                QCFLAGVALS combination, e.g.:
 *                  PARAMETER QCFLAGNAME = TEMP_quality_control
 *                  PARAMETER QCFLAGVALS = 1
 *                  PARAMETER QCFLAGNAME = DEPTH_quality_control
 *                  PARAMETER QCFLAGVALS = 1
 *                  PARAMETER QCFLAGNAME = LONGITUDE_quality_control
 *                  PARAMETER QCFLAGVALS = 1,8
 *                  PARAMETER QCFLAGNAME = LATITUDE_quality_control
 *                  PARAMETER QCFLAGVALS = 1,8
 *                An observation is considered valid if each of the specified
 *                flags takes a permitted value.
 *
 * Revisions:    HO 5/5/2018
 *                Add angle rotation logic and new entries.
 *               PS 6/7/2018
 *                Added parameters QCFLAGNAME and QCFLAGVALS. The latter is
 *                supposed to contain a list of allowed flag values.
 *               HO 10/1/2018
 *                Merged new code with rotation logics.
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

#define TYPE_DOUBLE 0
#define TYPE_SHORT 1

/**
 */
void reader_xy_gridded_imos_hfradar(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);

    /* Variables required to interpolate angles */
    int isperiodic_i = grid_isperiodic_i(g);
    int** numlevels = grid_getnumlevels(g);
    float** grid_angle = grid_getangle(g);
    int grid_ni = 0, grid_nj = 0;

    char* varname = NULL;

    /* U/V vars */
    char* u_varname = NULL;
    char* v_varname = NULL;

    char* lonname = NULL;
    char* latname = NULL;
    char* npointsname = NULL;
    char* stdname = NULL;
    char* estdname = NULL;
    char* timename = NULL;
    int ndim_var, ndim_xy;
    size_t dimlen_var[3], dimlen_xy[2];
    int nqcflags = 0;
    char** qcflagname = NULL;

    /* UV rotation*/
    char* rotation_flag = NULL;
    int need_rotation = -1;
    uint32_t* qcflagvals = 0;

    int ncid;
    float varshift = 0.0;

    /* UV varshifts */
    float u_varshift = 0.0;
    float v_varshift = 0.0;

    double mindepth = 0.0;
    char instrument[MAXSTRLEN];

    int iscurv = -1;
    size_t ni = 0, nj = 0, n = 0, n_var = 0;
    int varid_lon = -1, varid_lat = -1;
    double* lon = NULL;
    double* lat = NULL;
    int varid_var = -1, varid_npoints = -1, varid_std = -1, varid_estd = -1, varid_time = -1;
    float* var = NULL;
    float var_fill_value = NAN;
    float var_add_offset = NAN, var_scale_factor = NAN;
    double var_estd = NAN;
    short* npoints = NULL;
    float* std = NULL;
    float std_add_offset = NAN, std_scale_factor = NAN;
    float std_fill_value = NAN;
    float* estd = NULL;
    float estd_add_offset = NAN, estd_scale_factor = NAN;
    float estd_fill_value = NAN;
    /* Assumes any type instead of unsigned 32 int for qcflag */
    int** qcflag = NULL;


    int have_time = 1;
    int singletime = -1;
    float* time = NULL;
    float time_add_offset = NAN, time_scale_factor = NAN;
    float time_fill_value = NAN;
    char tunits[MAXSTRLEN];
    double tunits_multiple = NAN, tunits_offset = NAN;
    int i, nobs_read;

    /* UV requirements */
    int varid_u = -1, varid_v =-1;
    float* u_var = NULL;
    float* v_var = NULL;
    float u_var_fill_value = NAN;
    float u_var_add_offset = NAN, u_var_scale_factor = NAN;
    float v_var_fill_value = NAN;
    float v_var_add_offset = NAN, v_var_scale_factor = NAN;

    strcpy(instrument, meta->product);
    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "VARNAME") == 0)
            varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "TIMENAME") == 0)
            timename = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "NPOINTSNAME") == 0)
            npointsname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LONNAME") == 0)
            lonname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LATNAME") == 0)
            latname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "STDNAME") == 0)
            stdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ESTDNAME") == 0)
            estdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "VARSHIFT") == 0) {
            if (!str2float(meta->pars[i].value, &varshift))
                enkf_quit("%s: can not convert VARSHIFT = \"%s\" to float\n", meta->prmfname, meta->pars[i].value);
            enkf_printf("        VARSHIFT = %s\n", meta->pars[i].value);
        } else if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
            if (!str2double(meta->pars[i].value, &mindepth))
                enkf_quit("%s: can not convert MINDEPTH = \"%s\" to double\n", meta->prmfname, meta->pars[i].value);
            enkf_printf("        MINDEPTH = %f\n", mindepth);
        } else if (strcasecmp(meta->pars[i].name, "INSTRUMENT") == 0)
            strncpy(instrument, meta->pars[i].value, MAXSTRLEN);
        else if (strcasecmp(meta->pars[i].name, "QCFLAGNAME") == 0 || strcasecmp(meta->pars[i].name, "QCFLAGVALS") == 0)
            /*
             * QCFLAGNAME and QCFLAGVALS are dealt with separately
             */
            ;
        /* UV, rotation, requirements */
        else if (strcasecmp(meta->pars[i].name, "ROTATION") == 0) 
            rotation_flag = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "U_VARNAME") == 0)
        u_varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "V_VARNAME") == 0)
        v_varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "U_VARSHIFT") == 0) {
            if (!str2float(meta->pars[i].value, &u_varshift))
                enkf_quit("%s: can not convert U_VARSHIFT = \"%s\" to float\n", meta->prmfname, meta->pars[i].value);
            else
                enkf_printf("        U_VARSHIFT = %s\n", meta->pars[i].value);
        } else if (strcasecmp(meta->pars[i].name, "V_VARSHIFT") == 0) {
            if (!str2float(meta->pars[i].value, &v_varshift))
                enkf_quit("%s: can not convert V_VARSHIFT = \"%s\" to float\n", meta->prmfname, meta->pars[i].value);
            else
                enkf_printf("        V_VARSHIFT = %s\n", meta->pars[i].value);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    get_qcflags(meta, &nqcflags, &qcflagname, &qcflagvals);

    if (varname == NULL)
        enkf_quit("reader_xy_gridded_hfradar(): %s: VARNAME not specified", fname);
    else
        enkf_printf("        VARNAME = %s\n", varname);

    /* Check rotation parameter */
    if (rotation_flag != NULL) {
        int bool_rotation_flag, valid_rotation_flag;
        bool_rotation_flag = atoi(rotation_flag);
        valid_rotation_flag = (bool_rotation_flag == 0 || bool_rotation_flag == 1);

        if (valid_rotation_flag) {
            need_rotation = bool_rotation_flag == 1;
            grid_getdims(g, &grid_ni, &grid_nj, NULL);
            if (need_rotation == 1 && grid_angle == NULL)
                enkf_quit("reader_xy_gridded_hfradar(): %s cannot rotate vectors without ANGLE parameter defined within grid prm file",fname);
            else
                enkf_printf("\t#Rotation of vectors enabled\n");
        }
        else
            enkf_quit("reader_xy_gridded_hfradar(): %s invalid ROTATION parameter. Need to be a boolean digit");
    }

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid_var);
    ncw_inq_vardims(ncid, varid_var, 3, &ndim_var, dimlen_var);
    if (ndim_var == 3) {
        if (!ncw_var_hasunlimdim(ncid, varid_var))
            enkf_quit("reader_xy_gridded_hfradar(): %s: %s: depends on 3 dimensions, but has no unlimited dimension", fname, varname);
        if (dimlen_var[0] != 1)
            enkf_quit("reader_xy_gridded_hfradar(): %s: %s: %d records (currently only one is allowed)", fname, varname, dimlen_var[0]);
        n_var = dimlen_var[1] * dimlen_var[2];
    } else if (ndim_var == 2) {
        if (ncw_var_hasunlimdim(ncid, varid_var))
            enkf_quit("reader_xy_gridded_hfradar(): %s: %s: not enough spatial dimensions (must be 2)", fname, varname);
        n_var = dimlen_var[0] * dimlen_var[1];
    } else if (ndim_var != 2)
        enkf_quit("reader_xy_gridded_hfradar(): %s: %s: %d dimensions (must be 2 or 3 with only one record)", fname, varname, ndim_var);

    lonname = get_lonname(ncid, lonname);
    if (lonname != NULL) {
        enkf_printf("        LONNAME = %s\n", lonname);
        ncw_inq_varid(ncid, lonname, &varid_lon);
    } else {
        if (ncw_var_exists(ncid, "lon"))
            ncw_inq_varid(ncid, "lon", &varid_lon);
        else if (ncw_var_exists(ncid, "longitude"))
            ncw_inq_varid(ncid, "longitude", &varid_lon);
        else if (ncw_var_exists(ncid, "LONGITUDE"))
            ncw_inq_varid(ncid, "LONGITUDE", &varid_lon);
        else
            enkf_quit("reader_xy_gridded_hfradar(): %s: could not find longitude variable", fname);
    }

    ncw_inq_vardims(ncid, varid_lon, 2, &ndim_xy, dimlen_xy);
    if (ndim_xy == 1) {
        iscurv = 0;
        ni = dimlen_xy[0];
    } else if (ndim_xy == 2) {
        iscurv = 1;
        ni = dimlen_xy[1];
        nj = dimlen_xy[0];
    } else
        enkf_quit("reader_xy_gridded_hfradar(): %s: coordinate variable \"%s\" has neither 1 or 2 dimensions", fname, lonname);

    latname = get_latname(ncid, latname);
    if (latname != NULL) {
        enkf_printf("        LATNAME = %s\n", latname);
        ncw_inq_varid(ncid, latname, &varid_lat);
    } else {
        if (ncw_var_exists(ncid, "lat"))
            ncw_inq_varid(ncid, "lat", &varid_lat);
        else if (ncw_var_exists(ncid, "latitude"))
            ncw_inq_varid(ncid, "latitude", &varid_lat);
        else if (ncw_var_exists(ncid, "LATITUDE"))
            ncw_inq_varid(ncid, "LATITUDE", &varid_lat);
        else
            enkf_quit("reader_xy_gridded_hfradar(): %s: could not find latitude variable", fname);
    }
    if (ndim_xy == 1)
        ncw_inq_vardims(ncid, varid_lat, 1, &ndim_xy, &nj);

    if (iscurv == 0) {
        ncw_check_varndims(ncid, varid_lat, 1);
        ncw_inq_vardims(ncid, varid_lat, 1, NULL, &nj);
    } else
        ncw_check_vardims(ncid, varid_lat, 2, dimlen_xy);

    enkf_printf("        (ni, nj) = (%u, %u)\n", ni, nj);
    n = ni * nj;
    if (n != n_var)
        enkf_quit("reader_xy_gridded_hfradar(): %s: dimensions of variable \"%s\" do not match coordinate dimensions", fname, varname);
    if (dimlen_var[ndim_var - 1] != ni)
        enkf_quit("reader_xy_gridded_hfradar(): %s: %s: longitude must be the inner coordinate", fname, varname);

    if (iscurv == 0) {
        lon = malloc(ni * sizeof(double));
        lat = malloc(nj * sizeof(double));
    } else {
        lon = malloc(n * sizeof(double));
        lat = malloc(n * sizeof(double));
    }
    ncw_get_var_double(ncid, varid_lon, lon);
    ncw_get_var_double(ncid, varid_lat, lat);

    var = malloc(n * sizeof(float));
    ncw_get_var_float(ncid, varid_var, var);
    if (ncw_att_exists(ncid, varid_var, "add_offset")) {
        ncw_get_att_float(ncid, varid_var, "add_offset", &var_add_offset);
        ncw_get_att_float(ncid, varid_var, "scale_factor", &var_scale_factor);
    }
    if (ncw_att_exists(ncid, varid_var, "_FillValue"))
        ncw_get_att_float(ncid, varid_var, "_FillValue", &var_fill_value);

    if (npointsname != NULL)
        ncw_inq_varid(ncid, npointsname, &varid_npoints);
    else if (ncw_var_exists(ncid, "npoints"))
        ncw_inq_varid(ncid, "npoints", &varid_npoints);
    if (varid_npoints >= 0) {
        npoints = malloc(n * sizeof(short));
        ncw_get_var_short(ncid, varid_npoints, npoints);
    }

    if (stdname != NULL)
        ncw_inq_varid(ncid, stdname, &varid_std);
    else if (ncw_var_exists(ncid, "std"))
        ncw_inq_varid(ncid, "std", &varid_std);
    if (varid_std >= 0) {
        std = malloc(n * sizeof(float));
        ncw_get_var_float(ncid, varid_std, std);
        if (ncw_att_exists(ncid, varid_std, "_FillValue"))
            ncw_get_att_float(ncid, varid_std, "_FillValue", &std_fill_value);
        if (ncw_att_exists(ncid, varid_std, "add_offset")) {
            ncw_get_att_float(ncid, varid_std, "add_offset", &std_add_offset);
            ncw_get_att_float(ncid, varid_std, "scale_factor", &std_scale_factor);
        }
    }

    if (estdname != NULL)
        ncw_inq_varid(ncid, estdname, &varid_estd);
    else if (ncw_var_exists(ncid, "error_std"))
        ncw_inq_varid(ncid, "error_std", &varid_estd);
    if (varid_estd >= 0) {
        estd = malloc(n * sizeof(float));
        ncw_get_var_float(ncid, varid_estd, estd);
        if (ncw_att_exists(ncid, varid_estd, "_FillValue"))
            ncw_get_att_float(ncid, varid_estd, "_FillValue", &estd_fill_value);
        if (ncw_att_exists(ncid, varid_estd, "add_offset")) {
            ncw_get_att_float(ncid, varid_estd, "add_offset", &estd_add_offset);
            ncw_get_att_float(ncid, varid_estd, "scale_factor", &estd_scale_factor);
        }
    }

    if (std == NULL && estd == NULL)
        if (ncw_att_exists(ncid, varid_var, "error_std")) {
            ncw_check_attlen(ncid, varid_var, "error_std", 1);
            ncw_get_att_double(ncid, varid_var, "error_std", &var_estd);
        }

    if (nqcflags > 0) {
        int varid = -1;
        qcflag = alloc2d(nqcflags, n, sizeof(int));
        for (i = 0; i < nqcflags; ++i) {
                ncw_inq_varid(ncid, qcflagname[i], &varid);
                /* Change Calling all qcflagvals as unsigned integers can be a problem.
                 * type conversion somewhat failed when reading int8 with uint32 in netCDF.
                 * Assumes that all qcflagvals are of int */
                ncw_get_var_int(ncid, varid, qcflag[i]);
        }
    }

    timename = get_timename(ncid, timename);
    if (timename != NULL) {
        enkf_printf("        TIMENAME = %s\n", timename);
        ncw_inq_varid(ncid, timename, &varid_time);
    } else {
        enkf_printf("        reader_xy_gridded_hfradar(): %s: no TIME variable\n", fname);
        have_time = 0;
    }

    if (have_time) {
        int timendims;
        int timedimids[NC_MAX_DIMS];
        size_t timelen = 1;

        ncw_inq_varndims(ncid, varid_time, &timendims);
        ncw_inq_vardimid(ncid, varid_time, timedimids);
        for (i = 0; i < timendims; ++i) {
            size_t dimlen;

            ncw_inq_dimlen(ncid, timedimids[i], &dimlen);
            timelen *= dimlen;
        }

        if (timelen == 1) {
            singletime = 1;
            time = malloc(sizeof(float));
        } else {
            singletime = 0;
            assert(timelen == n);
            time = malloc(n * sizeof(float));
        }

        ncw_get_var_float(ncid, varid_time, time);
        if (ncw_att_exists(ncid, varid_time, "_FillValue"))
            ncw_get_att_float(ncid, varid_time, "_FillValue", &time_fill_value);
        if (ncw_att_exists(ncid, varid_time, "add_offset")) {
            ncw_get_att_float(ncid, varid_time, "add_offset", &time_add_offset);
            ncw_get_att_float(ncid, varid_time, "scale_factor", &time_scale_factor);
        }
        ncw_get_att_text(ncid, varid_time, "units", tunits);
        tunits_convert(tunits, &tunits_multiple, &tunits_offset);
    }

    /* Pre-load U/V if rotation is required*/
    if (need_rotation == 1) {

        if (u_varname != NULL)
            ncw_inq_varid(ncid, u_varname, &varid_u);
        else if (ncw_var_exists(ncid, "UCUR")) {
            ncw_inq_varid(ncid, "UCUR", &varid_u);
            u_varname = "UCUR";
        } else
            enkf_quit("reader_xy_gridded_hfradar(): %s: Variable %s not found.",fname, u_varname);

        u_var = malloc(n * sizeof(float));
        ncw_get_var_float(ncid,varid_u,u_var);
        if (ncw_att_exists(ncid,varid_u, "_FillValue"))
            ncw_get_att_float(ncid, varid_u, "_FillValue", &u_var_fill_value);
        if (ncw_att_exists(ncid,varid_u,"add_offset")) {
            ncw_get_att_float(ncid, varid_u, "add_offset", &u_var_add_offset);
            ncw_get_att_float(ncid, varid_u, "scale_factor", &u_var_scale_factor);
            }

        if (v_varname != NULL)
            ncw_inq_varid(ncid, v_varname, &varid_v);
        else if (ncw_var_exists(ncid, "VCUR")) {
            ncw_inq_varid(ncid, "VCUR", &varid_v);
            v_varname = "VCUR";
        } else
            enkf_quit("reader_xy_gridded_hfradar(): %s: Variable %s not found.",fname, v_varname);

        v_var = malloc(n * sizeof(float));
        ncw_get_var_float(ncid,varid_v,v_var);
        if (ncw_att_exists(ncid,varid_v, "_FillValue"))
            ncw_get_att_float(ncid,varid_v, "_FillValue", &v_var_fill_value);
        if (ncw_att_exists(ncid,varid_v,"add_offset")) {
            ncw_get_att_float(ncid, varid_v, "add_offset", &v_var_add_offset);
            ncw_get_att_float(ncid, varid_v, "scale_factor", &v_var_scale_factor);
            }
    }

    ncw_close(ncid);

    nobs_read = 0;
    for (i = 0; i < (int) n; ++i) {
        observation* o;
        obstype* ot;
        int ii;
        int invalid_npoints, invalid_var, invalid_std, invalid_estd, invalid_time;
        int invalid_u_var, invalid_v_var;
        
        invalid_npoints = (npoints != NULL && npoints[i] == 0);
        invalid_var = var[i] == var_fill_value || isnan(var[i]);
        invalid_std = (std != NULL && (std[i] == std_fill_value || isnan(std[i])));
        invalid_estd = (estd != NULL && (estd[i] == estd_fill_value || isnan(estd[i])));
        invalid_time = (have_time && !singletime && (time[i] == time_fill_value || isnan(time[i])));
        invalid_u_var = (u_var != NULL && (u_var[i] == u_var_fill_value || isnan(u_var[i])));
        invalid_v_var = (v_var != NULL && (v_var[i] == v_var_fill_value || isnan(v_var[i])));

        if (invalid_npoints || invalid_var || invalid_std || invalid_estd || invalid_time || invalid_u_var || invalid_v_var)
            continue;

        int invalid_qc = 1;
        char bitflag;
        for (ii = 0; ii < nqcflags; ++ii) {
            bitflag = 0;
            bitflag |= 1 << qcflag[ii][i];
            if ((bitflag & qcflagvals[ii])) {
                invalid_qc = 0;
                break;
            }
        }
        if (invalid_qc)
            continue;

        nobs_read++;
        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        ot = &obs->obstypes[o->type];
        o->instrument = st_add_ifabsent(obs->instruments, instrument, -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = 0;
        if (iscurv == 0) {
            o->lon = lon[i % ni];
            o->lat = lat[i % nj];
        } else {
            o->lon = lon[i];
            o->lat = lat[i];
        }
        o->depth = 0.0;
        o->fk = (double) ksurf;
        o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        if (!obs->allobs && o->status == STATUS_LAND)
            continue;
        if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax))
            o->status = STATUS_OUTSIDEOBSDOMAIN;
        o->model_depth = NAN;   /* set in obs_add() */
        if (have_time) {
            double t = (singletime) ? time[0] : time[i];

            if (!isnan(time_add_offset))
                o->date = (double) (t * time_scale_factor + time_add_offset) * tunits_multiple + tunits_offset;
            else
                o->date = (double) t* tunits_multiple + tunits_offset;
        } else
            o->date = NAN;

        /* Now set values since we got fi,fj */
        if (need_rotation == 1) {
            double zonal, meridional, oangle;
            if (!isnan(u_var_add_offset))
                zonal = u_var[i] * u_var_scale_factor + u_var_add_offset + varshift + u_varshift;
            else
                zonal = u_var[i];

            if (!isnan(v_var_add_offset))
                meridional = v_var[i] * v_var_scale_factor + v_var_add_offset + varshift + u_varshift;
            else
                meridional = v_var[i];

            oangle = interpolate2d(o->fi, o-> fj, grid_ni, grid_nj, grid_angle, numlevels, isperiodic_i);

            if (strcmp(varname, u_varname) == 0)
                o->value = zonal*cos(oangle) - meridional*sin(oangle);
            if (strcmp(varname, v_varname) == 0)
                o->value = zonal*sin(oangle) + meridional*cos(oangle);
        }
        else {
            if (!isnan(var_add_offset))
                o->value = (double) (var[i] * var_scale_factor + var_add_offset + varshift);
            else
                o->value = (double) (var[i] + varshift);
        }

        if (estd == NULL && std == NULL){
            if (!isnan(var_estd))
                o->std = var_estd;
        }
        else {
            if (std == NULL)
                o->std = 0.0;
            else {
                if (!isnan(std_add_offset))
                    o->std = (double) (std[i] * std_scale_factor + std_add_offset);
                else
                    o->std = (double) std[i];
            }
            if (estd != NULL) {
                if (!isnan(estd_add_offset)) {
                    double std2 = (double) (estd[i] * estd_scale_factor + estd_add_offset);
                    o->std = (o->std > std2) ? o->std : std2;
                } else
                    o->std = (o->std > estd[i]) ? o->std : estd[i];
            }
        }
        o->aux = -1;

        obs->nobs++;
    }
    enkf_printf("        nobs = %d\n", nobs_read);

    free(lon);
    free(lat);
    free(var);
    if (std != NULL)
        free(std);
    if (estd != NULL)
        free(estd);
    if (npoints != NULL)
        free(npoints);
    if (time != NULL)
        free(time);
    if (nqcflags > 0) {
        free(qcflagname);
        free(qcflagvals);
        free(qcflag);
    }
    /* Free rotation related vars */
    if (u_var != NULL)
        free(u_var);
    if (v_var != NULL)
        free(v_var);
}
