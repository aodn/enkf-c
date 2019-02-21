/******************************************************************************
 *
 * File:        reader_z_xbt_profile.c        
 *
 * Created:     06/02/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *              Hugo Oliveira
 *              AODN
 *
 * Description: Generic reader for soop xbt data.
 *                There are a number of parameters that must (++) or can be
 *              specified if they differ from the default value (+). Some
 *              parameters are optional (-):
 *              - VARNAME (++)
 *              - TIMENAME ("time") (+)
 *              - NPOINTSNAME ("npoints") (-)
 *                  number of collated points for each datum; used basically as
 *                  a data mask n = 0
 *              - LONNAME ("lon" | "longitude") (+)
 *              - LATNAME ("lat" | "latitude") (+)
 *              - ZNAME ("z") (+)
 *              - STDNAME ("std") (-)
 *                  internal variability of the collated data
 *              - ESTDNAME ("error_std") (-)
 *                  error STD; if absent then needs to be specified externally
 *                  in the oobservation data parameter file
 *              - VARSHIFT (-)
 *                  data offset to be added
 *              - MINDEPTH (-)
 *                  minimal allowed depth
 *              - INSTRUMENT (-)
 *                  instrument string that will be used for calculating
 *                  instrument stats
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
 * Revisions:   PS 4/7/2018
 *                Added parameters QCFLAGNAME and QCFLAGVALS. The latter is
 *                supposed to contain a list of allowed flag values.
 *
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

/**
 */
void reader_z_xbt_profile(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    char* varname = NULL;
    char* lonname = NULL;
    char* latname = NULL;
    char* zname = NULL;
    char* stdname = NULL;
    char* estdname = NULL;
    char* timename = NULL;
    double varshift = 0.0;
    double mindepth = 0.0;
    char instrument[MAXSTRLEN];
    int nqcflags = 0;
    char** qcflagname = NULL;
    uint32_t* qcflagvals = NULL;

    int ncid;
    size_t nobs;
    int varid_var = -1, varid_lon = -1, varid_lat = -1, varid_z = -1, varid_std = -1, varid_estd = -1, varid_time = -1;
    double* lon = NULL;
    double lon_add_offset, lon_scale_factor;
    double lon_fill_value = NAN;
    double* lat = NULL;
    double lat_add_offset, lat_scale_factor;
    double lat_fill_value = NAN;
    double* z = NULL;
    double z_add_offset, z_scale_factor;
    double z_fill_value = NAN;
    double* var = NULL;
    double var_fill_value = NAN;
    double var_add_offset = NAN, var_scale_factor = NAN;
    double var_estd = NAN;
    double* std = NULL;
    double std_add_offset = NAN, std_scale_factor = NAN;
    double std_fill_value = NAN;
    double* estd = NULL;
    double estd_add_offset = NAN, estd_scale_factor = NAN;
    double estd_fill_value = NAN;
    uint32_t** qcflag = NULL;
    int have_time = 1;
    int singletime = -1;
    double* time = NULL;
    double time_add_offset = NAN, time_scale_factor = NAN;
    double time_fill_value = NAN;
    char tunits[MAXSTRLEN];
    double tunits_multiple = NAN, tunits_offset = NAN;
    int i, nobs_read;
    size_t singleton = 1;
    int ndim_in_var = 2;
    int ndim_in_gcoord = 1;
    int ndim_in_zcoord = 1;

    strcpy(instrument, meta->product);
    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "VARNAME") == 0)
            varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "TIMENAME") == 0)
            timename = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LONNAME") == 0)
            lonname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LATNAME") == 0)
            latname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ZNAME") == 0)
            zname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "STDNAME") == 0)
            stdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ESTDNAME") == 0)
            estdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "VARSHIFT") == 0) {
            if (!str2double(meta->pars[i].value, &varshift))
                enkf_quit("%s: can not convert VARSHIFT = \"%s\" to double\n", meta->prmfname, meta->pars[i].value);
            enkf_printf("        VARSHIFT = %f\n", varshift);
        } else if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
            if (!str2double(meta->pars[i].value, &mindepth))
                enkf_quit("%s: can not convert MINDEPTH = \"%s\" to double\n", meta->prmfname, meta->pars[i].value);
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
        else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }
    get_qcflags(meta, &nqcflags, &qcflagname, &qcflagvals);

    if (varname == NULL)
        enkf_quit("reader_z_xbt_profile(): %s VARNAME not specified", fname);
    else
        enkf_printf("        VARNAME = %s\n", varname);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid_var);
    ncw_inq_vardims(ncid, varid_var, ndim_in_var, NULL, &nobs);


    lonname = get_lonname(ncid, lonname);
    if (lonname != NULL) {
        enkf_printf("        LONNAME = %s\n", lonname);
        ncw_inq_varid(ncid, lonname, &varid_lon);
        ncw_check_vardims(ncid, varid_lon, ndim_in_gcoord, &singleton);
    } else
        enkf_quit("reader_z_xbt_profile(): %s: could not find longitude variable", fname);

    latname = get_latname(ncid, latname);
    if (latname != NULL) {
        enkf_printf("        LATNAME = %s\n", latname);
        ncw_inq_varid(ncid, latname, &varid_lat);
        ncw_check_vardims(ncid, varid_lat, ndim_in_gcoord, &singleton);
    } else
        enkf_quit("reader_z_xbt_profile(): %s: could not find latitude variable", fname);

    zname = get_zname(ncid, zname);
    if (zname != NULL) {
        enkf_printf("        ZNAME = %s\n", zname);
        ncw_inq_varid(ncid, zname, &varid_z);
        ncw_check_vardims(ncid, varid_z, ndim_in_zcoord, &nobs);
    } else
        enkf_quit("reader_xzy_scattered(): %s: could not find Z variable", fname);

    lon = malloc(sizeof(double));
    lat = malloc(sizeof(double));
    z = malloc(nobs * sizeof(double));
    ncw_get_var_double(ncid, varid_lon, lon);
    ncw_get_var_double(ncid, varid_lat, lat);
    ncw_get_var_double(ncid, varid_z, z);
    
    if (ncw_att_exists(ncid, varid_lon, "_FillValue")) 
        ncw_get_att_double(ncid, varid_lon, "_FillValue", &lon_fill_value);
    
    if (lon[0] != lon_fill_value) {
        if (ncw_att_exists(ncid, varid_lon, "add_offset")) {
            ncw_get_att_double(ncid, varid_lon, "add_offset", &lon_add_offset);
            ncw_get_att_double(ncid, varid_lon, "scale_factor", &lon_scale_factor);
            lon[0] = lon[0] * lon_scale_factor + lon_add_offset;
        }
    }

    if (ncw_att_exists(ncid, varid_lat, "_FillValue")) 
        ncw_get_att_double(ncid, varid_lat, "_FillValue", &lat_fill_value);
    
    if (lat[0] != lat_fill_value) {
        if (ncw_att_exists(ncid, varid_lat, "add_offset")) {
            ncw_get_att_double(ncid, varid_lat, "add_offset", &lat_add_offset);
            ncw_get_att_double(ncid, varid_lat, "scale_factor", &lat_scale_factor);
            lat[0] = lat[0] * lat_scale_factor + lat_add_offset;
        }
    }

    if (ncw_att_exists(ncid, varid_z, "_FillValue")) 
        ncw_get_att_double(ncid, varid_z, "_FillValue", &z_fill_value);
    
    if (ncw_att_exists(ncid, varid_z, "add_offset")) {
        ncw_get_att_double(ncid, varid_z, "add_offset", &z_add_offset);
        ncw_get_att_double(ncid, varid_z, "scale_factor", &z_scale_factor);
        for (i = 0; i < nobs; ++i) {
            if (z[i] != z_fill_value)
                z[i] = z[i] * z_scale_factor + z_add_offset;
        }
    }
    // PRocess a positive/negative ZNAME depth variable 
    int move_inside_water = 1;
    int status;
    int is_positive;
    char positive_info[MAXSTRLEN];
    char *strinside;

    status = nc_get_att_text(ncid, varid_z, "positive", positive_info);
    strinside = strstr(positive_info, "down");

    is_positive = ((status == 0) && (strinside != NULL));
    if (is_positive)
      move_inside_water = -1;
    else {
      float valid_min, valid_max;
      status = nc_get_att_float(ncid, varid_z, "valid_min", &valid_min);
      status += nc_get_att_float(ncid, varid_z, "valid_max", &valid_max);
      if (status == 0) {
        is_positive = ((valid_min < 0) && (valid_max > 0) &&
                       (abs(valid_min) < abs(valid_max)));
        if (is_positive)
          move_inside_water = -1;
      } else {
        // Assumes positive measurments are in the ocean - AODN data.
        enkf_printf("Warning: Assuming ZNAME variable is positive down.\n");
        move_inside_water = -1;
      }
    }

    var = malloc(nobs * sizeof(double));
    ncw_get_var_double(ncid, varid_var, var);
    if (ncw_att_exists(ncid, varid_var, "_FillValue")) 
        ncw_get_att_double(ncid, varid_var, "_FillValue", &var_fill_value);
    if (ncw_att_exists(ncid, varid_var, "add_offset")) {
        ncw_get_att_double(ncid, varid_var, "add_offset", &var_add_offset);
        ncw_get_att_double(ncid, varid_var, "scale_factor", &var_scale_factor);

        for (i = 0; i < nobs; ++i) {
            if (var[i] != var_fill_value)
                var[i] = var[i] * var_scale_factor + var_add_offset;
        }
    }
    // Put data in the water
    for (i = 0; i < nobs; ++i)
      z[i] *= move_inside_water;


    if (stdname != NULL)
        ncw_inq_varid(ncid, stdname, &varid_std);
    else if (ncw_var_exists(ncid, "std"))
        ncw_inq_varid(ncid, "std", &varid_std);
    if (varid_std >= 0) {
        ncw_check_vardims(ncid, varid_std, ndim_in_var, &nobs);
        std = malloc(nobs * sizeof(double));
        ncw_get_var_double(ncid, varid_std, std);
        if (ncw_att_exists(ncid, varid_std, "_FillValue"))
            ncw_get_att_double(ncid, varid_std, "_FillValue", &std_fill_value);
        if (ncw_att_exists(ncid, varid_std, "add_offset")) {
            ncw_get_att_double(ncid, varid_std, "add_offset", &std_add_offset);
            ncw_get_att_double(ncid, varid_std, "scale_factor", &std_scale_factor);

            for (i = 0; i < nobs; ++i)
                if (std[i] != std_fill_value)
                    std[i] = std[i] * std_scale_factor + std_add_offset;
        }
    }

    if (estdname != NULL)
        ncw_inq_varid(ncid, estdname, &varid_estd);
    else if (ncw_var_exists(ncid, "error_std"))
        ncw_inq_varid(ncid, "error_std", &varid_estd);
    if (varid_estd >= 0) {
        ncw_check_vardims(ncid, varid_estd, ndim_in_var, &nobs);
        estd = malloc(nobs * sizeof(double));
        ncw_get_var_double(ncid, varid_estd, estd);
        if (ncw_att_exists(ncid, varid_estd, "_FillValue"))
            ncw_get_att_double(ncid, varid_estd, "_FillValue", &estd_fill_value);
        if (ncw_att_exists(ncid, varid_estd, "add_offset")) {
            ncw_get_att_double(ncid, varid_estd, "add_offset", &estd_add_offset);
            ncw_get_att_double(ncid, varid_estd, "scale_factor", &estd_scale_factor);

            for (i = 0; i < nobs; ++i)
                if (estd[i] != estd_fill_value)
                    estd[i] = estd[i] * estd_scale_factor + estd_add_offset;
        }
    }

    if (std == NULL && estd == NULL)
        if (ncw_att_exists(ncid, varid_var, "error_std")) {
            ncw_check_attlen(ncid, varid_var, "error_std", 1);
            ncw_get_att_double(ncid, varid_var, "error_std", &var_estd);
        }

    if (nqcflags > 0) {
        int varid = -1;

        qcflag = alloc2d(nqcflags, nobs, sizeof(int32_t));
        for (i = 0; i < nqcflags; ++i) {
            ncw_inq_varid(ncid, qcflagname[i], &varid);
            ncw_check_vardims(ncid, varid, ndim_in_var, &nobs);
            ncw_get_var_uint(ncid, varid, qcflag[i]);
        }
    }

    timename = get_timename(ncid, timename);
    if (timename != NULL) {
        enkf_printf("        TIMENAME = %s\n", timename);
        ncw_inq_varid(ncid, timename, &varid_time);
    } else {
        enkf_printf("        reader_z_xbt_profile(): %s: no TIME variable\n", fname);
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
            time = malloc(sizeof(double));
        } else {
            singletime = 0;
            assert(timelen == nobs);
            time = malloc(nobs * sizeof(double));
        }

        ncw_get_var_double(ncid, varid_time, time);
        if (ncw_att_exists(ncid, varid_time, "_FillValue"))
            ncw_get_att_double(ncid, varid_time, "_FillValue", &time_fill_value);
        if (ncw_att_exists(ncid, varid_time, "add_offset")) {
            ncw_get_att_double(ncid, varid_time, "add_offset", &time_add_offset);
            ncw_get_att_double(ncid, varid_time, "scale_factor", &time_scale_factor);

            for (i = 0; i < nobs; ++i)
                if (time[i] != time_fill_value)
                    time[i] = time[i] * time_scale_factor + time_add_offset;
        }
        ncw_get_att_text(ncid, varid_time, "units", tunits);
        tunits_convert(tunits, &tunits_multiple, &tunits_offset);
    }

    ncw_close(ncid);

    nobs_read = 0;
    for (i = 0; i < nobs; ++i) {
        observation* o;
        obstype* ot;
        int ii;

        if (lon[0] == lon_fill_value || isnan(lon[0]) || lat[0] == lat_fill_value || isnan(lat[0]))
            continue;
        if (z[i] == z_fill_value || isnan(z[i]) || var[i] == var_fill_value || isnan(var[i]) || (std != NULL && (std[i] == std_fill_value || isnan(std[i]))) || (estd != NULL && (estd[i] == estd_fill_value || isnan(estd[i]))) || (have_time && !singletime && (time[i] == time_fill_value || isnan(time[i]))))
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
        o->value = var[i] + varshift;
        if (estd == NULL)
            o->std = var_estd;
        else {
            if (std == NULL)
                o->std = estd[i];
            else
                o->std = (std[i] > estd[i]) ? std[i] : estd[i];
        }
        o->lon = lon[0];
        o->lat = lat[0];
        o->depth = z[i];
        o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax))
            o->status = STATUS_OUTSIDEOBSDOMAIN;
        if (o->status == STATUS_OK)
            o->status = grid_z2fk(g, o->fi, o->fj, o->depth, &o->fk);
        else
            o->fk = NAN;
        o->model_depth = NAN;   /* set in obs_add() */
        if (have_time)
            o->date = ((singletime) ? time[0] : time[i]) * tunits_multiple + tunits_offset;
        else
            o->date = NAN;
        o->aux = -1;

        obs->nobs++;
    }
    enkf_printf("        nobs = %d\n", nobs_read);

    free(lon);
    free(lat);
    free(z);
    free(var);
    if (std != NULL)
        free(std);
    if (estd != NULL)
        free(estd);
    if (time != NULL)
        free(time);
    if (nqcflags > 0) {
        free(qcflagname);
        free(qcflagvals);
        free(qcflag);
    }
}
