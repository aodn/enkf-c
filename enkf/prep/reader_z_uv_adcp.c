/******************************************************************************
 *
 * File:        reader_z_uv_adcp.c
 *
 * Created:     08/03/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *              Hugo Oliveira
 *              AODN
 *
 * Description: A Generic reader for Velocities measured by ADCP from IMOS.
 *                It currently assumes the following:
 *              - There is a rotation mode set by the user.
 *              - The variable can be of arbitrary dimensions as long as only two dimensions are non-singleton - TIME and DEPTH above the instrument.
 *                There are a number of parameters that must (++) or can be
 *              specified if they differ from the default value (+). Some
 *              parameters are optional (-):
 *              - VARNAME (++)
 *              - ZNAME (-)
 *              - HNAME (-)
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
 *               HO 30/1/2019
 *                Adapted from hf-radar reader.
 *               HO 19/2/2019
 *                Fixed general problems in the DEVL service.
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
void reader_z_uv_adcp(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
//    int ksurf = grid_getsurflayerid(g);

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
    char *zname = NULL;
    char *hname = NULL;

    char* stdname = NULL;
    char* estdname = NULL;
    char* timename = NULL;
    int nqcflags = 0;
    char** qcflagname = NULL;

    /* UV rotation*/
    char* rotation_flag = NULL;
    int need_rotation = -1;
    uint32_t* qcflagvals = 0;

    int ncid;
    int nc;
    int nt,nz;
    float varshift = 0.0;

    /* UV varshifts */
    float u_varshift = 0.0;
    float v_varshift = 0.0;

    double mindepth = 0.0;
    char instrument[MAXSTRLEN];
    float instrument_depth = -9999;
    int missing_depth = 0;

    size_t dimlen_z[NC_MAX_VAR_DIMS], dimlen_h[NC_MAX_VAR_DIMS], dimlen_var[NC_MAX_VAR_DIMS];
    int varid_lon = -1, varid_lat = -1;
    double* lon = NULL;
    double lon_add_offset, lon_scale_factor;
    double lon_fill_value = NAN;
    double* lat = NULL;
    double lat_add_offset, lat_scale_factor;
    double lat_fill_value = NAN;
    int varid_var = -1, varid_std = -1, varid_estd = -1, varid_time = -1;
    int varid_z = -1, varid_h = -1;

    

    double *z = NULL;
    double *rz = NULL;
    double z_add_offset, z_scale_factor;
    double *rh = NULL;
    double h_add_offset, h_scale_factor;


    float* var = NULL;
    float var_fill_value = NAN;
    float var_add_offset = NAN, var_scale_factor = NAN;
    double var_estd = NAN;
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
    double* time = NULL;
    double* rtime = NULL;
    double time_add_offset = NAN, time_scale_factor = NAN;
    double time_fill_value = NAN;
    char tunits[MAXSTRLEN];
    double tunits_multiple = NAN, tunits_offset = NAN;
    int i, nobs, nobs_read;

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
        else if (strcasecmp(meta->pars[i].name, "LONNAME") == 0)
            lonname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LATNAME") == 0)
            latname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ZNAME") == 0)
            zname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "HNAME") == 0)
            hname = meta->pars[i].value;
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
        enkf_quit("reader_z_uv_adcp(): %s: VARNAME not specified", fname);
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
                enkf_quit("reader_z_uv_adcp(): %s cannot rotate vectors without ANGLE parameter defined within grid prm file",fname);
            else
                enkf_printf("\t#Rotation of vectors enabled\n");
        }
        else
            enkf_quit("reader_z_uv_adcp(): %s invalid ROTATION parameter. Need to be a boolean digit");
    }

    ncw_open(fname, NC_NOWRITE, &ncid);
    // Rules here are: Vel variables are two dimensional but can ignore singleton dimensions
    //
    //
    //
    int ndim_in_lon = 0;
    int ndim_in_lat = 0;
    int ndim_in_z = 0;
    int ndim_in_h = 0;
    int ndim_in_var = 0;

    int lon_dimids[NC_MAX_VAR_DIMS];
    int lat_dimids[NC_MAX_VAR_DIMS];
    int z_dimids[NC_MAX_VAR_DIMS];
    int h_dimids[NC_MAX_VAR_DIMS];
    int var_dimids[NC_MAX_VAR_DIMS];

    for (i = 0; i < NC_MAX_VAR_DIMS; ++i) {
        lon_dimids[i] = -1;
        lat_dimids[i] = -1;
        z_dimids[i] = -1;
        h_dimids[i] = -1;
        var_dimids[i] = -1;
    }

    lonname = get_lonname(ncid, lonname);
    if (lonname != NULL) {
        enkf_printf("        LONNAME = %s\n", lonname);
        ncw_inq_varid(ncid, lonname, &varid_lon);
        ncw_inq_vardimid(ncid, varid_lon, lon_dimids);
        for (i = 0; i < NC_MAX_VAR_DIMS; ++i) {
            if (lon_dimids[i] == -1)
                break;
            else
                ndim_in_lon += 1;
        }
    } else
        enkf_quit("reader_z_uv_adcp(): %s: could not find longitude variable\n",
                fname);

    latname = get_latname(ncid, latname);
    if (latname != NULL) {
        enkf_printf("        LATNAME = %s\n", latname);
        ncw_inq_varid(ncid, latname, &varid_lat);
        ncw_inq_vardimid(ncid, varid_lat, lat_dimids);
        for (i = 0; i < NC_MAX_VAR_DIMS; ++i) {
            if (lat_dimids[i] == -1)
                break;
            else
                ndim_in_lat += 1;
        }
    } else
        enkf_quit("reader_z_uv_adcp(): %s: could not find latitude variable\n",
                fname);

    if (zname != NULL) {
        enkf_printf("        ZNAME = %s\n", zname);
        int status = nc_inq_varid(ncid, zname, &varid_z);
        if (status != NC_NOERR) {
          enkf_printf("        Warning: ZNAME is missing.\n");
          enkf_printf("        reading ZNAME as \"instrument_nominal_depth\"...\n");
          ncw_get_att_float(ncid, NC_GLOBAL, "instrument_nominal_depth",
                  &instrument_depth);
          if (instrument_depth < 0.0)
              enkf_printf("        WARNING: Global Attribute "
                      "\"instrument_nominal_depth\" is negative\n");
          else
              instrument_depth = instrument_depth * -1;
          ndim_in_z = 1;
          missing_depth = 1;
        }
        else {
            zname = get_zname(ncid, zname);
            ncw_inq_varid(ncid, zname, &varid_z);
            ncw_inq_vardimid(ncid, varid_z, z_dimids);
        }
    } else 
        enkf_quit("Could not discover ZNAME.\n");
    for (i = 0; i < NC_MAX_VAR_DIMS; ++i) {
        if (z_dimids[i] == -1)
            break;
        else
            ndim_in_z += 1;
    }


    hname = get_hname(ncid, hname);
    if (hname != NULL) {
        enkf_printf("        HNAME = %s\n", hname);
        ncw_inq_varid(ncid, hname, &varid_h);
        ncw_inq_vardimid(ncid, varid_h, h_dimids);
    } 
    else 
        enkf_quit("Could not discover HNAME.\n");
    for (i = 0; i < NC_MAX_VAR_DIMS; ++i) {
        if (h_dimids[i] == -1)
            break;
        else
            ndim_in_h += 1;
    }

    ncw_inq_varid(ncid, varname, &varid_var);
    ncw_inq_vardimid(ncid, varid_var, var_dimids);
    for (i = 0; i < NC_MAX_VAR_DIMS; ++i) {
        if (var_dimids[i] == -1)
            break;
        else
            ndim_in_var += 1;
    }
    nobs = 1;
    ncw_inq_vardims(ncid, varid_var, ndim_in_var, NULL, dimlen_var);
    for (i = 0; i < ndim_in_var; ++i) {
        if (dimlen_var[i] < 0)
            break;
        else
            nobs *= dimlen_var[i];
    }
    if (nobs == 1)
        enkf_quit("reader_z_uv_adcp(): %s: No observations\n");


    // End of dimension gathering

    int gcoord_mismatch = ndim_in_lon != ndim_in_lat;
    if (gcoord_mismatch)
        enkf_quit("reader_z_uv_adcp(): %s: dimension number mismatch between lon/lat",
                fname);

    int non_dim_gcoord = (ndim_in_lon == 0 && ndim_in_lat == 0);
    int singleton_dim_gcoord = (ndim_in_lon == 1 && ndim_in_lat == 1);
    int gcoord_is_singleton = (non_dim_gcoord || singleton_dim_gcoord);
    if (gcoord_is_singleton) {
        lon = malloc(sizeof(double));
        lat = malloc(sizeof(double));
    } else
        enkf_quit("reader_z_uv_adcp(): %s: Longitude has %s and latitude "
                "has %s number of dimensions. This reader accepts only "
                "dimensionless or singleton coordinate variables.\n",
                fname, ndim_in_lon, ndim_in_lat);

    ncw_get_var_double(ncid, varid_lon, lon);
    if (ncw_att_exists(ncid, varid_lon, "_FillValue"))
        ncw_get_att_double(ncid, varid_lon, "_FillValue", &lon_fill_value);

    if (lon[0] != lon_fill_value) {
        if (ncw_att_exists(ncid, varid_lon, "add_offset")) {
            ncw_get_att_double(ncid, varid_lon, "add_offset", &lon_add_offset);
            ncw_get_att_double(ncid, varid_lon, "scale_factor", &lon_scale_factor);
            lon[0] = lon[0] * lon_scale_factor + lon_add_offset;
        }
    }

    ncw_get_var_double(ncid, varid_lat, lat);
    if (ncw_att_exists(ncid, varid_lat, "_FillValue"))
        ncw_get_att_double(ncid, varid_lat, "_FillValue", &lat_fill_value);

    if (lat[0] != lat_fill_value) {
        if (ncw_att_exists(ncid, varid_lat, "add_offset")) {
            ncw_get_att_double(ncid, varid_lat, "add_offset", &lat_add_offset);
            ncw_get_att_double(ncid, varid_lat, "scale_factor", &lat_scale_factor);
            lat[0] = lat[0] * lat_scale_factor + lat_add_offset;
        }
    }

    int invalid_zh_dims = ((ndim_in_z != 1) || (ndim_in_h != 1));
    if (invalid_zh_dims) {
        if (ndim_in_z > 1) {
            int counter = 0;
            for (i =0; i < ndim_in_z ; ++i) {
                if (dimlen_z[i] > 1)
                    ++counter;
                if (counter > 1)

                    enkf_quit("reader_z_uv_adcp(): %s: Variable %s got [%d] dimensions with more than 1 non-singleton dimension",fname, zname, ndim_in_z);
            }
        }
        if (ndim_in_h > 1) {
            int counter = 0;
            for (i =0; i < ndim_in_h ; ++i) {
                if (dimlen_h[i] > 1)
                    ++counter;
                if (counter > 1)
                    enkf_quit("reader_z_uv_adcp(): %s: Variable %s got [%d] dimensions with more than 1 non-singleton dimension",fname, hname, ndim_in_h);
            }
        }
    }


    if (!missing_depth)
        ncw_inq_vardims(ncid, varid_z, NC_MAX_VAR_DIMS, NULL, dimlen_z);
    else {
        enkf_printf("        WARNING: assuming first variable dimension is along Z.\n");
        dimlen_z[0] = dimlen_var[0] ;
    }



    ncw_inq_vardims(ncid, varid_h, NC_MAX_VAR_DIMS, NULL, dimlen_h);
    ncw_inq_vardims(ncid, varid_var, NC_MAX_VAR_DIMS, NULL, dimlen_var);

    int construct_z = dimlen_z[0] != dimlen_h[0];
    if (construct_z) {
        int invalid_h_coord = dimlen_h[0] > dimlen_z[0];
        if (invalid_h_coord)
            enkf_quit("reader_z_uv_adcp(): %s: Invalid Z/H sizes. %s[%d] > %s[%d]",fname, hname, ndim_in_h, zname, ndim_in_z);

        int move_inside_water = 1;
        rz = malloc(dimlen_z[0]*sizeof(double));
        if (!missing_depth) {
            int status;
            int is_positive;
            char positive_info[MAXSTRLEN];
            char *strinside;

            status = nc_get_att_text(ncid, varid_z, "positive", positive_info);
            strinside = strstr(positive_info, "down");
            is_positive = ((status == 0) && (strinside != NULL));
            if (is_positive)
                ;
            else {
                float valid_min,valid_max;
                status = nc_get_att_float(ncid, varid_z, "valid_min", &valid_min);
                status += nc_get_att_float(ncid, varid_z, "valid_max", &valid_max);
                if (status == 0) {
                    is_positive = ((valid_min < 0) && (valid_max > 0) &&
                            (abs(valid_min) < abs(valid_max)));
                    if (is_positive)
                        ;
                } else {
                    // Assumes positive measurments are in the ocean - AODN data.
                    enkf_printf("         Warning: Assuming ZNAME variable is positive down.\n");
                    move_inside_water = -1;
                }
            }
            ncw_get_var_double(ncid, varid_z, rz);
        } else
            for (i = 0; i < dimlen_z[0]; ++i) 
                rz[i] = instrument_depth;



        rh = malloc(dimlen_h[0]*sizeof(double));
        ncw_get_var_double(ncid, varid_h, rh); 

        int z_transform = 0;
        if (ncw_att_exists(ncid, varid_z, "add_offset")) {
            ncw_get_att_double(ncid, varid_z, "add_offset", &z_add_offset);
            ncw_get_att_double(ncid, varid_z, "scale_factor", &z_scale_factor);
            z_transform = 1;
        }

        int h_transform = 0;
        if (ncw_att_exists(ncid, varid_h, "add_offset")) {
            ncw_get_att_double(ncid, varid_h, "add_offset", &h_add_offset);
            ncw_get_att_double(ncid, varid_h, "scale_factor", &h_scale_factor);
            h_transform = 1;
        }

        double zt,ht;
        int j;
        nc = 0;
        z = malloc(nobs * sizeof(double));
        for (i = 0; i < dimlen_z[0]; ++i) {
            for (j = 0; j < dimlen_h[0]; ++j) {
                if (z_transform)
                    zt = rz[i]*z_scale_factor + z_add_offset;
                else
                    zt = rz[i];

                if (h_transform)
                    ht = rh[j]*h_scale_factor + h_add_offset;
                else
                    ht = rh[j];

                z[nc] = zt*move_inside_water - ht;
                ++nc;
            }
        }
        int valid_obs = (nc == nobs);
        if (valid_obs) {
            nt = dimlen_z[0];
            nz = dimlen_h[0];
        }
        else
            enkf_quit("reader_z_uv_adcp(): %s: NOBS mismatch.",fname);
    } else 
        enkf_quit("reader_z_uv_adcp(): %s: Could not defined a Z/H coordinate. %s[%d], %s[%d], %s[%d] are not valid triplet.",fname, varname, ndim_in_var, zname, ndim_in_z, hname, ndim_in_h);


    // postpone scale/offset conversions to the loop
    var = malloc(nobs * sizeof(float));
    ncw_get_var_float(ncid, varid_var, var);
    if (ncw_att_exists(ncid, varid_var, "add_offset")) {
        ncw_get_att_float(ncid, varid_var, "add_offset", &var_add_offset);
        ncw_get_att_float(ncid, varid_var, "scale_factor", &var_scale_factor);
    }
    if (ncw_att_exists(ncid, varid_var, "_FillValue"))
        ncw_get_att_float(ncid, varid_var, "_FillValue", &var_fill_value);

    if (stdname != NULL)
        ncw_inq_varid(ncid, stdname, &varid_std);
    else if (ncw_var_exists(ncid, "std"))
        ncw_inq_varid(ncid, "std", &varid_std);
    if (varid_std >= 0) {
        std = malloc(nobs * sizeof(float));
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
        estd = malloc(nobs * sizeof(float));
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
        qcflag = alloc2d(nqcflags, nobs, sizeof(int));
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
        enkf_quit("reader_z_uv_adcp(): %s: no TIME variable\n", fname);
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
            assert(timelen == nt);
            time = malloc(nobs * sizeof(double));
            rtime = malloc(nt * sizeof(double));
        }

        ncw_get_var_double(ncid, varid_time, rtime);
        if (ncw_att_exists(ncid, varid_time, "_FillValue"))
            ncw_get_att_double(ncid, varid_time, "_FillValue", &time_fill_value);
        if (ncw_att_exists(ncid, varid_time, "add_offset")) {
            ncw_get_att_double(ncid, varid_time, "add_offset", &time_add_offset);
            ncw_get_att_double(ncid, varid_time, "scale_factor", &time_scale_factor);
        }
        ncw_get_att_text(ncid, varid_time, "units", tunits);
        tunits_convert(tunits, &tunits_multiple, &tunits_offset);
        int loc=-1;
        int j;
        for (i = 0 ; i < nt; ++i) {
            for (j = 0; j < nz; ++j) {
                loc += 1;
                time[loc] = rtime[i];
            }
        }
    }

    /* Perform loads required for rotation */
    if (need_rotation == 1) {

        if (u_varname != NULL)
            ncw_inq_varid(ncid, u_varname, &varid_u);
        else if (ncw_var_exists(ncid, "UCUR")) {
            ncw_inq_varid(ncid, "UCUR", &varid_u);
            u_varname = "UCUR";
        } else
            enkf_quit("reader_z_uv_adcp(): %s: Variable %s not found.",fname, u_varname);

        u_var = malloc(nobs * sizeof(float));
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
            enkf_quit("reader_z_uv_adcp(): %s: Variable %s not found.",fname, v_varname);

        v_var = malloc(nobs * sizeof(float));
        ncw_get_var_float(ncid,varid_v,v_var);
        if (ncw_att_exists(ncid,varid_v, "_FillValue"))
            ncw_get_att_float(ncid,varid_v, "_FillValue", &v_var_fill_value);
        if (ncw_att_exists(ncid,varid_v,"add_offset")) {
            ncw_get_att_float(ncid, varid_v, "add_offset", &v_var_add_offset);
            ncw_get_att_float(ncid, varid_v, "scale_factor", &v_var_scale_factor);
            }
    }

    enkf_printf("        (nt, nz) = (%u, %u)\n", nt, nz);
    ncw_close(ncid);

    nobs_read = 0;
    for (i = 0; i < (int) nobs; ++i) {
        observation* o;
        obstype* ot;
        int ii;
        int invalid_var, invalid_std, invalid_estd, invalid_time;
        int invalid_u_var, invalid_v_var;
        
        invalid_var = var[i] == var_fill_value || isnan(var[i]);
        invalid_std = (std != NULL && (std[i] == std_fill_value || isnan(std[i])));
        invalid_estd = (estd != NULL && (estd[i] == estd_fill_value || isnan(estd[i])));
        invalid_time = (have_time && !singletime && (time[i] == time_fill_value || isnan(time[i])));
        invalid_u_var = (u_var != NULL && (u_var[i] == u_var_fill_value || isnan(u_var[i])));
        invalid_v_var = (v_var != NULL && (v_var[i] == v_var_fill_value || isnan(v_var[i])));

        if (invalid_var || invalid_std || invalid_estd || invalid_time || invalid_u_var || invalid_v_var)
            continue;

        if (qcflagvals != NULL) {
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
        }

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
            o->fk =NAN;
        o->model_depth = NAN;   /* set in obs_add() */
        if (have_time) {
            double t = (singletime) ? time[0] : time[i];
            if (!isnan(time_add_offset))
                o->date = (double) ((t * time_scale_factor + time_add_offset) * tunits_multiple + tunits_offset);
            else
                o->date = (double) t*tunits_multiple + tunits_offset;
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
    if (time != NULL)
        free(time);
    if (nqcflags > 0) {
        free(qcflagname);
        free(qcflagvals);
        free(qcflag);
    }
    if (construct_z) {
        free(rz);
        free(rh);
        free(z);
    }
    /* Free rotation related vars */
    if (u_var != NULL)
        free(u_var);
    if (v_var != NULL)
        free(v_var);
}
