/******************************************************************************
 *
 * File:        reader_xy_gridded_hfradar.c        
 *
 * Created:     10/07/2018
 *
 * Author:      Pavel Sakov & Hugo Oliveira
 *
 * Description: A Generic reader for IMOS gridded surface velocity observations from both WERA HF-RADAR and SEASONDE. 
 *                It currently assumes the following:
 *              - there is only one data record (2D field);
 *              - longitude is the inner ("fast") coordinate of the variable apart from the singleton record dimension (time).
 *                There are a number of parameters that must (++) or can be
 *              specified if they differ from the default value (+). Some
 *              parameters are optional (-):
 *              - TIMENAME ("time") (+)
 *              - VARNAME ("UVARNAME | VVARNAME") (++)
 *              - ROTATION (True 1 | False 0)  (+)
 *                  a boolean to rotate the vectors to the current grid when reading.
 *                  if True, UVARNAME,VVARNAME are required, as well as a variable in the grid.prm [ANGLE]
 *              - U_VARNAME ("UCUR") (-)
 *                  Name of the variable for the zonal velocity
 *              - V_VARNAME ("VCUR") (-)
 *                  Name of the variable for the meridional velocity
 *              - U_STD ("UCUR_sd") (-)
 *                  Standard deviation of zonal velocity
 *              - V_STD ("VCUR_sd") (-)
 *                  Standard deviation of meridonal velocity
 *              - GDOPNAME(-_)
 *                  name of the radar beam intersection angle
 *              TODO: remove this!?
 *              - NPOINTSNAME_1,NPOINTSNAME_2 ("npoints") (-)
 *                  number of collated observational points from each antenna at locations
 *              TODO: remove this!?
 *              - LONNAME ("lon" | "longitude") (+)
 *              - LATNAME ("lat" | "latitude") (+)
 *              - STDNAME ("std") (-)
 *                  internal variability of the collated data if U_STD and V_STD are not provided
 *              - ESTDNAME ("error_std") (-)
 *                  error STD; if absent then needs to be specified externally
 *                  in the observation data parameter file.
 *                  If both  STDNAME or ESTDNAME are absent, the code will try
 *                  to read an VARNAME.'error_std' attribute value.
 *              - VARSHIFT (-)
 *                  data offset to be added
 *              - U_VARSHIFT (-)
 *                  u data offset to be added to u values
 *              - V_VARSHIFT (-)
 *                  u data offset to be added to v values
 *              - MINDEPTH (-)
 *                  minimal allowed depth
 *              - INSTRUMENT (-)
 *                  instrument string that will be used for calculating
 *                  instrument stats
 *              - QCFLAGNAME ("UCUR_quality_control" | "VCUR_quality_control") (-)
 *                  flags of the data
 *              - QCFLAGVALS (0,1,2,3) (-)
 *                  valid flags to include as valid obs
 *              - U_QCFLAGNAME ("UCUR_quality_control") (-)
 *                  flags of the data
 *              - V_QCFLAGNAME ("UCUR_quality_control") (-)
 *                  flags of the data

 *
 * Revisions:   HO 15/7/2018
 *                Added angle rotation logics.
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

#define TYPE_DOUBLE 0
#define TYPE_SHORT 1


/**
 */
void reader_xy_gridded_hfradar(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    int isperiodic_i = grid_isperiodic_i(g);
    int** numlevels = grid_getnumlevels(g);
    float** grid_angle = grid_getangle(g);
    int grid_ni = 0, grid_nj = 0; 

    char* varname = NULL;
    char* u_varname = NULL;
    char* v_varname = NULL;
    char* u_stdvarname = NULL;
    char* v_stdvarname = NULL;
    char* lonname = NULL;
    char* latname = NULL;
    char* gdopname = NULL; 
    char* npointsname_1 = NULL;
    char* npointsname_2 = NULL;
    char* stdname = NULL;
    char* estdname = NULL;
    char* timename = NULL;
    char* qcflagname = NULL;
    char* u_qcflagname = NULL;
    char* v_qcflagname = NULL;
    char* rotation_flag = NULL;

    int ncid;
    int ndim_var, ndim_xy;
    size_t dimlen_var[3], dimlen_xy[2];

    uint32_t qcflagvals = 0;
    float varshift = 0.0;
    float u_varshift = 0.0;
    float v_varshift = 0.0;
    double mindepth = 0.0;
    char instrument[MAXSTRLEN];

    int iscurv = -1;
    int need_rotation = -1;
    size_t ni = 0, nj = 0, n = 0, n_var = 0;
    int varid_lon = -1, varid_lat = -1;
    double* lon = NULL;
    double* lat = NULL;
    int varid_gdop = -1, varid_npoints_1 = -1, varid_npoints_2 = -1;
    int varid_var = -1, varid_std = -1, varid_estd = -1, varid_qcflag = -1, varid_time = -1;
    int varid_u = -1, varid_v =-1, varid_u_std = -1, varid_v_std = -1, varid_u_qcflag = -1, varid_v_qcflag = -1;
    float* var = NULL;
    float* u_var = NULL;
    float* v_var = NULL;
    float* u_std = NULL;
    float* v_std = NULL;
    float var_fill_value = NAN;
    float var_add_offset = NAN, var_scale_factor = NAN;
    
    float u_var_fill_value = NAN;
    float u_var_add_offset = NAN, u_var_scale_factor = NAN;
    float v_var_fill_value = NAN;
    float v_var_add_offset = NAN, v_var_scale_factor = NAN;
    
    float u_std_fill_value = NAN;
    float u_std_add_offset = NAN, u_std_scale_factor = NAN;
    float v_std_fill_value = NAN;
    float v_std_add_offset = NAN, v_std_scale_factor = NAN;

    double var_estd = NAN;
    float* gdop = NULL;
    short* npoints_1 = NULL;
    short* npoints_2 = NULL;
    float* std = NULL;
    float std_add_offset = NAN, std_scale_factor = NAN;
    float std_fill_value = NAN;
    float* estd = NULL;
    float estd_add_offset = NAN, estd_scale_factor = NAN;
    float estd_fill_value = NAN;
    int32_t* qcflag = NULL;
    int32_t* u_qcflag = NULL;
    int32_t* v_qcflag = NULL;

    int have_time = 1;
    int singletime = -1;
    float* time = NULL;
    float time_add_offset = NAN, time_scale_factor = NAN;
    float time_fill_value = NAN;
    char tunits[MAXSTRLEN];
    double tunits_multiple = NAN, tunits_offset = NAN;
    int i, nobs_read;

    strcpy(instrument, meta->product);
    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "VARNAME") == 0)
            varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ROTATION") == 0)
            rotation_flag = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "U_VARNAME") == 0)
            u_varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "V_VARNAME") == 0)
            v_varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "U_STD") == 0)
            u_stdvarname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "V_STD") == 0)
            v_stdvarname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "TIMENAME") == 0)
            timename = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "GDOPNAME") == 0)
            gdopname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "NPOINTSNAME_1") == 0)
            npointsname_1 = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "NPOINTSNAME_2") == 0)
            npointsname_2 = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LONNAME") == 0)
            lonname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LATNAME") == 0)
            latname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "STDNAME") == 0)
            stdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ESTDNAME") == 0)
            estdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "QCFLAGNAME") == 0)
            qcflagname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "U_QCFLAGNAME") == 0)
            u_qcflagname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "V_QCFLAGNAME") == 0)
            v_qcflagname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "QCFLAGVALS") == 0) {
            char seps[] = " ,";
            char* line = meta->pars[i].value;
            char* token;
            int val;

            qcflagvals = 0;
            while ((token = strtok(line, seps)) != NULL) {
                if (!str2int(token, &val))
                    enkf_quit("%s: could not convert QCFLAGVALS entry \"%s\" to integer", meta->prmfname, token);
                if (val < 0 || val > 31)
                    enkf_quit("%s: QCFLAGVALS entry = %d (supposed to be in [0,31] interval", meta->prmfname, val);
                qcflagvals |= 1 << val;
                line = NULL;
            }
            if (qcflagvals == 0)
                enkf_quit("%s: no valid flag entries found after QCFLAGVALS\n", meta->prmfname);
        } else if (strcasecmp(meta->pars[i].name, "VARSHIFT") == 0) {
            if (!str2float(meta->pars[i].value, &varshift))
                enkf_quit("%s: can not convert VARSHIFT = \"%s\" to float\n", meta->prmfname, meta->pars[i].value);
            enkf_printf("        VARSHIFT = %s\n", meta->pars[i].value);
        } else if (strcasecmp(meta->pars[i].name, "U_VARSHIFT") == 0) {
            if (!str2float(meta->pars[i].value, &u_varshift))
                enkf_quit("%s: can not convert U_VARSHIFT = \"%s\" to float\n", meta->prmfname, meta->pars[i].value);
            enkf_printf("        U_VARSHIFT = %s\n", meta->pars[i].value);
        } else if (strcasecmp(meta->pars[i].name, "V_VARSHIFT") == 0) {
            if (!str2float(meta->pars[i].value, &v_varshift))
                enkf_quit("%s: can not convert V_VARSHIFT = \"%s\" to float\n", meta->prmfname, meta->pars[i].value);
            enkf_printf("        V_VARSHIFT = %s\n", meta->pars[i].value);

        } else if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
            if (!str2double(meta->pars[i].value, &mindepth))
                enkf_quit("%s: can not convert MINDEPTH = \"%s\" to double\n", meta->prmfname, meta->pars[i].value);
            enkf_printf("        MINDEPTH = %f\n", mindepth);
        } else if (strcasecmp(meta->pars[i].name, "INSTRUMENT") == 0) {
            strncpy(instrument, meta->pars[i].value, MAXSTRLEN);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    if (varname == NULL)
        enkf_quit("reader_xy_gridded(): %s variable name not specified", fname);

    if (rotation_flag != NULL) {
        int bool_rotation_flag, valid_rotation_flag;
        bool_rotation_flag = atoi(rotation_flag);
        valid_rotation_flag = (bool_rotation_flag == 0 || bool_rotation_flag == 1);
        
        if (valid_rotation_flag) {
            need_rotation = bool_rotation_flag == 1;
            grid_getdims(g, &grid_ni, &grid_nj, NULL);
            if (need_rotation == 1 && grid_angle == NULL)
                enkf_quit("reader_xy_gridded_hfradar(): %s cannot rotate vectors without ANGLE parameter defined within grid prm file",fname);    
            enkf_printf("\t#Rotation of vectors enabled\n");
        }
        else
            enkf_quit("reader_xy_gridded_hfradar(): %s invalid ROTATION parameter. Need to be a boolean digit");
    }

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid_var);
    ncw_inq_vardims(ncid, varid_var, 3, &ndim_var, dimlen_var);
    if (ndim_var == 3) {
        int dimid[3];
        size_t nr;

        ncw_inq_vardimid(ncid, varid_var, dimid);
        ncw_inq_dimlen(ncid, dimid[0], &nr);
        if (nr != 1)
            enkf_quit("reader_xy_gridded(): %d records (currently only one is allowed)", nr);
        n_var = dimlen_var[1] * dimlen_var[2];
    } else if (ndim_var == 2) {
        if (nc_hasunlimdim(ncid))
            enkf_quit("reader_xy_gridded(): %s: %s: not enough spatial dimensions (must be 2)", fname, varname);
        n_var = dimlen_var[0] * dimlen_var[1];
    } else if (ndim_var != 2)
        enkf_quit("reader_xy_gridded(): %s: # dimensions = %d (must be 2 or 3 with only one record)", fname, ndim_var);

    if (lonname != NULL)
        ncw_inq_varid(ncid, lonname, &varid_lon);
    else if (ncw_var_exists(ncid, "lon"))
        ncw_inq_varid(ncid, "lon", &varid_lon);
    else if (ncw_var_exists(ncid, "longitude"))
        ncw_inq_varid(ncid, "longitude", &varid_lon);
    else if (ncw_var_exists(ncid, "LONGITUDE"))
        ncw_inq_varid(ncid, "LONGITUDE", &varid_lon);
    else
        enkf_quit("reader_xy_gridded(): %s: could not find longitude variable", fname);
    ncw_inq_vardims(ncid, varid_lon, 2, &ndim_xy, dimlen_xy);
    if (ndim_xy == 1) {
        iscurv = 0;
        ni = dimlen_xy[0];
    } else if (ndim_xy == 2) {
        iscurv = 1;
        ni = dimlen_xy[1];
        nj = dimlen_xy[0];
    } else
        enkf_quit("reader_xy_gridded(): %s: coordinate variable \"%s\" has neither 1 or 2 dimensions", fname, lonname);

    if (latname != NULL)
        ncw_inq_varid(ncid, latname, &varid_lat);
    else if (ncw_var_exists(ncid, "lat"))
        ncw_inq_varid(ncid, "lat", &varid_lat);
    else if (ncw_var_exists(ncid, "latitude"))
        ncw_inq_varid(ncid, "latitude", &varid_lat);
    else if (ncw_var_exists(ncid, "LATITUDE"))
        ncw_inq_varid(ncid, "LATITUDE", &varid_lat);
    else
        enkf_quit("reader_xy_gridded(): %s: could not find latitude variable", fname);
    if (iscurv == 0) {
        ncw_check_varndims(ncid, varid_lat, 1);
        ncw_inq_vardims(ncid, varid_lat, 1, NULL, &nj);
    } else
        ncw_check_vardims(ncid, varid_lat, 2, dimlen_xy);

    enkf_printf("        (ni, nj) = (%u, %u)\n", ni, nj);
    n = ni * nj;
    if (n != n_var)
        enkf_quit("reader_xy_gridded(): %s: dimensions of variable \"%s\" do not match coordinate dimensions", fname, varname);
    if (dimlen_var[ndim_var - 1] != ni)
        enkf_quit("reader_xy_gridded(): %s: %s: longitude must be the inner coordinate", fname, varname);

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

    if (gdopname != NULL)
        ncw_inq_varid(ncid, gdopname, &varid_gdop);
    else if (ncw_var_exists(ncid, "GDOP"))
        ncw_inq_varid(ncid, "GDOP", &varid_gdop);

    if (varid_gdop >= 0) {
        gdop = malloc(n * sizeof(float));
        ncw_get_var_float(ncid, varid_gdop, gdop);
        float gdop_fill_value, gdop_add_offset, gdop_scale_factor;
        if (ncw_att_exists(ncid, varid_gdop, "_FillValue"))
            ncw_get_att_float(ncid, varid_gdop, "_FillValue", &gdop_fill_value);
        if (ncw_att_exists(ncid, varid_gdop, "add_offset")) {
            ncw_get_att_float(ncid, varid_gdop, "add_offset", &gdop_add_offset);
            ncw_get_att_float(ncid, varid_gdop, "scale_factor", &gdop_scale_factor);
            }
    }

    if (npointsname_1 != NULL)
        ncw_inq_varid(ncid, npointsname_1, &varid_npoints_1);
    else if (ncw_var_exists(ncid, "NOBS1"))
        ncw_inq_varid(ncid, "NOBS1", &varid_npoints_1);
    if (varid_npoints_1 >= 0) {
        npoints_1 = malloc(n * sizeof(short));
        ncw_get_var_short(ncid, varid_npoints_1, npoints_1);
    }

    if (npointsname_2 != NULL)
        ncw_inq_varid(ncid, npointsname_2, &varid_npoints_2);
    else if (ncw_var_exists(ncid, "NOBS2"))
        ncw_inq_varid(ncid, "NOBS2", &varid_npoints_2);
    if (varid_npoints_2 >= 0) {
        npoints_2 = malloc(n * sizeof(short));
        ncw_get_var_short(ncid, varid_npoints_2, npoints_2);
    }

    if (need_rotation == 1) {

        if (u_varname != NULL ) 
            ncw_inq_varid(ncid, u_varname, &varid_u);
        else if (ncw_var_exists(ncid, "UCUR")) {
            ncw_inq_varid(ncid, "UCUR", &varid_u);
            u_varname = "UCUR";
        }
        if (varid_u >= 0) {
            u_var = malloc(n * sizeof(float));
            ncw_get_var_float(ncid,varid_u,u_var);
            if (ncw_att_exists(ncid,varid_u, "_FillValue"))
                ncw_get_att_float(ncid, varid_u, "_FillValue", &u_var_fill_value);
            if (ncw_att_exists(ncid,varid_u,"add_offset")) {
                ncw_get_att_float(ncid, varid_u, "add_offset", &u_var_add_offset);
                ncw_get_att_float(ncid, varid_u, "scale_factor", &u_var_scale_factor);
                }
        }
        else
            enkf_quit("reader_xy_gridded(): %s: Variable %s not found.",fname, u_varname);
        
        if (v_varname != NULL) 
            ncw_inq_varid(ncid, v_varname, &varid_v);
        else if (ncw_var_exists(ncid, "VCUR")) {
            ncw_inq_varid(ncid, "VCUR", &varid_v);
            v_varname = "VCUR";
        }
        if (varid_v >= 0) {
            v_var = malloc(n * sizeof(float));
            ncw_get_var_float(ncid,varid_v,v_var);
            if (ncw_att_exists(ncid,varid_v, "_FillValue"))
                ncw_get_att_float(ncid,varid_v, "_FillValue", &v_var_fill_value);
            if (ncw_att_exists(ncid,varid_v,"add_offset")) {
                ncw_get_att_float(ncid, varid_v, "add_offset", &v_var_add_offset);
                ncw_get_att_float(ncid, varid_v, "scale_factor", &v_var_scale_factor);
                }
        }
        else
            enkf_quit("reader_xy_gridded(): %s: Variable %s not found.",fname, v_varname);
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

    if (u_stdvarname != NULL) 
        ncw_inq_varid(ncid, u_stdvarname, &varid_u_std);
    else if (ncw_var_exists(ncid, "UCUR_sd"))
        ncw_inq_varid(ncid, "UCUR_sd", &varid_u_std);
    if (varid_u_std >= 0) {
        u_std = malloc(n * sizeof(float));
        ncw_get_var_float(ncid, varid_u_std, u_std);
        if (ncw_att_exists(ncid, varid_u_std, "_FillValue"))
            ncw_get_att_float(ncid, varid_u_std, "_FillValue", &u_std_fill_value);
        if (ncw_att_exists(ncid, varid_u_std, "add_offset")) {
            ncw_get_att_float(ncid, varid_u_std, "add_offset", &u_std_add_offset);
            ncw_get_att_float(ncid, varid_u_std, "scale_factor", &u_std_scale_factor);
        }
    }

    if (v_stdvarname != NULL) 
        ncw_inq_varid(ncid, v_stdvarname, &varid_v_std);
    else if (ncw_var_exists(ncid, "vcur_sd"))
        ncw_inq_varid(ncid, "vcur_sd", &varid_v_std);
    if (varid_v_std >= 0) {
        v_std = malloc(n * sizeof(float));
        ncw_get_var_float(ncid, varid_v_std, v_std);
        if (ncw_att_exists(ncid, varid_v_std, "_FillValue"))
            ncw_get_att_float(ncid, varid_v_std, "_FillValue", &v_std_fill_value);
        if (ncw_att_exists(ncid, varid_v_std, "add_offset")) {
            ncw_get_att_float(ncid, varid_v_std, "add_offset", &v_std_add_offset);
            ncw_get_att_float(ncid, varid_v_std, "scale_factor", &v_std_scale_factor);
        }
    }

    if (qcflagname != NULL) {
        ncw_inq_varid(ncid, qcflagname, &varid_qcflag);
        qcflag = malloc(n * sizeof(int32_t));
        ncw_get_var_int(ncid, varid_qcflag, qcflag);
    }

    if (u_qcflagname != NULL) 
        ncw_inq_varid(ncid, u_qcflagname, &varid_u_qcflag);
    else if (ncw_var_exists(ncid, "UCUR_quality_control"))
        ncw_inq_varid(ncid, "UCUR_quality_control", &varid_u_qcflag);
    if (varid_u_qcflag >= 0) {
        u_qcflag = malloc(n * sizeof(int32_t));
        ncw_get_var_int(ncid, varid_u_qcflag, u_qcflag);
    }

    if (v_qcflagname != NULL) 
        ncw_inq_varid(ncid, v_qcflagname, &varid_v_qcflag);
    else if (ncw_var_exists(ncid, "VCUR_quality_control"))
        ncw_inq_varid(ncid, "VCUR_quality_control", &varid_v_qcflag);
    if (varid_v_qcflag >= 0) {
        v_qcflag = malloc(n * sizeof(int32_t));
        ncw_get_var_int(ncid, varid_v_qcflag, v_qcflag);
    }

    if (timename != NULL)
        ncw_inq_varid(ncid, timename, &varid_time);
    else if (ncw_var_exists(ncid, "time"))
        ncw_inq_varid(ncid, "time", &varid_time);
    else if (ncw_var_exists(ncid, "TIME"))
        ncw_inq_varid(ncid, "TIME", &varid_time);
    else {
        enkf_printf("        reader_xy_gridded(): %s: no TIME variable\n", fname);
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

    ncw_close(ncid);

    nobs_read = 0;
    for (i = 0; i < (int) n; ++i) {
        observation* o;
        obstype* ot;
          int invalid_gdop, invalid_npoints1, invalid_npoints2, invalid_var, invalid_std, invalid_estd, invalid_u_var, invalid_v_var, invalid_time;

          invalid_gdop = (gdop != NULL && (gdop[i] > 180.0 || gdop[i] < 0.0));
          invalid_npoints1 = (npoints_1 != NULL &&  (npoints_1[i] == 0 || npoints_1[i]< 0));
          invalid_npoints2 = (npoints_2 != NULL && (npoints_2[i] == 0 || npoints_2[i] < 0));
          invalid_var = var[i] == var_fill_value || isnan(var[i]);
          invalid_std = (std != NULL && (std[i] == std_fill_value || isnan(std[i])));
          invalid_estd = (estd != NULL && (estd[i] == estd_fill_value || isnan(estd[i])));
          invalid_u_var = (u_var != NULL && (u_var[i] == u_var_fill_value || isnan(u_var[i])));
          invalid_v_var = (v_var != NULL && (v_var[i] == v_var_fill_value || isnan(v_var[i])));
          invalid_time = (have_time && !singletime && (time[i] == time_fill_value || isnan(time[i])));
 
        if (invalid_gdop || invalid_npoints1 || invalid_npoints2 || invalid_var || invalid_std || invalid_estd || invalid_u_var || invalid_v_var || invalid_time)
             continue;

        if (qcflag != NULL)
            /* general qcflags has precedence in discarding obs*/
            if ( (qcflagvals != (qcflagvals | 1<<qcflag[i])) || (u_qcflag != NULL && (qcflagvals != (qcflagvals | 1<<u_qcflag[i]))) || (v_qcflag != NULL && (qcflagvals != (qcflagvals | 1<<v_qcflag[i]))) )
                continue;

        /* [u,v]_qcflags has precedence in values. */
        if (u_qcflagname != NULL && v_qcflagname != NULL) {
            if ( (u_qcflag != NULL && (qcflagvals != (qcflagvals | 1<<u_qcflag[i]))) || (v_qcflag != NULL && (qcflagvals != (qcflagvals | 1<<v_qcflag[i]))) )
                 continue;
        }
        else if (u_qcflagname != NULL) {
            if (u_qcflag != NULL && (qcflagvals != (qcflagvals | 1<<u_qcflag[i])))
                continue;
        }
        else if (v_qcflagname != NULL) {
            if (v_qcflag != NULL && (qcflagvals != (qcflagvals | 1<<v_qcflag[i])))
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
        
        /* Shifted the reading order for HF-Radar to obtain indexes */
        if (iscurv == 0) {
            o->lon = lon[i % ni];
            o->lat = lat[i / ni];
        } else {
            o->lon = lon[i];
            o->lat = lat[i];
        }
        o->depth = 0.0;
        o->fk = (double) ksurf;
        o->model_depth = NAN;   /* set in obs_add() */
        o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax))
            o->status = STATUS_OUTSIDEOBSDOMAIN;

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
                        
            oangle = interpolate2d(o->fi, o-> fj, grid_ni, grid_nj, grid_angle, numlevels, isperiodic_i ) * DEG2RAD;

            if (strcmp(varname,u_varname) == 0)
                o->value = zonal*cos(oangle) - meridional*sin(oangle);
            if (strcmp(varname,v_varname) == 0)
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
            else
                enkf_quit("Aborted. error_std is missing for for product %s in obs.prm.",meta->product);
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

        if (have_time) {
            float t = (singletime) ? time[0] : time[i];

            if (!isnan(time_add_offset))
                o->date = (double) (t * time_scale_factor + time_add_offset) * tunits_multiple + tunits_offset;
            else
                o->date = (double) t* tunits_multiple + tunits_offset;
        } else
            o->date = NAN;

        o->aux = -1;

        obs->nobs++;
    }
    enkf_printf("        nobs = %d\n", nobs_read);

    free(lon);
    free(lat);
    free(var);
    if (gdop != NULL)
        free(gdop);
    if (npoints_1 != NULL) 
        free(npoints_1);
    if (npoints_2 != NULL)
        free(npoints_2);
    if (u_var != NULL)
        free(u_var);
    if (v_var != NULL)
        free(v_var);
    if (std != NULL)
        free(std);
    if (estd != NULL)
        free(estd);
    if (u_std != NULL)
        free(u_std);
    if (v_std != NULL)
        free(v_std);
    if (qcflag != NULL)
        free(qcflag);
    if (u_qcflag != NULL)
        free(u_qcflag);
    if (v_qcflag != NULL)
        free(v_qcflag);
    if (time != NULL)
        free(time);
}
