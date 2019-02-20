/******************************************************************************
 *
 * File:        reader_z_adcp_temp.c
 *
 * Created:     06/02/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *              Hugo Oliveira
 *              AODN
 *
 * Description: Generic reader for IMOS ANMN ctd profile data.
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
void reader_z_adcp_temp(char *fname, int fid, obsmeta *meta, grid *g,
                              observations *obs) {
  char *varname = NULL;
  char *lonname = NULL;
  char *latname = NULL;
  char *zname = NULL;
  char *stdname = NULL;
  char *estdname = NULL;
  char *timename = NULL;
  double varshift = 0.0;
  double mindepth = 0.0;
  char instrument[MAXSTRLEN];
  int nqcflags = 0;
  char **qcflagname = NULL;
  uint32_t *qcflagvals = NULL;

  int ncid;
  size_t nobs;
  int varid_var = -1, varid_lon = -1, varid_lat = -1, varid_z = -1,
      varid_std = -1, varid_estd = -1, varid_time = -1;
  double *lon = NULL;
  double lon_add_offset, lon_scale_factor;
  double lon_fill_value = NAN;
  double *lat = NULL;
  double lat_add_offset, lat_scale_factor;
  double lat_fill_value = NAN;
  double *z = NULL;
  double z_add_offset, z_scale_factor;
  double z_fill_value = NAN;
  double *var = NULL;
  double var_fill_value = NAN;
  double var_add_offset = NAN, var_scale_factor = NAN;
  double var_estd = NAN;
  double *std = NULL;
  double std_add_offset = NAN, std_scale_factor = NAN;
  double std_fill_value = NAN;
  double *estd = NULL;
  double estd_add_offset = NAN, estd_scale_factor = NAN;
  double estd_fill_value = NAN;
  uint32_t **qcflag = NULL;
  int have_time = 1;
  int singletime = -1;
  double *time = NULL;
  double time_add_offset = NAN, time_scale_factor = NAN;
  double time_fill_value = NAN;
  char tunits[MAXSTRLEN];
  double tunits_multiple = NAN, tunits_offset = NAN;
  int i, nobs_read;

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
        enkf_quit("%s: can not convert VARSHIFT = \"%s\" to double\n",
                  meta->prmfname, meta->pars[i].value);
      enkf_printf("        VARSHIFT = %f\n", varshift);
    } else if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
      if (!str2double(meta->pars[i].value, &mindepth))
        enkf_quit("%s: can not convert MINDEPTH = \"%s\" to double\n",
                  meta->prmfname, meta->pars[i].value);
    } else if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
      if (!str2double(meta->pars[i].value, &mindepth))
        enkf_quit("%s: can not convert MINDEPTH = \"%s\" to double\n",
                  meta->prmfname, meta->pars[i].value);
      enkf_printf("        MINDEPTH = %f\n", mindepth);
    } else if (strcasecmp(meta->pars[i].name, "INSTRUMENT") == 0)
      strncpy(instrument, meta->pars[i].value, MAXSTRLEN);
    else if (strcasecmp(meta->pars[i].name, "QCFLAGNAME") == 0 ||
             strcasecmp(meta->pars[i].name, "QCFLAGVALS") == 0)
      /*
       * QCFLAGNAME and QCFLAGVALS are dealt with separately
       */
      ;
    else
      enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
  }
  get_qcflags(meta, &nqcflags, &qcflagname, &qcflagvals);

  if (varname == NULL)
    enkf_quit("reader_z_adcp_temp(): %s VARNAME not specified\n", fname);
  else
    enkf_printf("        VARNAME = %s\n", varname);

  ncw_open(fname, NC_NOWRITE, &ncid);
  ncw_inq_varid(ncid, varname, &varid_var);
  // The rules for this reader is the following:
  // lon/lat are non-dimension variables or singleton ones
  // variables have one or more dimensions but only one dimension can be
  // non-singleton

  int ndim_in_lon = 0;
  int ndim_in_lat = 0;
  int ndim_in_z = 0;
  int ndim_in_var = 0;

  int lon_dimids[NC_MAX_VAR_DIMS];
  int lat_dimids[NC_MAX_VAR_DIMS];
  int z_dimids[NC_MAX_VAR_DIMS];
  int var_dimids[NC_MAX_VAR_DIMS];

  for (i = 0; i < NC_MAX_VAR_DIMS; ++i) {
    lon_dimids[i] = -1;
    lat_dimids[i] = -1;
    z_dimids[i] = -1;
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
    enkf_quit(
        "reader_z_adcp_temp(): %s: could not find longitude variable\n",
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
    enkf_quit(
        "reader_z_adcp_temp(): %s: could not find latitude variable\n",
        fname);


  if (zname != NULL) {
    enkf_printf("        ZNAME = %s\n", zname);
    int status = nc_inq_varid(ncid,zname,&varid_z);
    if (status != NC_NOERR) {
      enkf_printf("        Warning: ZNAME is missing.\n");
      ndim_in_z = 0;
    } else {
      zname = get_zname(ncid, zname);
      ncw_inq_varid(ncid, zname, &varid_z);
      ncw_inq_vardimid(ncid, varid_z, z_dimids);
      for (i = 0; i < NC_MAX_VAR_DIMS; ++i) {
        if (z_dimids[i] == -1)
          break;
        else
          ndim_in_z += 1;
      }
      }
  } else
    enkf_quit("reader_z_adcp_temp(): %s: could not find Z variable\n",
              fname);

  if (ndim_in_lon != ndim_in_lat)
    enkf_quit("reader_z_adcp_temp(): %s: dimension number mismatch\n"
              "between lon/lat",
              fname);

  int non_dim_gcoord = (ndim_in_lon == 0 && ndim_in_lat == 0);
  int singleton_dim_gcoord = (ndim_in_lon == 1 && ndim_in_lat == 1);
  if (non_dim_gcoord || singleton_dim_gcoord) {
    lon = malloc(sizeof(double));
    lat = malloc(sizeof(double));
  } else
    enkf_quit("reader_z_adcp_temp(): %s: Longitude has %s and latitude "
              "has %s number of dimensions. This reader accepts only "
              "dimensionless or singleton coordinate variables.\n",
              fname, ndim_in_lon, ndim_in_lat);

  ncw_inq_vardimid(ncid, varid_var, var_dimids);
  for (i = 0; i < NC_MAX_VAR_DIMS; ++i) {
    if (var_dimids[i] == -1)
      break;
    else
      ndim_in_var += 1;
  }

  if (ndim_in_var != ndim_in_z) {
    int counter = 0;
    size_t test;
    for (i = 0; i < ndim_in_var; ++i) {
      ncw_inq_dimlen(ncid,var_dimids[i],&test);
      if (test != 1)
        ++counter;
      if (counter > 1)
        enkf_quit("reader_z_adcp_temp(): %s: Non-singleton dimension %d "
                  "in variable %s is not part of the xyz variables.\n",
                  fname, var_dimids[i], varname);
    }
  }

  ncw_inq_vardims(ncid, varid_var, ndim_in_var, NULL, &nobs);
  z = malloc(nobs * sizeof(double));


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

 
  int missing_depth;
  missing_depth = (varid_z == -1);
  if (missing_depth) {
    enkf_printf("Reading ZNAME as \"nominal_instrument_depth\"...\n");
    float instrument_depth;
    ncw_get_att_float(ncid,NC_GLOBAL,"instrument_nominal_depth", &instrument_depth);
    if (instrument_depth < 0.0)
      enkf_quit("Global Attribute \"instrument_nominal_depth\" is negative. Aborted\n");    
    for (i =0; i < nobs; i++) {
      z[i] = instrument_depth;
    }
  }
  else {
    ncw_get_var_double(ncid, varid_z, z);
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
    enkf_printf("        reader_z_adcp_temp(): %s: no TIME variable\n",
                fname);
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
    observation *o;
    obstype *ot;
    int ii;

    if (lon[0] == lon_fill_value || isnan(lon[0]) || lat[0] == lat_fill_value ||
        isnan(lat[0]))
      continue;
    if (z[i] == z_fill_value || isnan(z[i]) || var[i] == var_fill_value ||
        isnan(var[i]) ||
        (std != NULL && (std[i] == std_fill_value || isnan(std[i]))) ||
        (estd != NULL && (estd[i] == estd_fill_value || isnan(estd[i]))) ||
        (have_time && !singletime &&
         (time[i] == time_fill_value || isnan(time[i]))))
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
    if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax ||
                                     o->lat <= ot->ymin || o->lat >= ot->ymax))
      o->status = STATUS_OUTSIDEOBSDOMAIN;
    if (o->status == STATUS_OK)
      o->status = grid_z2fk(g, o->fi, o->fj, o->depth, &o->fk);
    else
      o->fk = NAN;
    o->model_depth = NAN; /* set in obs_add() */
    if (have_time)
      o->date =
          ((singletime) ? time[0] : time[i]) * tunits_multiple + tunits_offset;
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
