#!/usr/bin/env Rscript
# --~- BLISS_PREC1d.R -~--
# Bayesian statisticaL Interpolation for Spatial analySis of precipitation
# See the software repository here: 
#..............................................................................
#Copyright and license
# Copyright (C) 2018 MET Norway. The software is licensed under GPL version 3 
# or (at your option) any later version.
# https://www.gnu.org/licenses/gpl-3.0.en.html
# 
# History:
# 05.10.2018 - Cristian Lussana. Original code.
# -----------------------------------------------------------------------------
#
rm(list=ls())
#
# -----------------------------------------------------------------------------
# Libraries
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("sp"))
suppressPackageStartupMessages(library("raster"))
suppressPackageStartupMessages(library("igraph"))
suppressPackageStartupMessages(library("rgdal"))
suppressPackageStartupMessages(library("ncdf4"))
suppressPackageStartupMessages(library("dotnc"))
#options(warn = 2, scipen = 999)
options(scipen = 999)
# 
# -----------------------------------------------------------------------------
# Constants
# CRS strings (always useful for copy-and-paste)
proj4.llwgs84<-"+proj=longlat +datum=WGS84"
#proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
#proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
#proj4.aeqd<-"+proj=aeqd +lat_0=59.458925 +lon_0=10.564472 +R=6371000 +datum=WGS84"
#
#..............................................................................
# Functions
# + manage fatal error
boom<-function(str=NULL) {
  print("Fatal Error:")
  if (!is.null(str)) print(str)
  quit(status=1)
}

#+ check if two rasters match
rasters_match<-function(r1,r2) {
  return ( projection(r1)==projection(r2) &
           extent(r1)==extent(r2) &
           res(r1)[1]==res(r2)[1] & 
           res(r1)[2]==res(r2)[2])
}

#+
set_NAs_to_NULL<-function(x) {
  if (!is.null(x)) {
    if (is.na(x)) x<-NULL
  }
  x
}

#
debug_plots<-function() {
  png(file=paste0("rb_",formatC(l,width=2,flag="0"),".png"),
      width=800,height=800)
  image(rb,breaks=c(0,seq(0.001,0.8,length=25)),
        col=c("gray",rev(rainbow(24))))
  points(VecX,VecY,pch=19,col="black",cex=0.5)
  dev.off()
  png(file=paste0("innoh_",formatC(l,width=2,flag="0"),".png"),
      width=800,height=800)
  d<-yo_relan[ixwet]-yb[ixwet]
  hist(d,breaks=seq(-1000.025,1000,by=.05),xlim=c(-1,1))
  dev.off()
  png(file=paste0("plot_",formatC(l,width=2,flag="0"),".png"),
      width=800,height=800)
  plot(yo_relan,yb,pch=19,col="blue",xlim=c(0,1),ylim=c(0,1))
  ix9<-which(prId==9)
  points(yo_relan[ix9],yb[ix9],pch=19,col="red")
  lines(-1000:1000,-1000:1000,col="red",lwd=2)
  dev.off()
}

#
debug_01_plots<-function() {
  png(file=paste0("ra.png"),width=800,height=800)
  image(ra,breaks=c(0,seq(0.1,50,length=45)),col=c("gray",rev(rainbow(44))))
  dev.off()
  png(file=paste0("ra1.png"),width=800,height=800)
  image(ra,breaks=c(0,0.09999,1000),col=c("gray","blue"))
  points(VecX[ixwet],VecY[ixwet],pch=19,col="cyan",cex=0.5)
  points(VecX[ixdry],VecY[ixdry],pch=19,col="red",cex=0.5)
  dev.off()
}

#+ find the n-th largest element from each matrix row 
findRow <- function (x,n) {   
# order by row number then by value
  y<-t(x)
  array(y[order(col(y), y)], dim(y))[nrow(y) - (n-1), ]
}

#+ mean radial distance between an observation and its k-th closest obs
dobs_fun<-function(obs,k) {
# NOTE: k=1 will return 0 everywhere (1st closest obs is the obs itself)
  nobs<-length(obs$x)
  if (nobs<k) return(NA)
  # distance matrices
  disth<-(outer(obs$x,obs$x,FUN="-")**2.+
          outer(obs$y,obs$y,FUN="-")**2.)**0.5
  dobsRow<-findRow(x=disth,n=(nobs-k+1))
  mean(dobsRow)
}

# + replace elements of a string with date-time elements
`replaceDate`<-function(string=NULL,
                        date.str=NULL,
                        format="%Y-%m-%d %H:%M:%S") {
#------------------------------------------------------------------------------
  if (is.null(string) | is.null(date.str)) return(NULL)
  Rdate<-as.POSIXlt(str2Rdate(date.str,format=format))
  yyyy<-Rdate$year+1900
  mm<-formatC(Rdate$mon+1,width=2,flag="0")
  dd<-formatC(Rdate$mday,width=2,flag="0")
  hh<-formatC(Rdate$hour,width=2,flag="0")
  out<-gsub("yyyy",yyyy,string)
  out<-gsub("mm",formatC(mm,width=2,flag="0"),out)
  out<-gsub("dd",formatC(dd,width=2,flag="0"),out)
  out<-gsub("hh",formatC(hh,width=2,flag="0"),out)
  out
}

#+
`OI_RR_fast`<-function(yo,
                       yb,
                       xb,
                       xgrid,
                       ygrid,
                       zgrid=NULL,
                       VecX,
                       VecY,
                       VecZ=NULL,
                       Dh,
                       Dz=NULL) {
#------------------------------------------------------------------------------
  if (is.null(Dz)) {
    Dz<-100000
    VecZ<-rep(0,length(VecX))
    zgrid<-rep(0,length(xgrid))
  }
  no<-length(yo)
  ng<-length(xb)
  xa<-vector(mode="numeric",length=ng)
  vec<-vector(mode="numeric",length=no)
  d<-yo-yb
  out<-.C("oi_rr_first",no=as.integer(no), 
                        innov=as.double(d),
                        SRinv=as.numeric(InvD),
                        vec=as.double(vec) ) 
  vec[1:no]<-out$vec[1:no]
  rm(out)
  out<-.C("oi_rr_fast",ng=as.integer(ng),
                       no=as.integer(no),
                       xg=as.double(xgrid),
                       yg=as.double(ygrid),
                       zg=as.double(zgrid),
                       xo=as.double(VecX),
                       yo=as.double(VecY),
                       zo=as.double(VecZ),
                       Dh=as.double(Dh),
                       Dz=as.double(Dz),
                       xb=as.double(xb),
                       vec=as.double(vec),
                       xa=as.double(xa) )
  out$xa[1:ng]
}

#+
`OI_RR_var`<-function(yo,
                      yb,
                      xb,
                      gx,
                      gy,
                      ox,
                      oy,
                      Dh,
                      eps2
                      ) {
#------------------------------------------------------------------------------
  no<-length(yo)
  ng<-length(gx)
  xa<-vector(mode="numeric",length=ng)
  ya<-vector(mode="numeric",length=no)
  xa_errvar<-vector(mode="numeric",length=ng)
  ya_errvar<-vector(mode="numeric",length=no)
  o_errvar<-0
  xa[]<-0
  ya[]<-0
  xa_errvar[]<-0
  ya_errvar[]<-0
  out<-.C("oi_rr_var",ng=as.integer(ng),
                      no=as.integer(no),
                      SRinv=as.numeric(InvD),
                      eps2=as.double(eps2),
                      Dh=as.double(Dh),
                      gx=as.double(gx),
                      gy=as.double(gy),
                      ox=as.double(ox),
                      oy=as.double(oy),
                      yo=as.double(yo),
                      yb=as.double(yb),
                      xb=as.double(xb),
                      xa=as.double(xa),
                      ya=as.double(ya),
                      xa_errvar=as.double(xa_errvar),
                      ya_errvar=as.double(ya_errvar),
                      o_errvar=as.double(o_errvar))
  return( list(xa= out$xa[1:ng],
               ya= out$ya[1:no],
               xa_errvar= out$xa_errvar[1:ng],
               ya_errvar= out$ya_errvar[1:no],
               o_errvar= out$o_errvar))
}

#
boxcox<-function(x,lambda) {
  if (lambda==0) {
    return(log(x))
  } else {
    return((x**lambda-1)/lambda)
  }
}

#+
tboxcox4pdf_apply<-function(par,
                            lambda,
                            brrinf,
                            int=400) {
  mean<-par[1]
  sd<-par[2]
  if (mean<brrinf) return(0)
  if (sd==0) return(tboxcox(mean,lambda))
  aux<-qnorm(mean=mean,sd=sd,p=seq(1/int,(1-1/int),length=(int-1)))
  if (lambda!=0) aux[aux<(-1/lambda)]<-(-1/lambda)
  mean(tboxcox(aux,lambda=lambda))
}

#+
tboxcox<-function(x,lambda) {
  if (lambda==0) {
    return(exp(x))
  } else {
    return((1+lambda*x)**(1./lambda))
  }
}


#
#==============================================================================
# MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN -
#==============================================================================
t0<-Sys.time()
# [] Read command line arguments and/or set parameters to default
# create parser object
p <- arg_parser("snrr")
#------------------------------------------------------------------------------
# miscellaneous 
p <- add_argument(p, "--verbose",help="debug mode",flag=T,short="-v")
p <- add_argument(p, "--rrinf",
                  help="precipitation yes/no threshold",
                  type="numeric",
                  default=0.1)
#------------------------------------------------------------------------------
# cross-validation mode
p <- add_argument(p, "--cv_mode",
                  help="standard cross-validation mode",
                  flag=T)
p <- add_argument(p, "--loocv_mode",
                  help="leave-one-out cross-validation mode",
                  type="logical",
                  default=F)
p <- add_argument(p, "--idiv_instead_of_elev",
                  help="leave-one-out cross-validation mode",
                  type="logical",
                  default=F)
#------------------------------------------------------------------------------
# statistical interpolation mode
p <- add_argument(p, "--mode",
                  help="statistical interpolation scheme (\"OI_multiscale\",\"OI_firstguess\")",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# time-related variables
p <- add_argument(p, "--date_out",
                  help="date to write in output files (nc, %Y%m%d%H%M)",
                  type="character",
                  default=NULL)
p <- add_argument(p, "--date_out_fmt",
                  help="date out format",
                  type="character",
                  default="%Y%m%d%H%M")
p <- add_argument(p, "--time_bnds_string",
                  help="time bounds with respect to date_out (e.g., \"-1 day\" \"-1 min\")",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# dataset labels for IO
p <- add_argument(p, "--dset",
                  help="dataset title (e.g., NGCD)",
                  type="character",
                  default="nolab")
p <- add_argument(p, "--var",
                  help="variable abbreviation (e.g., RR)",
                  type="character",
                  default="noab")
#------------------------------------------------------------------------------
# OI_multiscale / OI_firstguess parameters
p <- add_argument(p, "--eps2",
                  help="ratio of observation to background error covariance",
                  type="numeric",
                  default=1)
p <- add_argument(p, "--Dh",
                  help="horizontal de-corellation length scale (km)",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--ovarc",
                  help="observation error variance correction factor",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--ovarc.prId",
                  help="observation error variance correction factor (prId)",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--prId.exclude",
                  help="observation provider identifiers to exclude",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--prId.cv",
                  help="observation provider identifiers to reserve for cross-validation",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
#
#
#------------------------------------------------------------------------------
# paths
p <- add_argument(p, "--path2src",
                  help="path to the shared objects (.so files)",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# input files
p <- add_argument(p, "--config.file",
                  help="configuration file",
                  type="character",
                  default=NULL,
                  short="cf")
p <- add_argument(p, "--iff_obs",
                  help="full file name for the observations (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg",
                  help="full file name for the first-guess field (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf",
                  help="full file name for the rescaling factor (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask",
                  help="full file name for the mask (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_black",
                  help="blacklist file",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# output files
p <- add_argument(p, "--off_stn",
                  help="full file name for output at station locations (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_ver",
                  help="full file name for output at station locations analysis vs obs - verif format (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_ver_fg",
                  help="full file name for output at station locations fist-guess vs obs - verif format (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd",
                  help="full file name for output at gridpoints (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_cvstn",
                  help="full file name for output at gridpoints (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_cvver",
                  help="full file name for output at gridpoints (nc)",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# customization of the observation file
p <- add_argument(p, "--iff_obs.sep",
                  help="separator character",
                  type="character",
                  default=";")
p <- add_argument(p, "--iff_obs.x",
                  help="easting coordinate name",
                  type="character",
                  default="x")
p <- add_argument(p, "--iff_obs.y",
                  help="northing coordinate name",
                  type="character",
                  default="y")
p <- add_argument(p, "--iff_obs.z",
                  help="elevation name",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.value",
                  help="variable name",
                  type="character",
                  default="value")
p <- add_argument(p, "--iff_obs.proj4",
                  help="proj4 string for the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.sourceId",
                  help="station identifier",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.prId",
                  help="provider identifier",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.dqc",
                  help="data quality control",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# Master grid definition
p <- add_argument(p, "--grid_master.x1",
                  help="easting coordinate of the first gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.xn",
                  help="easting coordinate of the last gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.y1",
                  help="northing coordinate of the first gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.yn",
                  help="northing coordinate of the last gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.resx",
                  help="grid spacing along the easting coordinate",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.resy",
                  help="grid spacing along the northing coordinate",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.proj4",
                  help="proj4 string for the master grid",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# Mask file netcdf parameters
p <- add_argument(p, "--iff_mask.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_mask.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_mask.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_mask.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_mask.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_mask.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
#------------------------------------------------------------------------------
# Rescaling factor file netcdf parameters
p <- add_argument(p, "--iff_rf.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.varname_lat",
                  help="name of the latitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.varname_lon",
                  help="name of the longitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_rf.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_rf.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
#------------------------------------------------------------------------------
# first-guess file netcdf parameters
p <- add_argument(p, "--iff_fg.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_fg.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_fg.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_fg.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_fg.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_fg.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
#------------------------------------------------------------------------------
# output file netcdf parameters
p <- add_argument(p, "--off_grd.grid",
                  help="grid type",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.varname",
                  help="variable name",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.varlongname",
                  help="variable long name",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.standardname",
                  help="variable standard name",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.varversion",
                  help="variable version",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.varunit",
                  help="variable unit",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.timesunit",
                  help="time unit",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.reference",
                  help="references",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.write_lonlat",
                  help="add latitude and longitude variables",
                  type="logical",
                  default=F)
p <- add_argument(p, "--off_grd.diground",
                  help="rounding digits",
                  type="numeric",
                  default=3)
p <- add_argument(p, "--off_grd.summary",
                  help="summary",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.sourcestring",
                  help="source string",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.title",
                  help="title",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.comment",
                  help="title",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.cell_methods",
                  help="title",
                  type="cell aggregation method",
                  default="none")
#------------------------------------------------------------------------------
# gaussian anamorphosis
p <- add_argument(p, "--transf",
                  help="transformation used in the gaussian anamorphosis (\"none\",\"Box-Cox\")",
                  type="character",
                  default="none")
p <- add_argument(p, "--transf.boxcox_lambda",
                  help="Box-Cox power parameter",
                  type="numeric",
                  default=NA)

#------------------------------------------------------------------------------
#
argv <- parse_args(p)
#
#------------------------------------------------------------------------------
# read configuration file
if (!is.na(argv$config.file)) {
  if (file.exists(argv$config.file)) {
    source(argv$config.file)
    argv_tmp<-append(argv,conf)
    names_argv_tmp<-names(argv_tmp)
    argv_def<-list()
    names_argv_def<-integer(0)
    k<-0
    for (i in 1:length(argv_tmp)) {
      if (names_argv_tmp[i] %in% names_argv_def) next
      k<-k+1
      j<-which(names_argv_tmp==names_argv_tmp[i])
      argv_def[[k]]<-argv_tmp[[j[length(j)]]]
      names_argv_def<-c(names_argv_def,names_argv_tmp[i])
    }
    names(argv_def)<-names_argv_def
    rm(argv_tmp,names_argv_tmp,names_argv_def)
    rm(argv)
    argv<-argv_def
    rm(argv_def)
  } else {
    print("WARNING: config file not found")
    print(argv$config.file)
  }
}
#
#------------------------------------------------------------------------------
# check input arguments
#
if (!(file.exists(argv$iff_obs))) boom(paste0("file not found ",argv$iff_obs))
if (argv$mode=="OI_multiscale") {
  if (argv$verbose) {
    if (!(file.exists(argv$iff_rf))) print("warning: file not found",argv$iff_rf)
  }
} else if (argv$mode=="OI_firstguess") {
  if (!(file.exists(argv$iff_fg))) boom(paste0("file not found ",argv$iff_fg))
} else {
  boom("error statistical interpolation scheme undefined")
}
# define/check paths and load external functions
if ( !(file.exists(argv$path2src)) ) 
  ext<-boom("path not found")
#
# load external C functions
dyn.load(file.path(argv$path2src,"oi_rr_first.so"))
dyn.load(file.path(argv$path2src,"oi_rr_fast.so"))
dyn.load(file.path(argv$path2src,"oi_rr_var.so"))
#
#------------------------------------------------------------------------------
# Create master grid
xmn<-as.numeric(argv$grid_master.x1)-as.numeric(argv$grid_master.resx)/2
xmx<-as.numeric(argv$grid_master.xn)+as.numeric(argv$grid_master.resx)/2
ymn<-as.numeric(argv$grid_master.y1)-as.numeric(argv$grid_master.resy)/2
ymx<-as.numeric(argv$grid_master.yn)+as.numeric(argv$grid_master.resy)/2
rmaster<-raster(extent(xmn,xmx,ymn,ymx),
                res=c(as.numeric(argv$grid_master.resx),
                      as.numeric(argv$grid_master.resy)),
                crs=argv$grid_master.proj4)
rmaster[]<-1
# use mask if provided
if (file.exists(argv$iff_mask)) {
  argv$iff_mask.epos<-set_NAs_to_NULL(argv$iff_mask.epos)
  argv$iff_mask.tpos<-set_NAs_to_NULL(argv$iff_mask.tpos)
  argv$iff_mask.e<-set_NAs_to_NULL(argv$iff_mask.e)
  if (argv$iff_mask.t=="none") argv$iff_mask.t<-nc4.getTime(argv$iff_mask)[1]
  raux<-try(read_dotnc(nc.file=argv$iff_mask,
                       nc.varname=argv$iff_mask.varname,
                       topdown=argv$iff_mask.topdown,
                       out.dim=list(ndim=argv$iff_mask.ndim,
                                    tpos=argv$iff_mask.tpos,
                                    epos=argv$iff_mask.epos,
                                    names=argv$iff_mask.names),
                       proj4=argv$iff_mask.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_mask.t,
                                      format=argv$iff_mask.tfmt,
                                      e=argv$iff_mask.e)))
  if (is.null(raux)) 
    boom("error reading the mask file")
  rmask<-raux$stack; rm(raux)
  if (!rasters_match(rmask,rmaster)) rmask<-projectRaster(rmask,rmaster)
  rmaster<-mask(rmaster,rmask); rm(rmask)
}
#
nx<-ncol(rmaster)
ny<-nrow(rmaster)
xmn<-xmin(rmaster)
xmx<-xmax(rmaster)
ymn<-ymin(rmaster)
ymx<-ymax(rmaster)
#
# extract all the cell values: zvalues[1] contains the rmaster[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
zvalues<-getValues(rmaster)
storage.mode(zvalues)<-"numeric"
xy<-xyFromCell(rmaster,1:ncell(rmaster))
x<-sort(unique(xy[,1]))
y<-sort(unique(xy[,2]),decreasing=T)
xgrid<-xy[,1]
ygrid<-xy[,2]
mask<-which(!is.na(zvalues))
ngrid<-length(mask)
# clean memory
rm(zvalues,xy)
# debug info
if (argv$verbose) {
  print("+---------------------------------------------------------------+")
  print("+ grid parameters")
  print(paste("nx ny dx dy",
    as.integer(nx),as.integer(ny),round(xres(rmaster),2),round(yres(rmaster),2)))
  print(paste("xmn xmx ymn ymx",
    round(xmn,2),round(xmx,2),round(ymn,2),round(ymx,2)))
  print(paste("# grid points (master/unmasked)=",as.integer(ngrid)))
}
#
#------------------------------------------------------------------------------
# read rescaling factor
if (file.exists(argv$iff_rf)) {
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
    print(paste("+ rescaling factor",argv$iff_rf))
  }
  argv$iff_rf.epos<-set_NAs_to_NULL(argv$iff_rf.epos)
  argv$iff_rf.tpos<-set_NAs_to_NULL(argv$iff_rf.tpos)
  argv$iff_rf.e<-set_NAs_to_NULL(argv$iff_rf.e)
  if (argv$iff_rf.t=="none") argv$iff_rf.t<-nc4.getTime(argv$iff_rf)[1]
  raux<-try(read_dotnc(nc.file=argv$iff_rf,
                       nc.varname=argv$iff_rf.varname,
                       topdown=argv$iff_rf.topdown,
                       out.dim=list(ndim=argv$iff_rf.ndim,
                                    tpos=argv$iff_rf.tpos,
                                    epos=argv$iff_rf.epos,
                                    names=argv$iff_rf.names),
                       proj4=argv$iff_rf.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_rf.t,
                                      format=argv$iff_rf.tfmt,
                                      e=argv$iff_rf.e)))
  if (is.null(raux)) 
    boom("error reading the rescaling file")
  rrf<-raux$stack; rm(raux)
  if (!rasters_match(rrf,rmaster)) {
    if (argv$iff_rf.varname_lat=="none") {
      rrf<-projectRaster(rrf,rmaster)
    } else {
      raux<-try(read_dotnc(nc.file=argv$iff_rf,
                       nc.varname=argv$iff_rf.varname_lat,
                       topdown=argv$iff_rf.topdown,
                       out.dim=list(ndim=argv$iff_rf.ndim,
                                    tpos=argv$iff_rf.tpos,
                                    epos=argv$iff_rf.epos,
                                    names=argv$iff_rf.names),
                       proj4=argv$iff_rf.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_rf.t,
                                      format=argv$iff_rf.tfmt,
                                      e=argv$iff_rf.e)))
      if (is.null(raux)) 
        boom("error reading the rescaling file (lat)")
      lat<-raux$data; rm(raux)
      raux<-try(read_dotnc(nc.file=argv$iff_rf,
                           nc.varname=argv$iff_rf.varname_lon,
                           topdown=argv$iff_rf.topdown,
                           out.dim=list(ndim=argv$iff_rf.ndim,
                                        tpos=argv$iff_rf.tpos,
                                        epos=argv$iff_rf.epos,
                                        names=argv$iff_rf.names),
                           proj4=argv$iff_rf.proj4,
                           nc.proj4=list(var=NULL,
                                         att=NULL),
                           selection=list(t=argv$iff_rf.t,
                                          format=argv$iff_rf.tfmt,
                                          e=argv$iff_rf.e)))
      if (is.null(raux)) 
        boom("error reading the rescaling file (lon)")
      lon<-raux$data; rm(raux)
      coord.new<-spTransform(SpatialPoints(cbind(lon,lat),
                                           proj4string=argv$iff_rf.proj4),
                             CRS(argv$grid_master.proj4))
      #
      rf<-getValues(rrf)
      rrfagg<-rasterize(coord.new,
                        aggregate(rmaster,fact=4),
                        rf)
      rrf<-mask( crop( disaggregate(rrfagg,fact=4,method="bilinear"),
                       rmaster),
                 rmaster)
      rf<-getValues(rrf)
    }
  }
}
#
#------------------------------------------------------------------------------
# read first guess
if (file.exists(argv$iff_fg)) {
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
    print(paste("+ first guess ",argv$iff_fg))
  }
  argv$iff_fg.epos<-set_NAs_to_NULL(argv$iff_fg.epos)
  argv$iff_fg.tpos<-set_NAs_to_NULL(argv$iff_fg.tpos)
  argv$iff_fg.e<-set_NAs_to_NULL(argv$iff_fg.e)
  if (argv$iff_fg.t=="none") argv$iff_fg.t<-nc4.getTime(argv$iff_fg)[1]
  raux<-try(read_dotnc(nc.file=argv$iff_fg,
                       nc.varname=argv$iff_fg.varname,
                       topdown=argv$iff_fg.topdown,
                       out.dim=list(ndim=argv$iff_fg.ndim,
                                    tpos=argv$iff_fg.tpos,
                                    epos=argv$iff_fg.epos,
                                    names=argv$iff_fg.names),
                       proj4=argv$iff_fg.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_fg.t,
                                      format=argv$iff_fg.tfmt,
                                      e=argv$iff_fg.e)))
  if (is.null(raux)) 
    boom("error reading the rescaling file")
  rfg<-raux$stack; rm(raux)
  if (!rasters_match(rfg,rmaster)) rfg<-projectRaster(rfg,rmaster)
  rfg<-mask(rfg,rmaster)
  xb0<-getValues(rfg)
  aix<-which(!is.na(xb0))
  xb<-xb0[aix]
  rm(xb0)
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
  }
}
#
#------------------------------------------------------------------------------
# Read observations
dat<-read.table(file=argv$iff_obs,
                header=T,
                sep=argv$iff_obs.sep,
                stringsAsFactors=F,
                strip.white=T)
# read x_orig,y_orig,value & set x,y 
varidxtmp<-match(c(argv$iff_obs.x,
                   argv$iff_obs.y,
                   argv$iff_obs.value),
                   names(dat))
if (any(is.na(varidxtmp))) {
  print("ERROR in the specification of the variable names")
  print(paste("        x=",argv$iff_obs.x))
  print(paste("        y=",argv$iff_obs.y))
  print(paste("    value=",argv$iff_obs.value))
  print("header of input file:")
  print(argv$argv$iff_obs)
  print(names(dat))
  quit(status=1)
}
data<-data.frame(dat[,varidxtmp])
names(data)<-c("x_orig","y_orig","value")
data$x_orig<-suppressWarnings(as.numeric(data$x_orig))
data$y_orig<-suppressWarnings(as.numeric(data$y_orig))
data$value<-suppressWarnings(as.numeric(data$value))
if (argv$iff_obs.proj4!=argv$grid_master.proj4) {
  xymaster<-spTransform(SpatialPoints(cbind(data$x_orig,data$y_orig),
                                      proj4string=CRS(argv$iff_obs.proj4)) ,
                        CRS(argv$grid_master.proj4))
  data$x<-attr(xymaster,"coords")[,1]
  data$y<-attr(xymaster,"coords")[,2]
  rm(xymaster)
} else {
  data$x<-data$x_orig
  data$y<-data$y_orig
}
# lat-lon coords are required by verif
if (argv$iff_obs.proj4!=proj4.llwgs84) {
  xyll<-spTransform(SpatialPoints(cbind(data$x_orig,data$y_orig),
                                  proj4string=CRS(argv$iff_obs.proj4)) ,
                    CRS(proj4.llwgs84))
  data$lon<-attr(xyll,"coords")[,1]
  data$lat<-attr(xyll,"coords")[,2]
  rm(xyll)
} else {
  data$lon<-data$x_orig
  data$lat<-data$y_orig
}
# read z 
if (argv$iff_obs.z!="none") {
  varidxtmp<-which(argv$iff_obs.z==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     z=",argv$iff_obs.z))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$z<-dat[,varidxtmp]
} else {
  data$z<-rep(0,length(data$x))
}
# read sourceId 
if (argv$iff_obs.sourceId!="none") {
  varidxtmp<-which(argv$iff_obs.sourceId==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     sourceId=",argv$iff_obs.sourceId))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$sourceId<-dat[,varidxtmp]
} else {
  data$sourceId<-1:length(data$x)
}
# read prId 
if (argv$iff_obs.prId!="none") {
  varidxtmp<-which(argv$iff_obs.prId==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     prId=",argv$iff_obs.prId))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$prId<-dat[,varidxtmp]
} else {
  data$prId<-rep(0,length(data$x))
}
# read dqc flag 
if (argv$iff_obs.dqc!="none") {
  varidxtmp<-which(argv$iff_obs.dqc==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     dqc=",argv$iff_obs.dqc))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$dqc<-dat[,varidxtmp]
} else {
  data$dqc<-rep(0,length(data$x))
}
# data, dataframe: x,y,z,x_orig,y_orig,value,sourceId,prId,dqc
#
if (file.exists(argv$iff_black)) {
  bstid<-read.csv(argv$iff_black,header=T,stringsAsFactors=F,strip.white=T)
} else {
  bstid<-integer(0)
}
# select only observation within the master grid
flag_in_master<-!is.na(extract(rmaster,cbind(data$x,data$y)))
flag_in_fg<-rep(T,length(data$x))
if (file.exists(argv$iff_fg)) 
  flag_in_fg<-!is.na(extract(rfg,cbind(data$x,data$y))) 
# on-the-fly dqc, used for testing
#  flag_in_fg<-!is.na(extract(rfg,cbind(data$x,data$y))) &
#              data$value > (extract(rfg,cbind(data$x,data$y))-0.2*extract(rfg,cbind(data$x,data$y)))
#CVmode
if (argv$cv_mode) {
# prId=1 MET-WMO stations
  ixcv<-which( data$dqc==0 & 
               data$prId %in% argv$prId.cv & 
               !is.na(data$value) &
               flag_in_master &
               flag_in_fg   )
  if (length(ixcv)==0) boom("ERROR \"cv_mode\" running without CV-observations ")
  VecX_cv<-data$x[ixcv]
  VecY_cv<-data$y[ixcv]
  VecXorig_cv<-data$x_orig[ixcv]
  VecYorig_cv<-data$y_orig[ixcv]
  VecLat_cv<-data$lat[ixcv]
  VecLon_cv<-data$lon[ixcv]
  VecZ_cv<-data$z[ixcv]
  VecS_cv<-data$sourceId[ixcv]
  yo_cv<-data$value[ixcv]
  if (exists("rrf")) yrf_cv<-extract(rrf,cbind(VecX_cv,VecY_cv),na.rm=T)
  data$value[ixcv]<-NA
  data$dqc[ixcv]<-999
  ncv<-length(ixcv)
}
if (any(!is.na(argv$prId.exclude))) {
  ix0<-which(data$dqc==0 & 
             !(data$prId %in% argv$prId.exclude) &
             flag_in_master &
             flag_in_fg   )
} else {
  ix0<-which(data$dqc==0 &
             flag_in_master &
             flag_in_fg   )
}
n0<-length(ix0)
# definitive station list
if (n0==0) boom("No observations. Stop here.")
VecX<-data$x[ix0]
VecY<-data$y[ix0]
VecXorig<-data$x_orig[ix0]
VecYorig<-data$y_orig[ix0]
VecLat<-data$lat[ix0]
VecLon<-data$lon[ix0]
VecZ<-data$z[ix0]
VecS<-data$sourceId[ix0]
yo<-data$value[ix0]
if (exists("rrf"))  yrf<-extract(rrf,cbind(VecX,VecY),na.rm=T)
prId<-data$prId[ix0]
ydqc.flag<-rep(0,length=n0)
rm(data)
ixwet<-which(yo>=argv$rrinf)
ixdry<-which(yo< argv$rrinf)
nwet<-length(ixwet)
ndry<-length(ixdry)
# Easy-peasy
#if (nwet==0) {
#  if (argv$verbose) print("no rain over the whole domain")
#  writeIfNoPrec()
#  if (argv$verbose) print("Success exit")
#  quit(status=0)
#}
# observation error variance correction factor
ovarc<-rep(1,n0)
if (any(!is.na(argv$ovarc.prId))) {
  for (i in 1:length(argv$ovarc.prId)) {
    if (any(prId==argv$ovarc.prId[i])) 
      ovarc[which(prId==argv$ovarc.prId[i])]<-argv$ovarc[i]
  }
}
if (argv$verbose) { 
  print("+---------------------------------------------------------------+")
  print(paste("#observations (wet/dry) =",n0,"(",nwet,"/",ndry,")"))
}
#
#------------------------------------------------------------------------------
# Set the OI multi-scale parameters
if (argv$mode=="OI_multiscale") {
  if (n0<100) {
    kseq<-c(2,3,4,seq(5,n0,by=5),n0)
  } else {
    kseq<-c(2,3,4,seq(5,100,by=5),seq(100,n0,by=200),n0)
  }
  kseq<-rev(unique(kseq))
  vecd_tmp<-vector()
  vecf<-vector()
  for (i in 1:length(kseq)) 
    vecd_tmp[i]<-round(dobs_fun(obs=data.frame(x=VecX,y=VecY),k=kseq[i])/1000,0)
  if (vecd_tmp[length(kseq)]>=15) {
    vecd_tmp[c(i+1,i+2,i+3)]<-c(12,8,5)
  } else if (vecd_tmp[length(kseq)]>=10) {
    vecd_tmp[c(i+1,i+3)]<-c(10,5)
  }
  vecf_tmp<-pmin(min(c(nx,ny))/5,pmax(1,round(vecd_tmp/5,0)))
  vecd<-vecd_tmp[which(!duplicated(vecf_tmp,fromLast=T))]
  kseq<-kseq[which(!duplicated(vecf_tmp,fromLast=T))]
  vecf<-round(vecd/5,0)
  vece<-rep(argv$eps2,length(vecf))
  nl<-length(vecd)
  if (length(which(vecf<=2))>0) {
    vece[which(vecf<=2)]<-rep(0.1,length(which(vecf<=2)))
  } else {
    vece[length(vece)]<-0.1
  }
  rm(vecf_tmp,vecd_tmp)
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
    print("kseq vecd vecf vece")
    print(" #    km    -    -")
    print(cbind(kseq,vecd,vecf,vece))
    print("+---------------------------------------------------------------+")
  }
}
#
#------------------------------------------------------------------------------
# compute Disth (symmetric) matrix: 
#  Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
Disth<-matrix(ncol=n0,nrow=n0,data=0.)
Disth<-(outer(VecY,VecY,FUN="-")**2.+
        outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
#
#------------------------------------------------------------------------------
# ANALYSIS
# ===> OI with background  <===
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
    print("Analysis")
  }
if (argv$mode=="OI_firstguess") {
  yb<-extract(rfg,
              cbind(VecX,VecY),
              method="bilinear")
  D<-exp(-0.5*(Disth/argv$Dh)**2.)
  diag(D)<-diag(D)+argv$eps2
  InvD<-chol2inv(chol(D))
  if (argv$transf=="none") {
    xa<-OI_RR_fast(yo=yo,
                   yb=yb,
                   xb=xb,
                   xgrid=xgrid[aix],
                   ygrid=ygrid[aix],
                   VecX=VecX,
                   VecY=VecY,
                   Dh=argv$Dh)
  } else if (argv$transf=="Box-Cox") {
    t00<-Sys.time()
    res<-OI_RR_var(yo=boxcox(yo,argv$transf.boxcox_lambda),
                   yb=boxcox(yb,argv$transf.boxcox_lambda),
                   xb=boxcox(xb,argv$transf.boxcox_lambda),
                   gx=xgrid[aix],
                   gy=ygrid[aix],
                   ox=VecX,
                   oy=VecY,
                   Dh=argv$Dh,
                   eps2=argv$eps2)
    t11<-Sys.time()
    if (argv$verbose) print(paste("OI_RR, time=",round(t11-t00,1),attr(t11-t00,"unit")))
    t00<-Sys.time()
    xa<-apply(cbind(res$xa,sqrt(abs(res$xa_errvar))),
              MARGIN=1,
              FUN=tboxcox4pdf_apply,
                  lambda=argv$transf.boxcox_lambda,
                  brrinf=boxcox(argv$rrinf,argv$transf.boxcox_lambda))
    t11<-Sys.time()
    if (argv$verbose) print(paste("backtransf, time=",round(t11-t00,1),attr(t11-t00,"unit")))
  } else {
    boom("transformation not defined")
  }
  ra<-rmaster
  ra[]<-NA
  ra[aix]<-xa
  rm(rfg,xb,xa,aix,xgrid,ygrid)
  if (argv$loocv_mode) {
    W<-tcrossprod((D-argv$eps2*diag(n0)),InvD)
    if (argv$transf=="none") {
      ya<-OI_RR_fast(yo=yo,
                     yb=yb,
                     xb=yb,
                     xgrid=VecX,
                     ygrid=VecY,
                     VecX=VecX,
                     VecY=VecY,
                     Dh=argv$Dh)
      yav<-yo + 1./(1.-diag(W)) * (ya-yo)
    } else if (argv$transf=="Box-Cox") {
      ya<-apply(cbind(res$ya,sqrt(abs(res$ya_errvar))),
                MARGIN=1,
                FUN=tboxcox4pdf_apply,
                    lambda=argv$transf.boxcox_lambda,
                    brrinf=boxcox(argv$rrinf,argv$transf.boxcox_lambda))
      yavbc<-boxcox(yo,argv$transf.boxcox_lambda) + 
             1./(1.-diag(W)) * (res$ya-boxcox(yo,argv$transf.boxcox_lambda))
      # this is from Lussana et al (2010), Eq.(19) 
      yavbc_errvar<-res$o_errvar/argv$eps2*(diag(InvD)-res$o_errvar)
      yav<-apply(cbind(yavbc,
                       sqrt(abs(yavbc_errvar))),
                       MARGIN=1,
                       FUN=tboxcox4pdf_apply,
                           lambda=argv$transf.boxcox_lambda,
                           brrinf=boxcox(argv$rrinf,argv$transf.boxcox_lambda))
    }
  # loo-crossvalidation not required
  } else {
    if (argv$transf=="none") {
      ya<-extract(ra,cbind(VecX,VecY),method="bilinear")
    } else if (argv$transf=="Box-Cox") {
      ya<-apply(cbind(res$ya,sqrt(abs(res$ya_errvar))),
                MARGIN=1,
                FUN=tboxcox4pdf_apply,
                    lambda=argv$transf.boxcox_lambda,
                    brrinf=boxcox(argv$rrinf,argv$transf.boxcox_lambda))
    }
    yav<-rep(-9999,length(ya))
  }
  if (argv$idiv_instead_of_elev) {
    if (!exists("W")) W<-tcrossprod((D-argv$eps2*diag(n0)),InvD)
    # this is the cross-validation integral data influence ("yidiv")
    elev_for_verif<-rep(1,n0) + 1./(1.-diag(W)) * (rowSums(W)-rep(1,n0))
  } else {
    elev_for_verif<-VecZ
  }
# ===>  OI multiscale (without background)   <===
} else if (argv$mode=="OI_multiscale") {
  # multi-scale OI operates on relative anomalies (if rescaling factor is available)
  if (exists("yrf")) {
    if (any(yrf==0)) yrf[which(yrf==0)]<-1
    yo_relan<-yo/yrf
  #  arrinf<-mean(argv$rrinf/rf[mask])
  } else {
    yo_relan<-yo
  #  arrinf<-mean(argv$rrinf/rf[mask])
  }
  # multi-scale OI
  for (l in 1:nl) {
    if (argv$verbose) 
      print(paste("scale # ",formatC(l,width=2,flag="0"),
                  " of ",nl,
                  " (",formatC(vecd[l],width=4,flag="0"),"km ",
                  "fact=",formatC(vecf[l],width=3,flag="0"),")",sep=""))
    D<-exp(-0.5*(Disth/vecd[l])**2.)
    diag(D)<-diag(D)+vece[l]*ovarc
    InvD<-chol2inv(chol(D))
    if (l==nl | vecf[l]==1) {
      r<-rmaster
    } else {
      r<-aggregate(rmaster, fact=vecf[l], expand=T, na.rm=T)
    }
    zvalues.l<-getValues(r)
    storage.mode(zvalues.l)<-"numeric"
    xy.l<-xyFromCell(r,1:ncell(r))
    x.l<-sort(unique(xy.l[,1]))
    y.l<-sort(unique(xy.l[,2]),decreasing=T)
    mask.l<-which(!is.na(zvalues.l))
    zgrid.l<-zvalues.l[mask.l]
    xgrid.l<-xy.l[mask.l,1]
    ygrid.l<-xy.l[mask.l,2]
    rm(xy.l,zvalues.l)
    if (l==1) {
      yb<-rep(mean(yo_relan),length=n0)
      xb<-rep(mean(yo_relan),length=length(xgrid.l))
    } else {
      rb<-resample(ra,r,method="bilinear")
      xb<-getValues(rb)[mask.l]
      count<-0
      while (any(is.na(xb))) {
        count<-count+1
        buffer_length<-round(vecd[l]/(10-(count-1)),0)*1000
        if (!is.finite(buffer_length)) break 
        ib<-which(is.na(xb))
        aux<-extract(rb,cbind(xgrid.l[ib],ygrid.l[ib]),
                     na.rm=T,
                     buffer=buffer_length)
        for (ll in 1:length(aux)) xb[ib[ll]]<-mean(aux[[ll]],na.rm=T)
        rb[mask.l]<-xb
        rm(aux,ib)
      }
      yb<-extract(rb,cbind(VecX,VecY),method="bilinear")
#      debug_plots()
      rm(rb)
    }
    if (any(is.na(xb))) print("xb is NA")
    if (any(is.na(yb))) print("yb is NA")
    xa.l<-OI_RR_fast(yo=yo_relan,
                     yb=yb,
                     xb=xb,
                     xgrid=xgrid.l,
                     ygrid=ygrid.l,
                     VecX=VecX,
                     VecY=VecY,
                     Dh=vecd[l]) 
    xidiw.l<-OI_RR_fast(yo=yo[ixwet],
                        yb=rep(0,nwet),
                        xb=rep(0,length(xgrid.l)),
                        xgrid=xgrid.l,
                        ygrid=ygrid.l,
                        VecX=VecX[ixwet],
                        VecY=VecY[ixwet],
                        Dh=vecd[l]) 
    xidid.l<-OI_RR_fast(yo=yo[ixdry],
                        yb=rep(0,ndry),
                        xb=rep(0,length(xgrid.l)),
                        xgrid=xgrid.l,
                        ygrid=ygrid.l,
                        VecX=VecX[ixdry],
                        VecY=VecY[ixdry],
                        Dh=vecd[l]) 
    if (any(xidid.l>xidiw.l)) xa.l[which(xidid.l>xidiw.l)]<-0
  #  if (any(xidid.l<xidiw.l & xa.l<arrinf)) 
  #    xa.l[which(xidid.l<xidiw.l & xa.l<arrinf)]<-arrinf
    ra<-r
    ra[]<-NA
    ra[mask.l]<-xa.l
  } # end of multi-scale OI
  # back to precipitation values (from relative anomalies)
  if (exists("rf")) xa<-round(xa.l*rf[mask.l],2)
  ra[mask.l]<-xa
  rm(mask.l,xa)
  ya<-extract(ra,cbind(VecX,VecY),method="bilinear")
  yb<-rep(-9999,length(ya))
  yav<-rep(-9999,length(ya))
}
#
#------------------------------------------------------------------------------
# set dry regions to no-rain
jump<-T
if (!jump) {
  print("+---------------------------------------------------------------+")
  print("adjust for unlikely wet regions left by multi OI")
  Dh<-vecd[which(kseq==10 & !is.na(kseq))]
  D<-exp(-0.5*(Disth[ixwet,ixwet]/Dh)**2.)
  diag(D)<-diag(D)+ovarc[ixwet]
  InvD<-chol2inv(chol(D))
  xidiw.l<-OI_RR_fast(yo=rep(1,nwet),
                      yb=rep(0,nwet),
                      xb=rep(0,length(xgrid.l)),
                      xgrid=xgrid.l,
                      ygrid=ygrid.l,
                      VecX=VecX[ixwet],
                      VecY=VecY[ixwet],
                      Dh=Dh) 
  D<-exp(-0.5*(Disth[ixdry,ixdry]/Dh)**2.)
  diag(D)<-diag(D)+ovarc[ixdry]
  InvD<-chol2inv(chol(D))
  xidid.l<-OI_RR_fast(yo=rep(1,ndry),
                      yb=rep(0,ndry),
                      xb=rep(0,length(xgrid.l)),
                      xgrid=xgrid.l,
                      ygrid=ygrid.l,
                      VecX=VecX[ixdry],
                      VecY=VecY[ixdry],
                      Dh=Dh)
  ixd<-which(((0.6*xidid.l)>xidiw.l) & xa<1)
  if (length(ixd)>0) xa[ixd]<-0
  #debug_01_plots()
  ra[mask.l]<-xa
  rm(r,xb,yb,xa.l)
  #
  print("remove wet regions with no obs in them")
  xa_aux<-getValues(ra)
  xa_aux[which(!is.na(xa_aux) & xa_aux<argv$rrinf)]<-NA
  ra[]<-xa_aux
  rclump<-clump(ra)
  oclump<-extract(rclump,cbind(VecX[ixwet],VecY[ixwet]))
  fr<-freq(rclump)
  # remove clumps of YESprec cells less than (4x4)km^2 or not including wet obs
  ix<-which(!is.na(fr[,1]) & !is.na(fr[,2]) & ( (fr[,2]<=16) | !(fr[,1] %in% oclump)) )
  xa[which(getValues(rclump)[mask.l] %in% fr[ix,1])]<-0
  rm(xa_aux,rclump,oclump,fr,ix)
  #
  print("substitute dry regions with no obs in them with wet regions")
  ra[mask.l]<-xa
  xa_aux<-getValues(ra)
  xa_aux[which(!is.na(xa_aux) & xa_aux>=argv$rrinf)]<-NA
  xa_aux[which(!is.na(xa_aux))]<-1
  ra[]<-xa_aux
  rclump<-clump(ra)
  oclump<-extract(rclump,cbind(VecX[ixdry],VecY[ixdry]))
  fr<-freq(rclump)
  # remove clumps of YESprec cells less than (4x4)km^2 or not including wet obs
  ix<-which(!is.na(fr[,1]) & !is.na(fr[,2]) & ( (fr[,2]<=16) | !(fr[,1] %in% oclump)) )
  xa[which((getValues(rclump)[mask.l] %in% fr[ix,1]) & (xidiw.l>xidid.l))]<-argv$rrinf
  rm(xa_aux,rclump,oclump,fr,ix,xidiw.l,xidid.l)
  ## final field
}
# 
#CVmode
if (argv$cv_mode) {
  ya_cv<-extract(ra,cbind(VecX_cv,VecY_cv),method="bilinear")
  if (argv$idiv_instead_of_elev) {
    Disth<-matrix(ncol=ncv,nrow=ncv,data=0.)
    Disth<-(outer(VecY_cv,VecY_cv,FUN="-")**2.+
            outer(VecX_cv,VecX_cv,FUN="-")**2.)**0.5/1000.
    D<-exp(-0.5*(Disth/argv$Dh)**2.)
    diag(D)<-diag(D)+argv$eps2
    InvD<-chol2inv(chol(D))
    W<-tcrossprod((D-argv$eps2*diag(ncv)),InvD)
    # this is the cross-validation integral data influence ("yidiv")
    elev_for_verif_cv<-rep(1,ncv) + 1./(1.-diag(W)) * (rowSums(W)-rep(1,ncv))
  } else {
    elev_for_verif_cv<-VecZ_cv
  }
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if (argv$verbose) print("++ Output")
# Station Points
# CVmode
if (argv$cv_mode) {
  cat("# variable: Precipitation\n",file=argv$off_cvver,append=F)
  cat("# units: $mm$\n",file=argv$off_cvver,append=T)
  cat("unixtime     leadtime location  lat     lon      altitude obs      fcst\n",
      file=argv$off_cvver,append=T)
  ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(yav))
  cat(paste( as.numeric(as.POSIXct(argv$date_out, format=argv$date_out_fmt))," ",
             rep(0,length(VecS_cv))," ",
             VecS_cv," ",
             VecLat_cv," ",
             VecLon_cv," ",
             elev_for_verif_cv," ",
             round(yo_cv,2)," ",
             round(ya_cv,2),"\n",
             sep=""),
      file=argv$off_cvver,append=T)
  print(paste("data saved on file",argv$off_cvver))
  cat("date;sourceId;x;y;z;yo;yb;ya;yav;dqc;\n",
      file=argv$off_cvstn,append=F)
  cat(paste(argv$date_out,
            formatC(VecS_cv,format="f",digits=0),
            formatC(VecX_cv,format="f",digits=0),
            formatC(VecY_cv,format="f",digits=0),
            formatC(VecZ_cv,format="f",digits=0),
            formatC(yo_cv,format="f",digits=1),
            formatC(yb_cv,format="f",digits=1),
            formatC(ya_cv,format="f",digits=1),
            formatC(ya_cv,format="f",digits=1),
            rep(0,length(VecS_cv)),
            "\n",sep=";"),
      file=argv$off_cvstn,append=T)
  print(paste("data saved on file",argv$off_cvstn))
# non-CVmode
} else {
  cat("# variable: Precipitation\n",file=argv$off_ver,append=F)
  cat("# units: $mm$\n",file=argv$off_ver,append=T)
  cat("unixtime     leadtime location  lat     lon      altitude obs      fcst\n",
      file=argv$off_ver,append=T)
  if (argv$off_ver_fg!="none") {
    ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(ya) & !is.na(yb))
  } else {
    ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(ya))
  }
  cat(paste( as.numeric(as.POSIXct(argv$date_out, format=argv$date_out_fmt))," ",
             rep(0,length(ix))," ",
             VecS[ix]," ",
             VecLat[ix]," ",
             VecLon[ix]," ",
             elev_for_verif[ix]," ",
             round(yo[ix],1)," ",
             round(ya[ix],1),"\n",
             sep=""),
      file=argv$off_ver,append=T)
  print(paste("data saved on file",argv$off_ver))
  if (argv$off_ver_fg!="none") {
    cat("# variable: Precipitation\n",file=argv$off_ver_fg,append=F)
    cat("# units: $mm$\n",file=argv$off_ver_fg,append=T)
    cat("unixtime     leadtime location  lat     lon      altitude obs      fcst\n",
        file=argv$off_ver_fg,append=T)
    cat(paste( as.numeric(as.POSIXct(argv$date_out, format=argv$date_out_fmt))," ",
               rep(0,length(ix))," ",
               VecS[ix]," ",
               VecLat[ix]," ",
               VecLon[ix]," ",
               elev_for_verif[ix]," ",
               round(yo[ix],1)," ",
               round(yb[ix],1),"\n",
               sep=""),
        file=argv$off_ver_fg,append=T)
    print(paste("data saved on file",argv$off_ver_fg))
  }
  cat("date;sourceId;x;y;z;yo;yb;ya;yav;dqc;\n",
      file=argv$off_stn,append=F)
  digxy<-0
  if (argv$grid_master.proj4==proj4.llwgs84) digxy<-6
  cat(paste(argv$date_out,
            formatC(VecS,format="f",digits=0),
            formatC(VecX,format="f",digits=digxy),
            formatC(VecY,format="f",digits=digxy),
            formatC(VecZ,format="f",digits=0),
            formatC(yo,format="f",digits=1),
            formatC(yb,format="f",digits=1),
            formatC(ya,format="f",digits=1),
            formatC(yav,format="f",digits=1),
            formatC(ydqc.flag,format="f",digits=0),
            "\n",sep=";"),
      file=argv$off_stn,append=T)
  print(paste("data saved on file",argv$off_stn))
  # grid
  r.list<-list()
  r.list[[1]]<-matrix(data=getValues(ra),
                      ncol=length(y),
                      nrow=length(x))
  # define time for output
  tstamp_nc<-format(strptime(argv$date_out,argv$date_out_fmt),
                    format="%Y%m%d%H%M",tz="GMT")
  time_bnds<-array(format(rev(seq(strptime(argv$date_out,argv$date_out_fmt),
                                           length=2,by=argv$time_bnds_string)),
                   format="%Y%m%d%H%M",tz="GMT"),dim=c(1,2))
  out<-write_dotnc(grid.list=r.list,
                   file.name=argv$off_grd,
                   grid.type=argv$off_grd.grid,
                   x=x,
                   y=y,
                   var.name=argv$off_grd.varname,
                   var.longname=argv$off_grd.varlongname,
                   var.standardname=argv$off_grd.varstandardname,
                   var.version=argv$off_grd.varversion,
                   var.unit=argv$off_grd.varunit,
                   times=tstamp_nc,
                   times.unit=argv$off_grd.timesunit,
                   reference=argv$off_grd.reference,
                   proj4.string=argv$grid_master.proj4,
                   lonlat.out=argv$off_grd.write_lonlat,
                   round.dig=argv$off_grd.diground,
                   summary=argv$off_grd.summary,
                   source.string=argv$off_grd.sourcestring,
                   title=argv$off_grd.title,
                   comment=argv$off_grd.comment,
                   atts.var.add=NULL,
                   var.cell_methods=argv$off_grd.cell_methods,
                   time_bnds=time_bnds,
                   cf_1.7=T)
  print(paste("data saved on file",argv$off_grd))
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# normal exit
t1<-Sys.time()
print(paste("Normal exit, time=",round(t1-t0,1),attr(t1-t0,"unit")))
quit(status=0)
