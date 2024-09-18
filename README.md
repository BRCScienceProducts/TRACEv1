# TRACEv1 (for MATLAB)
Tracer-based Rapid Anthropogenic Carbon Estimation (TRACE)
Version 1.0 (Version at time of manuscript submission, will be updated to 1.1 upon publication)

This code generates estimates of ocean anthropogenic carbon content from
user-supplied inputs of coordinates (lat, lon, depth), salinity,
temperature, and date. Information is also needed about the historical
and/or future CO2 trajectory.  This information can be provided or
default values can be assumed.  

**Setup:** download the TRACEv1 archive to your computer and extract to a location 
on your MATLAB path (or add the location to your path).  Ensure the TRACEv1.m 
function and the /private directory are in the same folder.  Call the function 
using the guidelines below (and in the function description).

A detailed description of the methodology will be provided in a yet-to-be submitted 
manuscript.  The following is a summary.  The approach takes several steps:

**Steps taken during the formulation of this code: **
(1) transit time distributions are fit to CFC-11, CFC-12, and SF6
measurements made since ~2000 c.e. predominantly on repeat hydrography
cruises. These transit time distributions are inverse gaussians with one
scale factor that sets the mean age and width of the distribution. (2) A
neural network is constructed that relates the scale factor to
coordinate information and measured salinity and temperature. (3) A
separate neural network is uses an identical approach to relate the same
information to preformed property distributions estimated in a previous
research effort (Carter et al. 2021a:Preformed Properties for Marine
Organic Matter and Carbonate Mineral Cycling Quantification)

**Steps taken when the code is called: **
(1) Both neural networks are called, generating information that is used 
herein to construct preformed property information along with a transit 
time distribution.  (2) The user-provided dates are then combined with 
default, built-in, or user-provided atmospheric CO2 trajectories to 
produce records of the atmospheric CO2 exposure for each parcel of water. 
(3) This information is combined with estimated equilibrium values for 
the given CO2 exposure along with preindustrial (i.e., 280 uatm pCO2) 
CO2 exposure.  (4) The difference between the two is computed as the 
anthropogenic carbon estimate.

Updated 2024.08.16

Citation information: 
TRACE, ESPER_SF_NN, and ESPER_PP_NN Carter et al. submitted

Related citations: (related to seawater property estimation)
ESPER_LIR and ESPER_NN, Carter et al., 2021: https://doi.org/10.1002/lom3.10461
LIARv1*: Carter et al., 2016, doi: 10.1002/lom3.10087
LIARv2*, LIPHR*, LINR* citation: Carter et al., 2018, https://doi.org/10.1002/lom3.10232
LIPR*, LISIR*, LIOR*, first described/used: Carter et al., 2021, https://doi.org/10.1029/2020GB006623
* deprecated in favor of ESPERs

ESPER_NN and TRACE are inspired by CANYON-B, which also uses neural networks: 
Bittig et al. 2018: https://doi.org/10.3389/fmars.2018.00328

This function needs the CSIRO seawater package to run.  Scale
differences from TEOS-10 are a negligible component of estimate errors
as formulated.

*************************************************************************
Input/Output dimensions:
.........................................................................
p: Integer number of desired  estimate types. For TRACE this is always 1.
n: Integer number of desired estimate locations
e: Integer number of equations used at each location (up to 4)
y: Integer number of parameter measurement types provided by the user.
n*e: Total number of estimates returned as an n by e array

*************************************************************************
Required Inputs:
%
OutputCoordinates (required n by 3 array): 
    Coordinates at which estimates are desired.  The first column should
    be longitude (degrees E), the second column should be latitude
    (degrees N), and the third column should be depth (m).

Dates (required n by 1 array):
    Year c.e.  Decimals are allowed, but the information contained in the
    decimal is disregarded because the technique uncertainties are too
    great to resolve the impacts of fractional years. 
    
PredictorMeasurements (required n by y array): 
    Parameter measurements that will be used to estimate alkalinity.  The
    column order (y columns) is arbitrary, but specified by
    PredictorTypes. Temperature should be expressed as degrees C and
    salinity should be specified with the unitless convention.  NaN
    inputs are acceptable, but will lead to NaN estimates for any
    equations that depend on that parameter.  If temperature is not
    provided then it will be estimated from salinity (not recommended).
    
PredictorTypes (required 1 by y vector): 
    Vector indicating which parameter is placed in each column of the
    'PredictorMeasurements' input.  Note that salinity is required for
    all equations.  
    
    Input Parameter Key: 
    1. Salinity
    2. Temperature
   
AtmCO2Trajectory (integer between 1 and 9):
    This specifies the atmospheric CO2Trajectory.
    There are several options between Historical and SSPs.
    Historical are derived from this website:
    https://www.ncei.noaa.gov/access/paleo-search/study/10437
    supplemented from here after 1959
    https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_annmean_mlo.txt
    Historical is currently updated through 2022 and assumes a linear
    increase in the years that follow (recommended that you update the
    underlying data file as needed).
    SSPs are from the supplement to Meinshausen et al. 2020:
    https://gmd.copernicus.org/articles/13/3571/2020/
    RCPs are not currently supported, but can be added to the file if the
    user desires
   
    Options:
    1. Historical/Linear (modify historical CO2 file for updates)
    2. SSP1_1.9	
    3. SSP1_2.6	
    4. SSP2_4.5	
    5. SSP3_7.0	
    6. SSP3_7.0_lowNTCF	
    7. SSP4_3.4	
    8. SSP4_6.0	
    9. SSP5_3.4_over

*************************************************************************
Optional inputs:  All remaining inputs must be specified as sequential
input argument pairs (e.g. "..., 'Equations',[1:4], etc...")

PreindustrialxCO2 (optional 1 by 1 scalar, default 280 ppm): 
    Allows the user to specify an arbitrary reference xCO2 for the
    preindustrial era in units of ppm.  Higher numbers will result in
    lower Canth estimates including negative values.  Small negative
    estimates are also possible with the default input.
   
Equations (optional 1 by e vector, default [1]):
    Vector indicating which predictors will be used to estimate
    properties. In TRACEv1, this input should always be omitted because
    there is only one possible equation, but this functionality is
    retained in the code to allow for eventual generalization of the
    TRACE NN to operate with more predictor combinations.
    
    (S=salinity, Theta=potential temperature, N=nitrate, Si=silicate,
    T=temperature, O2=dissolved oxygen molecule... see
    'PredictorMeasurements' for units).
    ...........................................................
    1.  S, T,
   
MeasUncerts (Optional n by y array or 1 by y vector, default: [0.003 S,
    0.003 degrees C (T or theta), 1AOU]: Array of measurement
    uncertainties (see 'PredictorMeasurements' for units). Uncertainties
    should be presented in order indicated by PredictorTypes. Providing
    these estimates will improve estimate uncertainties in
    'UncertaintyEstimates'. Measurement uncertainties are a small part of
    TRACE estimate uncertainties for WOCE-quality measurements. However,
    estimate uncertainty scales with measurement uncertainty, so it is
    recommended that measurement uncertainties be specified for sensor
    measurements.  If this optional input argument is not provided, the
    default WOCE-quality uncertainty is assumed.  If a 1 by y array is
    provided then the uncertainty estimates are assumed to apply
    uniformly to all input parameter measurements.
   
PerKgSwTF (Optional boolean, default true): 
    Retained for future development (allowing for flexible units for
    currently-unsupported predictors).
   
VerboseTF (Optional boolean, default true): 
    Setting this to false will make TRACE stop printing updates to
    the command line.  This behavior can also be permanently disabled
    below. Warnings and errors, if any, will be given regardless.
   
*************************************************************************
Outputs:

OutputEstimates: 
    A n by e array of TRACE estimates specific to the coordinates and
    parameter measurements provided as inputs.  Units are micromoles per
    kg.
	
UncertaintyEstimates: 
    A n by e array of uncertainty estimates specific to the coordinates,
    parameter measurements, and parameter uncertainties provided.
    
*************************************************************************
Missing data: should be indicated with a NaN.  A NaN coordinate will
yield NaN estimates for all equations at that coordinate.  A NaN
parameter value will yield NaN estimates for all equations that require
that parameter.

*************************************************************************
Please send questions or related requests to brendan.carter@noaa.gov or
brendan.carter@gmail.com.
*************************************************************************

Example calls:

_This first example asks for estimates at the surface ocean at the equator/prime meridian in the years 2000 and 2200 assuming SSP5_3.4_over is followed_

[Canth]=TRACEv1([0 0 0;0 0 0],[2000;2200],[35 20;35 20],[1 2],[9],[0])

_Results in_
Canth =

   45.8414
   80.2198

_This second example demonstrates a function call performed without providing temperature information, which is not recommended and should result in a warning_

[Canth]=TRACEv1([0 0 0;0 0 0],[2000;2010],[35;35],[1],[1],[0])

_Results in_
Warning: TRACE was called either without providing temperature or without
specifying which column of PredictorMeasurements contains temperature.
Temperature is therefore being estimated from salinity and coordinate information,
but this is not optimal and the validation for TRACE should not be considered
appropriate for the estimates returned from this function call. 
> In TRACEv1 (line 318) 

Canth =

   53.8992
   64.0482
