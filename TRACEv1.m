function [Canth,PrefProps,Ages,Uncertainty]=TRACEv1(OutputCoordinates,Dates,...
    PredictorMeasurements,PredictorTypes,AtmCO2Trajectory,...
            varargin)                                     % Optional inputs
%  TRansit-time-distribution-based Anthropogenic Carbon Estimation (TRACE)
%  Version 1.0 (Version at time of submission)
%  
%  This code generates etimates of ocean anthropogenic carbon content from
%  user-supplied inputs of coordinates (lat, lon, depth), salinity,
%  temperature, and date. Information is also needed about the historical
%  and/or future CO2 trajectory.  This information can be provided or
%  default values can be asssumed.  The approach takes several steps:
%
%  Steps taken during the formulation of this code: 
%
%  (1) transit time distributions are fit to CFC-11, CFC-12, and SF6
%  measurements made since ~2000 c.e. predominantly on repeat hydrography
%  cruises. These transit time distributions are inverse gaussians with one
%  scale factor that sets the mean age and width of the distribution. (2) A
%  neural network is constructed that relates the scale factor to
%  coordinate information and measured salinity and temperature. (3) A
%  separate neural network is uses an identical approach to relate the same
%  information to preformed property distributions estimated in a previous
%  research effort (Carter et al. 2021a:Preformed Properties for Marine
%  Organic Matter and Carbonate Mineral Cycling Quantification)
%
%  Steps taken when the code is called: (1) Both neural networks are
%  called, generating information that is used herein to constuct preformed
%  property information along with a transit time distribution.  (2) The
%  user-provided dates are then combined with default, built-in, or
%  user-provided atmospheric CO2 trajectories to produce records of the
%  atmospheric CO2 exposure for each parcel of water. (3) This information
%  is combined with estimated equilibrium values for the given CO2 exposure
%  along with preindustrial (i.e., 280 uatm pCO2) CO2 exposure.  (4) The
%  difference between the two is computed as the anthropogenic carbon
%  estimate.
%
%  Updated 2024.08.16
%
%  Citation information: 
%  TRACE, ESPER_SF_NN, and ESPER_PP_NN Carter et al. submitted
%
%  Related citations: (related to seawater property estimation)
%  ESPER_LIR and ESPER_NN, Carter et al., 2021: https://doi.org/10.1002/lom3.10461
%  LIARv1*: Carter et al., 2016, doi: 10.1002/lom3.10087
%  LIARv2*, LIPHR*, LINR* citation: Carter et al., 2018, https://doi.org/10.1002/lom3.10232
%  LIPR*, LISIR*, LIOR*, first described/used: Carter et al., 2021, https://doi.org/10.1029/2020GB006623
%  * deprecated in favor of ESPERs
%
%  ESPER_NN and TRACE are inspired by CANYON-B, which also uses neural networks: 
%  Bittig et al. 2018: https://doi.org/10.3389/fmars.2018.00328
%
%  This function needs the CSIRO seawater package to run.  Scale
%  differences from TEOS-10 are a negligible component of estimate errors
%  as formulated.
% 
% *************************************************************************
% Input/Output dimensions:
% .........................................................................
% p: Integer number of desired  estimate types. For TRACE this is always 1.
% n: Integer number of desired estimate locations
% e: Integer number of equations used at each location (up to 4)
% y: Integer number of parameter measurement types provided by the user.
% n*e: Total number of estimates returned as an n by e array
% 
% *************************************************************************
% Required Inputs:
%
% OutputCoordinates (required n by 3 array): 
    % Coordinates at which estimates are desired.  The first column should
    % be longitude (degrees E), the second column should be latitude
    % (degrees N), and the third column should be depth (m).
% 
% Dates (required n by 1 array):
    % Year c.e.  Decimals are allowed, but the information contained in the
    % decimal is disregarded because the technique uncertainties are too
    % great to resolve the impacts of fractional years. 
    % 
% PredictorMeasurements (required n by y array): 
    % Parameter measurements that will be used to estimate alkalinity.  The
    % column order (y columns) is arbitrary, but specified by
    % PredictorTypes. Temperature should be expressed as degrees C and
    % salinity should be specified with the unitless convention.  NaN
    % inputs are acceptable, but will lead to NaN estimates for any
    % equations that depend on that parameter.  If temperature is not
    % provided then it will be estimated from salinity (not recommended).
    % 
% PredictorTypes (required 1 by y vector): 
    % Vector indicating which parameter is placed in each column of the
    % 'PredictorMeasurements' input.  Note that salinity is required for
    % all equations.  
    % 
    % Input Parameter Key: 
    % 1. Salinity
    % 2. Temperature
    %
% AtmCO2Trajectory (integer between 1 and 9):
    % This specifies the atmospheric CO2Trajectory.
    % There are several options between Historical and SSPs.
    % Historical are derived from this website:
    % https://www.ncei.noaa.gov/access/paleo-search/study/10437
    % supplemented from here after 1959
    % https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_annmean_mlo.txt
    % Historical is currently updated through 2022 and assumes a linear
    % increase in the years that follow (recommended that you update the
    % underlying data file as needed).
    % SSPs are from the supplement to Meinshausen et al. 2020:
    % https://gmd.copernicus.org/articles/13/3571/2020/
    % RCPs are not currently supported, but can be added to the file if the
    % user desires
    %
    % Options:
    % 1. Historical/Linear (modify historical CO2 file for updates)
    % 2. SSP1_1.9	
    % 3. SSP1_2.6	
    % 4. SSP2_4.5	
    % 5. SSP3_7.0	
    % 6. SSP3_7.0_lowNTCF	
    % 7. SSP4_3.4	
    % 8. SSP4_6.0	
    % 9. SSP5_3.4_over
%
% *************************************************************************
% Optional inputs:  All remaining inputs must be specified as sequential
% input argument pairs (e.g. "..., 'Equations',[1:4], etc...")
%
% PreindustrialxCO2 (optional 1 by 1 scalar, default 280 ppm): 
    % Allows the user to specify an arbitrary reference xCO2 for the
    % preindustrial era in units of ppm.  Higher numbers will result in
    % lower Canth estimates including negative values.  Small negative
    % estimates are also possible with the default input.
    %
% Equations (optional 1 by e vector, default [1]):
    % Vector indicating which predictors will be used to estimate
    % properties. In TRACEv1, this input should always be omitted because
    % there is only one possible equation, but this functionality is
    % retained in the code to allow for eventual generalization of the
    % TRACE NN to operate with more predictor combinations.
    % 
    % (S=salinity, Theta=potential temperature, N=nitrate, Si=silicate,
    % T=temperature, O2=dissolved oxygen molecule... see
    % 'PredictorMeasurements' for units).
    % ...........................................................
    % 1.  S, T,
    %
% MeasUncerts (Optional n by y array or 1 by y vector, default: [0.003 S,
    % 0.003 degrees C (T or theta), 1% AOU]: Array of measurement
    % uncertainties (see 'PredictorMeasurements' for units). Uncertainties
    % should be presented in order indicated by PredictorTypes. Providing
    % these estimates will improve estimate uncertainties in
    % 'UncertaintyEstimates'. Measurement uncertainties are a small part of
    % TRACE estimate uncertainties for WOCE-quality measurements. However,
    % estimate uncertainty scales with measurement uncertainty, so it is
    % recommended that measurement uncertainties be specified for sensor
    % measurements.  If this optional input argument is not provided, the
    % default WOCE-quality uncertainty is assumed.  If a 1 by y array is
    % provided then the uncertainty estimates are assumed to apply
    % uniformly to all input parameter measurements.
    %
% PerKgSwTF (Optional boolean, default true): 
    % Retained for future development (allowing for flexible units for
    % currently-unsupported predictors).
    %
% VerboseTF (Optional boolean, default true): 
    % Setting this to false will make TRACE stop printing updates to
    % the command line.  This behavior can also be permanently disabled
    % below. Warnings and errors, if any, will be given regardless.
    %
% *************************************************************************
% Outputs:
% 
% OutputEstimates: 
    % A n by e array of TRACE estimates specific to the coordinates and
    % parameter measurements provided as inputs.  Units are micromoles per
    % kg.
	%
% UncertaintyEstimates: 
    % A n by e array of uncertainty estimates specific to the coordinates,
    % parameter measurements, and parameter uncertainties provided.
    % 
% *************************************************************************
% Missing data: should be indicated with a NaN.  A NaN coordinate will
% yield NaN estimates for all equations at that coordinate.  A NaN
% parameter value will yield NaN estimates for all equations that require
% that parameter.
% 
% *************************************************************************
% Please send questions or related requests to brendan.carter@noaa.gov or
% brendan.carter@gmail.com.
% *************************************************************************
% Parsing inputs, setting defaults, and sanity checking inputs.
%
% Starting timer
tic

% Determining whether the user requested command-line update text.
a=strcmpi(varargin,'VerboseTF');
if any(a)
    VerboseTF=varargin{1,logical([0 a(1:end-1)])};
else
    VerboseTF=true;
end
% Uncomment following line beginning with "VerboseTF" and save the function
% if you want less command line spam and you don't want to have to keep
% telling the code to be quiet.
% VerboseTF=false;

% Verifying required inputs are provided
if nargin<4; error('TRACE called with too few input arguments.'); end

% Checking whether specific equations are specified.
a=strcmpi(varargin,'Equations');
if any(a)
    Equations=varargin{1,logical([0 a(1:end-1)])};
else
    Equations=1;
end
% Making [] argument for Equations equivalent to no argument.
if isempty(Equations);Equations=1; end
% Making 0 argument for Equations equivalent to no argument.
if Equations==0; Equations=1; end

% Checking for PerKgSwTF input and setting default if not given. Note, this
% input is not needed for TRACE, currently
a=strcmpi(varargin,'PerKgSwTF');
if any(a)
    PerKgSwTF=varargin{1,logical([0 a(1:end-1)])};
else
    PerKgSwTF=true;
end

% Checking for PreindustrialxCO2 input and setting default if not given
a=strcmpi(varargin,'PreindustrialxCO2');
if any(a)
    PreindustrialxCO2=varargin{1,logical([0 a(1:end-1)])};
else
    PreindustrialxCO2=280;                                                 % ppm
end

% Checking for MeasUncerts input and setting default if not given
a=strcmpi(varargin,'MeasUncerts');
if any(a)
    InputU=varargin{1,logical([0 a(1:end-1)])};
    UseDefaultUncertainties=false;
    % Sanity checking the MeasUncerts argument.  This also deals with the
    % possibility that the user has provided a single set of uncertainties
    % for all estimates.
    if not(max(size(MeasUncerts)==size(PredictorMeasurements))) && not(min(size(MeasUncerts)==size(PredictorMeasurements))) && not(max(size(MeasUncerts))==0)
        error('MeasUncerts must be undefined, a vector with the same number of elements as PredictorMeasurements has columns, [] (for default values), or an array of identical dimension to PredictorMeasurements.')
    elseif not(min(size(MeasUncerts)==size(PredictorMeasurements))) && not(max(size(MeasUncerts))==0)
        if ~(size(MeasUncerts,2)==size(PredictorMeasurements,2))
            error('There are different numbers of columns of input uncertainties and input measurements.')
        end
        InputU=ones(size(PredictorMeasurements(:,1)))*InputU;              % Copying uncertainty estimates for all estimates if only singular values are provided
    end
    if ~(size(PredictorTypes,2)==size(PredictorMeasurements,2))            % Making sure all provided predictors are identified.
        error('The PredictorTypes input does not have the same number of columns as the PredictorMeasurements input.  This means it is unclear which measurement is in which column.');
    end
else
    UseDefaultUncertainties=true;
end


% TRACE requires non-NaN coordinates to provide an estimate.  This step
% eliminates NaN coordinate combinations prior to estimation.  NaN
% estimates will be returned for these coordinates.
NaNGridCoords=max(isnan(OutputCoordinates),[],2) | isnan(Dates);

% This step checks for negative depths.  If found, it changes them to
% positive depths and issues a warning.
if min(OutputCoordinates(~NaNGridCoords,3))<0
    OutputCoordinates(~NaNGridCoords,3)=abs(OutputCoordinates(~NaNGridCoords,3));
    warning('TRACE: Negative depths were detected and changed to positive values.')
end

% Doing a size check for the coordinates.
if ~(size(OutputCoordinates,2)==3)
    error('OutputCoordinates has either too many or two few columns.  This version only allows 3 columns with the first being longitude (deg E), the second being latitude (deg N), and the third being depth (m).');
end

% Figuring out how many estimates are required
n=sum(~NaNGridCoords);

% Checking for common missing data indicator flags and warning if any are
% found.  Consider adding your commonly used flags here.
if max(ismember([-999 -9 -1*10^20],PredictorMeasurements))==1
    warning('TRACE: A common non-NaN missing data indicator (e.g. -999, -9, -1e20) was detected in the input measurements provided.  Missing data should be replaced with NaNs.  Otherwise, TRACE will interpret your inputs at face value and give terrible estimates.')
end

% Book-keeping with coordinate inputs and adjusting negative longitudes.
C=OutputCoordinates(~NaNGridCoords,:);
C(C(:,1)>360,1)=rem(C(C(:,1)>360,1),360);
C(C(:,1)<0,1)=rem(C(C(:,1)<0,1,1),360)+360;

% Throwing an error if latitudes are out of the expected range.
if max(abs(C(:,2)))>90
    error('TRACE: A latitude >90 degrees (N or S) has been detected.  Verify latitude is in the 2nd colum of the coordinate input.');
end

% Preparing full PredictorMeasurement uncertainty grid.  Some of this code
% is vestigial for TRACEv1
DefaultUncertainties=diag([1 1 0.02 0.02 0.02 0.01]);
DefaultUAll=zeros(size(PredictorMeasurements,1),6);
DefaultUAll(:,PredictorTypes)=PredictorMeasurements* ...
    DefaultUncertainties(PredictorTypes,PredictorTypes);                   % Setting multiplicative default uncertainties for P, N, O2, and Si.
DefaultUAll(:,ismember(PredictorTypes,[1 2]))=0.003;                       % Then setting additive default uncertainties for T and S.
DefaultUAll=DefaultUAll(~NaNGridCoords,:);
if UseDefaultUncertainties==false
    InputUAll=zeros(size(PredictorMeasurements));
    InputUAll(:,PredictorTypes)=InputU;
    InputUAll=max(cat(3,InputUAll, DefaultUAll),[],3);                     % Overriding user provided uncertainties that are smaller than the (minimum) default uncertainties
else
    InputUAll=DefaultUAll;
end  

% Making sure all provided predictors are identified.
if ~(size(PredictorTypes,2)==size(PredictorMeasurements,2))
    error('The PredictorTypes input does not have the same number of columns as the Measurements input.  This means it is unclear which measurement is in which column.');
end

% If temperature is not provided, estimate it from salinity
if ~ismember(2,PredictorTypes);
    Ests=ESPER_EstT_NN(1,OutputCoordinates,PredictorMeasurements(:,PredictorTypes==1),[1]);
    PredictorMeasurements=horzcat(PredictorMeasurements,Ests.Temperature);
    PredictorTypes=horzcat(PredictorTypes,2);
    warning('TRACE was called either without providing temperature or without specifying which column of PredictorMeasurements contains temperature.  Temperature is therefore being estimated from salinity and coordinate information, but this is not optimal and the validation for TRACE should not be considered appropriate for the estimates returned from this function call.')
end

% Putting all provided measurement inputs in standard order
MAll=NaN(n,2);                                                             % PredictorMeasurements
UAll=NaN(n,2);                                                             % Uncertainties
MAll(:,PredictorTypes)=PredictorMeasurements(~NaNGridCoords,:);            % Reordering and limiting to viable coordinates
UAll(:,PredictorTypes)=InputUAll(:,PredictorTypes);                        % This was already limited to viable coordinates for later use.
if size(Dates,1)==1
    Dates=Dates*ones(size(NaNGridCoords));                                 % Resizing if only a single date is provided
end
Dates=Dates(~NaNGridCoords,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remapping preformed property estimates using a neural network suite
[PrefPropsSub]=ESPER_PP_NN([1 2 4],C,MAll,[1 2],'VerboseTF','false');

% Remapping the scale factors using another neural network suite
[SFs]=ESPER_SF_NN(C,MAll,[1 2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading in the CO2 history.  Note, this history has been modified to
% reflect the values that would be expected in the surface ocean given the
% slow response of the surface ocean to a rapidly changing atmospheric
% value. "Adjusted" can be deleted in the following line to use the
% original atmospheric values.  If this approach is used, then users should
% consider altering CanthDiseq below to modulate the degree of equilibrium.
CO2Rec=load('CO2TrajectoriesAdjusted.txt');
CO2Rec=vertcat(CO2Rec(1,:),CO2Rec); % Constructing
CO2Rec(1,1)=-10^10; % Setting ancient CO2 to preindustrial placeholder

% The following line should eventually replace everything down to
% allocating section, but you need to verify that this has the expected
% behavior.  This will allow TRACE to run iteratively through 10000
% calculations.  Otherwise matrices start to get prohibitively large.

% [CanthSub]=CalcCanthFromSFsIGDiseq(Dates,horzcat(MAll(:,1:2),PrefPropsSub.PreformedSi,PrefPropsSub.PreformedP,PrefPropsSub.PreformedTA),horzcat(SFs.SFs(:,1,1),SFs.SFs(:,1,2)),CO2Set);

% Building the histograms
x=0.01:0.01:5;
pf=makedist('InverseGaussian','mu',1,'lambda',1/1.3); % lambda should perhaps be 1/1.3 from He et al Schwinger paper
y=pdf(pf,x);
NumVal=size(C,1);
Ventilation=y/sum(y);
SigVent=Ventilation>0.01; % Disregarding all ventilations that are less than 1% for computational efficiency
CO2Set=interp1(CO2Rec(:,1),CO2Rec(:,1+AtmCO2Trajectory),Dates-SFs.SFs(:,1,1)*[1:500])';
CO2Set=CO2Set'*Ventilation';
Age=(SFs.SFs(:,1,1)*[1:500])*Ventilation';

% Vapor correction term
VPWP = exp(24.4543 - 67.4509.*(100./(293.15+MAll(:,2))) - 4.8489.*log((293.15+MAll(:,2))./100));
VPCorrWP = exp(-0.000544.*MAll(:,1));
VPSWWP = VPWP.*VPCorrWP;
VPFac=1-VPSWWP;
%
% VPFac=1; % This overrides the commented code above
CanthDiseq=1; % This allows the user to arbitrarily change the degree of equilibration with the anthropogenic transient (recommended value of 1 if using adjusted CO2 trajectories, which is default).

Out=CO2SYSImpatient(PrefPropsSub.Preformed_TA,VPFac.*(CanthDiseq*(CO2Set-280)+280),1,4,MAll(:,1),MAll(:,2),MAll(:,2),0,0,PrefPropsSub.Preformed_Si,PrefPropsSub.Preformed_P,1,10,3); % Calculating equilibrium DIC with Canth
OutRef=CO2SYSImpatient(PrefPropsSub.Preformed_TA,PreindustrialxCO2*VPFac,1,4,MAll(:,1),MAll(:,2),MAll(:,2),0,0,PrefPropsSub.Preformed_Si,PrefPropsSub.Preformed_P,1,10,3); % Calculating equilibrium DIC with Preindustrial values
CanthSub=(Out(:,2)-OutRef(:,2)); % Removing disequilibrium term
% Preallocating
Ages=NaN(size(OutputCoordinates,1),1);
Canth=NaN(size(OutputCoordinates,1),1);
PrefProps.Preformed_TA=NaN(size(OutputCoordinates,1),1);
PrefProps.Preformed_Si=NaN(size(OutputCoordinates,1),1);
PrefProps.Preformed_P=NaN(size(OutputCoordinates,1),1);
% Allocating
Ages(~NaNGridCoords)=Age;
Canth(~NaNGridCoords)=CanthSub;
Uncertainty=(4.4^2+2^2+(0.15*Canth).^2).^(0.5);
PrefProps.Preformed_TA(~NaNGridCoords,1)=PrefPropsSub.Preformed_TA;
PrefProps.Preformed_Si(~NaNGridCoords,1)=PrefPropsSub.Preformed_Si;
PrefProps.Preformed_P(~NaNGridCoords,1)=PrefPropsSub.Preformed_P;
end