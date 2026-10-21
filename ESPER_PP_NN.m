function [Estimates]= ...
    ESPER_PP_NN(DesiredVariables,OutputCoordinates,PredictorMeasurements,PredictorTypes, ... % Required inputs
            varargin)                                     % Optional inputs
%  ESPER for Preformed Properties.  This approach is a computationally more
%  efficient form of generating interpolated preformed property estimates
%  than triangulation/interpolation, and has the advantage that it
%  transitions between preformed and observed properties within shallow
%  depths.
%  
%  Generally following ESPER_NN, though the estimated properties are
%  different and there is only one viable predictor combination (T and S).
%
%  Documentation and citations:
%  LIARv1: Carter et al., 2016, doi: 10.1002/lom3.10087
%  LIARv2, LIPHR, LINR citation: Carter et al., 2018, https://doi.org/10.1002/lom3.10232
%  LIPR, LISIR, LIOR, first described/used: Carter et al., 2021, https://doi.org/10.1029/2020GB006623
%  ESPER_PP and ESPER_NN, 2021.
%
%  ESPER_NN is inspired by CANYON-B, which also uses neural networks: 
%  Bittig et al. 2018: https://doi.org/10.3389/fmars.2018.00328
%
%  This function needs the CSIRO seawater package to run if measurements
%  are povided in molar units or if potential temperature or AOU are
%  needed but not provided by the user.  Scale differences from TEOS-10 are
%  a negligible component of  estimate error provided these
%  calculations are performed internally.
% 
% *************************************************************************
% Input/Output dimensions:
% .........................................................................
% p: Integer number of desired property estimate types (e.g., pTA, pN, pSi)
% n: Integer number of desired estimate locations
% e: Integer number of equations used at each location
% y: Integer number of parameter measurement types provided by the user.
% n*e: Total number of estimates returned as an n by e array
% 
% *************************************************************************
% Required Inputs:
%
% DesiredVariables:
    % Specifies which variables will be returned all outputs are in umol/kg
    % 1. pTA (preformed total titration seawater alkalinity)
    % 2. pP (preformed phosphate)
    % 3. pN (preformed nitrate)
    % 4. pSi (prefomred silicate)
    % 5. pO (preformed dissolved molecular oxygen,i.e., O2<aq>)
%
% OutputCoordinates (required n by 3 array): 
    % Coordinates at which estimates are desired.  The first column should
    % be longitude (degrees E), the second column should be latitude
    % (degrees N), and the third column should be depth (m).
% 
% PredictorMeasurements (required n by y array): 
    % Parameter measurements that will be used to estimate alkalinity.  The
    % column order (y columns) is arbitrary, but specified by
    % PredictorTypes. Concentrations should be expressed as micromol per kg
    % seawater unless PerKgSwTF is set to false in which case they should
    % be expressed as micromol per L, temperature should be expressed as
    % degrees C, and salinity should be specified with the unitless
    % convention.  NaN inputs are acceptable, but will lead to NaN
    % estimates for any equations that depend on that parameter.
    % 
% PredictorTypes (required 1 by y vector): 
    % Vector indicating which parameter is placed in each column of the
    % 'PredictorMeasurements' input.  Note that salinity is required for
    % all equations.  If O2 is provided, then temperature or potential
    % temperature must also be provided to convert O2 to AOU. For example,
    % if the first three columns of 'PredictorMeasurements' contain
    % salinity, silicate, and temperature, then PredictorTypes should equal
    % [1 5 7].
    % 
    % Input Parameter Key: 
    % 1. Salinity
    % 2. Temperature
%
% *************************************************************************
% Optional inputs:  All remaining inputs must be specified as sequential
% input argument pairs (e.g. "..., 'Equations',[1:16], 'OAAdjustTF', false,
% etc...")
%
% Equations (optional 1 by e vector, default []):
    % Vector indicating which equations will be used to estimate property.
    % There is only one option for ESPER_PP_NN, so this output should be
    % omitted. Some parts of the associated code are retained to throw an
    % error when an incorrect request is given to this code
    % 
    % (S=salinity, T=temperature, O2=dissolved oxygen molecule... see
    % 'PredictorMeasurements' for units).
    % ...........................................................
    % Output Equation Key (See below for explanation of A, B, and C):
    % 1.  S, T
    %
% MeasUncerts (Optional n by y array or 1 by y vector, default: [0.003 S,
    % 0.003 degrees C (T or theta), 1% P, 1% AOU or O2, 1% Si]: Array of
    % measurement uncertainties (see 'PredictorMeasurements' for units).
    % Uncertainties should be presented in order indicated by
    % PredictorTypes. Providing these estimates will improve ESPER_PP
    % estimate uncertainties in 'UncertaintyEstimates'. Measurement
    % uncertainties are a small part of ESPER_PP estimate uncertainties
    % for WOCE-quality measurements. However, estimate uncertainty scales
    % with measurement uncertainty, so it is recommended that measurement
    % uncertainties be specified for sensor measurements.  If this optional
    % input argument is not provided, the default WOCE-quality uncertainty
    % is assumed.  If a 1 by y array is provided then the uncertainty
    % estimates are assumed to apply uniformly to all input parameter
    % measurements.
    %
% PerKgSwTF (Optional boolean, default true): 
    % Many sensors provide measurements in micromol per L (molarity)
    % instead of micromol per kg seawater. Indicate false if provided
    % measurements are expressed in molar units (concentrations must be
    % micromol per L if so).  Outputs will remain in molal units
    % regardless.
    %
% VerboseTF (Optional boolean, default true): 
    % Setting this to false will make ESPER_PP stop printing updates to
    % the command line.  This behavior can also be permanently disabled
    % below. Warnings and errors, if any, will be given regardless.
    %
% *************************************************************************
% Outputs:
% 
% OutputEstimates: 
    % A n by e array of ESPER_PP estimates specific to the coordinates and
    % parameter measurements provided as inputs.  Units are micromoles per
    % kg (equivalent to the deprecated microeq per kg seawater).
	%
% UncertaintyEstimates: 
    % Uncertainties are not estimated for preformed properties. See the
    % preformed properties manuscript for uncertainty estimation for the
    % underlying (unmapped) values.  The uncertainties on the underlying
    % information should be larger than the uncertainty on the remapping
    % except for regions where the underlying information is not provided.
    % Users should use these routines in marginal seas at their own peril.
    %
% *************************************************************************
% Missing data: should be indicated with a NaN.  A NaN coordinate will
% yield NaN estimates for all equations at that coordinate.  A NaN
% parameter value will yield NaN estimates for all equations that require
% that parameter.
% 
% *************************************************************************
% Please send questions or related requests to brc.oceans@gmail.com.
% *************************************************************************

% Parsing inputs, setting defaults, and sanity checking inputs.

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
if nargin<4; error('ESPER_PP called with too few input arguments.'); end

% Checking whether specific equations are specified.
a=strcmpi(varargin,'Equations');
if any(a)
    Equations=varargin{1,logical([0 a(1:end-1)])};
else
    Equations=1;
end
% Making [] argument for Equations equivalent to no argument.
if isempty(Equations);Equations=1:4; end
% Making 0 argument for Equations equivalent to no argument.
if Equations==0; Equations=1:4; end

% Checking for PerKgSwTF input and setting default if not given
a=strcmpi(varargin,'PerKgSwTF');
if any(a)
    PerKgSwTF=varargin{1,logical([0 a(1:end-1)])};
else
    PerKgSwTF=true;
end

% ESPER_PP requires non-NaN coordinates to provide an estimate.  This step
% eliminates NaN coordinate combinations prior to estimation.  NaN
% estimates will be returned for these coordinates.
NaNGridCoords=max(isnan(OutputCoordinates),[],2);

% Doing a size check for the coordinates.
if ~(size(OutputCoordinates,2)==3)
    error('OutputCoordinates has either too many or two few columns.  This version only allows 3 columns with the first being longitude (deg E), the second being latitude (deg N), and the third being depth (m).');
end

% end input parsing
% *************************************************************************
% Figuring out how many estimates are required
n=sum(~NaNGridCoords);
e=size(Equations,2);
p=size(DesiredVariables,2);

% Checking for common missing data indicator flags and warning if any are
% found.  Consider adding your commonly used flags here.
if max(ismember([-999 -9 -1*10^20],PredictorMeasurements))==1
    warning('ESPER_PP: A common non-NaN missing data indicator (e.g. -999, -9, -1e20) was detected in the input measurements provided.  Missing data should be replaced with NaNs.  Otherwise, ESPER_PP will interpret your inputs at face value and give terrible estimates.')
end

% Book-keeping with coordinate inputs and adjusting negative longitudes.
C=OutputCoordinates(~NaNGridCoords,:);
C(C(:,1)>360,1)=rem(C(C(:,1)>360,1),360);
C(C(:,1)<0,1)=rem(C(C(:,1)<0,1,1),360)+360;

% Throwing an error if latitudes are out of the expected range.
if max(abs(C(:,2)))>90
    error('ESPER_PP: A latitude >90 degrees (N or S) has been detected.  Verify latitude is in the 2nd colum of the coordinate input.');
end

% Making sure all provided predictors are identified.
if ~(size(PredictorTypes,2)==size(PredictorMeasurements,2))
    error('The PredictorTypes input does not have the same number of columns as the Measurements input.  This means it is unclear which measurement is in which column.');
end

% Putting all provided measurement inputs in standard order
MAll=NaN(n,6);                                                             % PredictorMeasurements
% UAll=NaN(n,6);                                                             % Uncertainties
MAll(:,PredictorTypes)=PredictorMeasurements(~NaNGridCoords,:);            % Reordering and limiting to viable coordinates
% UAll(:,PredictorTypes)=InputUAll(:,PredictorTypes);                        % This was already limited to viable coordinates for later use.
YouHaveBeenWarnedCanth=false;                                              % Calculating Canth is slow, so this flag is used to make sure it only happens once.

% Beginning the iterations through the requested properties
for PIter=1:p                                                              
    Property=DesiredVariables(1,PIter);
    % Specifying which variables are required for each property estimate.
    % Each row is for a different property.
    NeededForProperty=[1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5];
    % Which of the 5 predictors are required for the equations specified?
    % Each row is for a different equation.
    VarVec=logical([1 1 1 0 0]); % depth sal and temp in NN
    NeedVars=any(VarVec(Equations,:),1);
    HaveVars=false(1,6);  HaveVars(1,PredictorTypes+1)=1; % Adding one because depth is provided
    HaveVars(1,1)=1; % This is depth in this version
    
    % Temperature is required if O2 is used as a predictor (to convert to
    % AOU)... equation 15 is therefore only valid under some circumstances.
    if ismember(6,PredictorTypes)==1; NeedVars(1,2)=true; end
    
    % Ignoring variables that aren't used for the current property
    HaveVarsSub=HaveVars(1,NeededForProperty(Property,:));

    % Temperature is required if measurements are provided in molar units
    % or if CO2sys will be used to adjust pH values
    if PerKgSwTF==false; NeedVars(1,2)=1; end 

    % Making sure all needed variables are present
    if ~all(HaveVarsSub(1,NeedVars)) && VerboseTF==true  % We are assuming only power users turn VerboseTF off, and hopefully they understand the function call
        disp('Warning: One or more regression equations for the current property require one or more input parameters that are either not provided or not labeled correctly. These equations will return NaN for all estimates.  All 16 equations are used by default unless specific equations are specified.  Temperature is also required when density or carbonate system calculations are called.'); % They should already be NaN
    end
    
    % Limiting measurements to the subset needed for this property estimate
    M=MAll(:,NeededForProperty(Property,:));
%     U=UAll(:,NeededForProperty(Property,:));
%     DefaultU=DefaultUAll(:,NeededForProperty(Property,:));

    % Checking to see whether temperature is needed, and subbing in
    % potential temp if it is.  As of the version 3 update the  routines
    % now request temperature, but convert it to potential temperature
    % before using it for predictions.
    if NeedVars(1,2)
        M(:,2)=sw_ptmp(M(:,1),M(:,2),sw_pres(C(:,3),C(:,2)),0);
    end
    % Checking to see whether O2 is needed. Defining AOU and subbing in for
    % O2 if yes (see above).
    if any(ismember(NeededForProperty(Property,:),6))
        M(:,4)=sw_satO2(M(:,1),M(:,2))*44.64-M(:,4);
    end
    % Converting units to molality if they are provided as molarity.
    if PerKgSwTF==false
        densities=sw_dens(M(:,1),M(:,2),sw_pres(C(:,3),C(:,2)))/1000;
        M(:,3)=M(:,3)./densities;
        M(:,4)=M(:,4)./densities;
        M(:,5)=M(:,5)./densities;
    end

    % *********************************************************************
    % Beginning treatment of inputs and calculations
    if     Property==1; VName='Preformed_TA';
    elseif Property==2; VName='Preformed_P';
    elseif Property==3; VName='Preformed_N';
    elseif Property==4; VName='Preformed_Si';
    elseif Property==5; VName='Preformed_O';
    else; error('A property identifier >8 or <1 was supplied, but this routine only has 2 possible property estimates.  The property identifier is the first input.')
    end
    
    % Loading the data, with an error message if not found
    FN='Polys.mat';
    % Making sure you downloaded the needed file and put it somewhere it
    % can be found
    if exist(FN,'file')<2; error('ESPER_PP could not find the file(s) needed to run.  These mandatory file(s) should be distributed from the same website as ESPER_PP.  Contact the corresponding author if you cannot find it there.  If you do have it then make sure all of the contents of the ESPER.zip extract are on the MATLAB path or in the active directory.  This will require adding several subfolders for ESPER.'); end
    L=load(FN);
  
    
    %Preallocating for speed
    OutputEstimates=NaN(size(OutputCoordinates,1),e);                      % Using size instead of n since we want to preallocate for NaN coordinate locations as well   
    Est=NaN(n,e);

    % Some of the equations use depth [m] as a predictor so this appends it
    % as a predictor.  This, so far, is unique to TRACE.
    M=horzcat(C(:,3),M);  


    % Disambiguation:
    % Eq... a counter for which equation ESPER_PP is on
    % e... the number of equations that will be used
    % Equation... the specific equation number currently being used
    % Equations... the user provided list of equations that will be used.
    for Eq=1:e                                                             % Iterating over the (up to 16) equations
        Equation=Equations(1,Eq);
        NeedVarsNow=any(VarVec(Equation,:),1);
        if all(HaveVarsSub(1,NeedVarsNow))                                 % Skipping if we don't have what we need
            P=horzcat(cosd(C(:,1)-20),sind(C(:,1)-20),C(:,2),M(:,VarVec(Equations(Eq),:)));
            for Net=1:4 % A committee of 4 neural networks is used.
                % Separate neural networks are used for the Arctic/Atlantic and
                % the rest of the ocean (more info below).  This functional
                % form of ESPER_NN is designed to avoid calling the nerual
                % network toolbox.  It will require the (many) ESPER_NN
                % functions to be on the MATLAB path.
                FunctionOther=str2func(horzcat('TRACE_',VName,'_',num2str(Equation),'_Other_',num2str(Net))); % Basically an eval statement.  Sorry Matlab power users.
                FunctionAtl=str2func(horzcat('TRACE_',VName,'_',num2str(Equation),'_Atl_',num2str(Net)));
                EstAtl(:,Eq,Net)=FunctionAtl(P')';
                EstOther(:,Eq,Net)=FunctionOther(P')';
            end
        end
    end
    % Averaging across neural network committee members
    EstAtl=nanmean(EstAtl,3);
    EstOther=nanmean(EstOther,3);

    % We do not want information to propagate across the Panama Canal (for
    % instance), so data is carved into two segments... the Atlantic/Arctic
    % (AA) and everything else.
    AAInds=or(inpolygon(C(:,1),C(:,2),L.Polys.LNAPoly(:,1),L.Polys.LNAPoly(:,2)), ...
        or(inpolygon(C(:,1),C(:,2),L.Polys.LSAPoly(:,1),L.Polys.LSAPoly(:,2)), ...
        or(inpolygon(C(:,1),C(:,2),L.Polys.LNAPolyExtra(:,1),L.Polys.LNAPolyExtra(:,2)), ...
        or(inpolygon(C(:,1),C(:,2),L.Polys.LSAPolyExtra(:,1),L.Polys.LSAPolyExtra(:,2)), ...
        inpolygon(C(:,1),C(:,2),L.Polys.LNOPoly(:,1),L.Polys.LNOPoly(:,2))))));
    % We'd like a smooth transition in the Bering Strait and in the South
    % Atlantic.  This linearly interpolates between the networks at their
    % boundaries.
    BeringInds=inpolygon(C(:,1),C(:,2),L.Polys.Bering(:,1),L.Polys.Bering(:,2));
    SoAtlInds=(C(:,1)>290 | C(:,1)<20) & (C(:,2)>-44) & (C(:,2)<-34);
    SoAfrInds=(C(:,1)>19 & C(:,1)<27) & (C(:,2)>-44) & (C(:,2)<-34);
    Est=EstOther;                                                          % Pulling out Indo-Pacific estimates
    Est(AAInds,:)=EstAtl(AAInds,:,:);                                      % Pulling out all other estimates
    Est(BeringInds,:)=EstAtl(BeringInds,:).* ...                           % Blending estimates in the Bering Strait vicinity
        repmat((C(BeringInds,2)-62.5)/(7.5),[1,size(Equations,2)]) ...
        +EstOther(BeringInds,:).* ...
        repmat((70-C(BeringInds,2))/(7.5),[1,size(Equations,2)]);
    Est(SoAtlInds,:)=EstAtl(SoAtlInds,:).* ...                             % Blending estimates in the Southern Atlantic to Southern Ocean transition
        repmat((C(SoAtlInds,2)+44)/(10),[1,size(Equations,2)])+ ...
        EstOther(SoAtlInds,:).* ...
        repmat((-34-C(SoAtlInds,2))/(10),[1,size(Equations,2)]);
    Est(SoAfrInds,:)=Est(SoAfrInds,:).* ...                                % Blending estimates South of South Africa toward Indo-Pac Estimates
        repmat((27-C(SoAfrInds,1))/(8),[1,size(Equations,2)])+ ...
        EstOther(SoAfrInds,:).* ...
        repmat((C(SoAfrInds,1)-19)/(8),[1,size(Equations,2)]);
    
    % Filling the outputs with the estimates at viable locations.
    OutputEstimates(~NaNGridCoords,:)=Est;
    Estimates.(VName)=OutputEstimates;
end