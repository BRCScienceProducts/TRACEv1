function [Y,Xf,Af] = TRACE_EstT_Temperature_1_Atl_4(X,~,~)
%TRACE_ESTT_TEMPERATURE_1_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 21-Aug-2024 10:13:41.
% 
% [Y] = TRACE_EstT_Temperature_1_Atl_4(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 5xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999961923064171;-0.999961923064171;-39.5;0;30.0000762939453];
x1_step1.gain = [1.00003807838574;1.00003807838574;0.0155038759689922;0.000363636363636364;0.206927881437043];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.0419886356995973209;2.2420934709461421974;-1.129637170316138528;-0.57423469379453773165;-1.1728266064269388114];
IW1_1 = [-0.0020046938928910732446 0.043096490870167050713 1.3996454060670586728 0.37939048928935681371 0.54722553218978209344;0.015929123487396477238 -0.0094254235567978583016 -0.13292936492990145569 0.50991599597213543849 -2.0043572016163961891;1.0513893261635411935 0.22085298208661852115 -3.3402504577177940526 -2.0601479422193764712 0.06873160888663795709;0.0071252190471554500278 -0.062226361586654614644 -1.3591350107633322075 -0.018003514740217493179 -0.45261584601920329396;-0.0061831941333087957896 -0.060477861034444024835 -1.3730034391238867819 -0.63593320044818257841 -0.45657218151418021979];

% Layer 2
b2 = [-12.149324813946176249;-52.230870116645498058;7.8358848404293759415;10.207900302548209126;-71.603823068536428309;-5.3785639188302205582;231.90419106902496083;13.132560966382749612;12.42540570656812271;-122.86796621059440326;-0.68685168625871739767;-3.1154909598785995861;0.65358749933637583229;-54.835798355053086084;133.43864901780401055];
LW2_1 = [-83.50350025258092046 14.256669167159111922 1.3430864577953876005 -37.29470539969435805 -52.57534175205331195;192.09936763820164174 10.049347181277868302 1.9445461465210425178 8.8638935789393897835 138.61089602535361109;1.8745150331066264826 -15.644342048436643111 5.4385709220681990317 -314.83825635548834043 305.92035186892951515;-53.977663837200474006 -11.276552303919757492 6.9431296807972708862 21.830819599674150311 -81.023541845147647678;18.258873306628025546 2.4182278891038833102 57.147715782392282335 16.402326192646594194 14.350379362568961739;96.216505711768192555 -73.637461778868200213 61.290437697421289442 -4.5859978837419550501 -48.700844197238907896;-495.34203435377833102 298.80361969253414145 -220.96996201998570086 -162.67678710119460561 354.4293380295051179;-0.30075210762732401681 -13.746784044524844148 0.0039019566001782684563 -0.76959385322126805384 0.6741820722220305262;2.9525728971969291692 -11.745509252477058126 -1.1087681766000077044 56.991462188859422611 -52.962075075735761231;-124.46830642949953472 96.808399450173084233 10.624027956819242391 -51.431189646301746166 -115.42643681023005797;-6.717700868198440034 0.75453568715616803342 0.12659718360410096571 18.781085272967068533 -25.819695630805334474;-59.238728765226390749 3.266432488281943769 -4.6807448756666776291 -4.6484030271656751054 -50.742935701893202349;6.7339192091635213799 -0.76780959540147020181 -0.14451125834598321784 -17.764572473193261715 24.861026403799307616;201.40588462515142965 43.758549904264086194 -8.117281293912538942 94.641311440326006732 97.068233114868561984;22.141894882484994156 -8.4706632581652776537 -13.516943659366431163 64.552430898072714172 108.69969960730514913];

% Layer 3
b3 = -13.50470928069444021;
LW3_2 = [-0.048324223424351282763 0.095119129571034938264 0.063019774601889497201 0.095091318449520456646 -14.172482377160273259 0.017480276492915034631 -0.0058230502888186918067 0.49877410992992798899 -0.58484821106798101642 0.049957499901467028458 8.3897029721558347148 0.034926938360461193267 8.8537334834704974895 0.025276528072659402296 0.08607732000683869511];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0596262898164022;
y1_step1.xoffset = -2.99982047080994;

% ===== SIMULATION ========

% Format Input Arguments
isCellX = iscell(X);
if ~isCellX
  X = {X};
end

% Dimensions
TS = size(X,2); % timesteps
if ~isempty(X)
  Q = size(X{1},2); % samples/series
else
  Q = 0;
end

% Allocate Outputs
Y = cell(1,TS);

% Time loop
for ts=1:TS

    % Input 1
    Xp1 = mapminmax_apply(X{1,ts},x1_step1);
    
    % Layer 1
    a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*Xp1);
    
    % Layer 2
    a2 = tansig_apply(repmat(b2,1,Q) + LW2_1*a1);
    
    % Layer 3
    a3 = repmat(b3,1,Q) + LW3_2*a2;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a3,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(3,0);

% Format Output Arguments
if ~isCellX
  Y = cell2mat(Y);
end
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end
