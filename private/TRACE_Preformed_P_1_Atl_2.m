function [Y,Xf,Af] = TRACE_Preformed_P_1_Atl_2(X,~,~)
%TRACE_PREFORMED_P_1_ATL_2 neural network simulation function.
%
% Auto-generated by MATLAB, 21-Aug-2024 10:13:38.
% 
% [Y] = TRACE_Preformed_P_1_Atl_2(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 6xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999961923064171;-0.999961923064171;-39.5;0;30.0316562652588;-2.96613502502441];
x1_step1.gain = [1.00003807838574;1.00003807838574;0.0155038759689922;0.000363636363636364;0.207805673664695;0.0607333297500273];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.5900084516265899914;2.7705142645288769465;-3.5795497308531492031;-6.4584052410234829011;1.2272317813982314139;3.1869245871482121402;-1.4614131574981195083;-3.1541431565105537693;0.9882986362094554611;-2.437538219793318639];
IW1_1 = [0.013299705865955049194 -0.54860690332126482716 -1.3994716860100027134 -0.1157521285615405876 -6.8921044796861297144 3.1012598997860369998;0.025605131242928780982 0.041218934091922843077 2.4846813250895434244 0.39850807521917330778 7.2534438265918286248 0.63521205245387346583;-0.3979542229896530503 0.14967313351293159984 -0.5132442252979466657 -3.0255249000132224957 0.34223982297446348122 -0.19917765989759983469;5.2698177603902616184 -2.8294256385036127277 -5.3783737115899130643 2.8329907693352787135 -4.2652843264178947891 8.1275109421951441391;-0.1085463136873630885 -0.037467802696458017697 -0.6339651869931408612 0.92830123468880054727 1.899055052204864813 0.48669148042801879672;0.73585137432372382804 2.1836315897961715571 6.3519801713452759273 -1.0860989598662837974 0.6010672493982655018 -0.58481376638685411695;0.75082535695771035478 -0.30947290945755101532 -0.32039557820976344304 -1.2258849981232577964 -1.3085714225433928704 -0.65473879589837358228;0.55002683958515630191 -0.57472031188047689731 -0.86900562977272677312 0.57090232768711801636 -0.56953534953976281408 -2.4519837146574183784;-0.063198989319134132869 -0.03027571739257974015 2.3154428901059538326 0.99184408608039642541 -1.9642142716262716018 -0.87025387164792256112;-0.54228660937988226731 0.29314730704670322048 0.018609326178891369419 -2.957957490724945071 1.3845173451734962278 0.69595368287300307664];

% Layer 2
b2 = [-12.104403099692358836;2.36479057304469098;-4.043405170550206762;3.2895306354206255151;0.60155799649706009902;1.4985905350370758438;-0.74780311417127687257;0.66045403309488204258;-5.7002954488220369811;1.4368551270203355763];
LW2_1 = [-1.4073759226693305369 -0.33732178991373451549 -4.6457203485159537237 0.73941595154312556026 9.0822994075379188672 -0.87656182682905314785 1.9585707760595063576 0.16742465100295034208 0.83787385285791249334 0.32588744701265282;0.42194201532535230426 -0.30937287841163407798 5.1234194635721648226 -0.84808426295250605165 -2.314647299516239709 -0.024052856588538631954 -0.38918499726680411133 -0.49723740508993929943 -0.24589412082700060069 -2.4635022140645239652;-0.81876645793259617001 0.036042461345548004736 -4.8948586109583445136 0.019997786887099491121 -0.0030400848504947568414 0.35223114627575086155 -0.56601446780413566806 -0.57880856123443591432 -0.1890044145027393041 0.98592263045959738221;0.27327913418084498476 0.68714980857289731375 1.3875768457447261728 -0.56962451586780071811 0.83490458989146576485 0.29188726687924537728 -0.041623993412997009611 2.4792362803280765249 0.71147416655281603415 -0.83651403737073182221;0.59322731383354043722 -10.796514918064200828 -1.2845173284793744983 -1.4615985289119037827 -2.7314624285741513354 -1.9704853591291526804 -6.0676693515244402022 -14.813360936780346933 -4.3065180724089326958 -10.811105574801802121;-0.24197032690771561048 0.61881976355086976138 -2.419548198477535994 1.9972008409852533006 1.92179313711214661 1.5757413743838390996 -5.5085337243307490596 0.4150806140878786632 -0.88311183228626499719 -5.7758783199011869414;-0.26917033903266296413 -0.07902335409196255267 -2.4766993545503317442 0.73073060344902984742 -0.084117960036181646521 -0.080726819643527056303 -0.0026253531201034336837 -0.10346973201154428812 1.2562924851938910287 1.2789457612001726705;0.28314286534666449624 0.1663228115736641588 0.58278278702476027373 -0.26346755710584213173 -0.72751776908068677141 -0.011164945216807510217 -0.37083200496370372079 -0.53021501280398464218 -0.5171043557169545668 -0.61720337166534844364;0.58279122789584347775 -1.2601581500929206481 -13.239449583909518537 1.1621006646394715922 0.44826236109290257259 1.0836473204762844258 -0.29897082515746042919 -0.96546163699709763861 -0.62239307419310252456 4.2475879296178442601;0.06376004703045577704 0.4620473144817484612 2.0285500379918377512 -0.60150797265908895195 6.4181966194088889074e-05 -0.028343308672519808289 0.22033791100364413795 0.77253726431850822109 -0.65112210536172532294 -0.54121025715861570404];

% Layer 3
b3 = -0.11715619093963801356;
LW3_2 = [-1.620678925928947578 0.40701236548745933774 -0.45999550812876788264 1.1891201682351644209 0.24217874306567568898 0.86961866213884742027 -1.3941540567627799252 -3.2284964057799245829 -0.91895209179355941131 -1.966804247038343556];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.965459001831992;
y1_step1.xoffset = -0.121618166565895;

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
