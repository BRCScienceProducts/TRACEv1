function [Y,Xf,Af] = TRACE_Preformed_Si_1_Other_3(X,~,~)
%TRACE_PREFORMED_SI_1_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 21-Aug-2024 10:13:39.
% 
% [Y] = TRACE_Preformed_Si_1_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999961923064171;-0.999961923064171;-77.5;0;30.1649227142334;-2.9811577796936];
x1_step1.gain = [1.00003807838574;1.00003807838574;0.0119760479041916;0.000363636363636364;0.22152510192902;0.0532708976942295];
x1_step1.ymin = -1;

% Layer 1
b1 = [28.835181575575017376;-0.56245351457895964575;3.577714649926821977;14.612743895576262076;0.92384639725801898091;1.3563550912139423765;-0.97478289507371829181;3.0895604202967534846;0.60409462805588820533;-27.234359064038372367;-1.40932919953278879;-1.6736646478170473795;25.875420086558015953;-0.54365964482705431138;1.0316164416913513602];
IW1_1 = [-14.758466859917138692 8.73462385544909381 -6.1258860851462930697 10.344255947425503805 1.3825499066460118502 4.4688873627883110018;0.19416259029185065077 -0.00087938699385753903159 0.2204414096226208486 0.087057779770638374361 7.6422142436621109951 0.056746851285764239592;-0.33140748612786941063 -0.71578918253667644933 0.19715514046452251029 -0.11801835231237574331 1.4644991188094704793 2.6060179415658741853;0.015756157280582359231 -0.12705057821256743322 1.1252064678555933064 10.830474221972085758 -8.4821430259326113799 3.1706687661853676374;0.29653521820252670693 0.028140363020322359949 1.6795151199486972882 0.25465446546375791748 -4.2015554719552241636 2.1669347846180824213;0.020047261876451220874 -0.070608199622377756932 0.4159446614399297526 0.0083837669629050273756 -1.4116303964051180042 0.82477171584020847472;0.30424507710561854434 -0.080709941826605721693 0.4286606761326712256 0.24771512718318858259 -13.493410360326731023 -0.37591387990650981799;-4.7032618427720169763 0.99401462095952053399 -3.4450946127743478087 -2.0748036841590788981 -21.924184432318714499 4.3967241187237844358;-12.553716877060264068 8.0800812391403304957 -19.818684666266104699 -0.35759262683871056998 -6.9798533314826549656 3.050525735140584338;4.0003793923142065836 -2.8165598297682525875 -19.209173155822014678 -9.4285797975297178652 0.41042122449896711878 5.7179278425990442969;3.8836338879663534485 1.4804332693232793439 3.3817519386142702587 -0.015903718499465990116 1.5872195243815943044 -4.8473573104938720846;-0.10952299914181556906 0.14231125185164533731 -0.74500359802806748188 -0.021317261686647284713 4.5104606330101448819 -1.3986697824628067099;-0.050777994405563527114 -0.23029553052905768484 2.2785349501518847859 21.363171256450392832 -12.242076011955273174 4.0256753252035704094;0.24914952592808184417 -0.14345613372484741777 1.3777492279693475119 0.45907929546380737973 -2.1291847147732019252 -1.7270332314829939246;-1.0783333742410994294 -0.17165976833436055737 0.9818849190627592094 1.4302817987007130363 -18.220288442762608838 6.6775716879968793904];

% Layer 2
b2 = [-15.879556754004195085;-20.532750843082908432;-0.62189889575644607955;-4.5088307670420473272;-16.421077432528939255];
LW2_1 = [11.776632211737076972 -2.23112383694729699 1.7129343827317025095 33.098391655796717714 3.1770975253316215969 -40.642218725360422127 -1.274656127366859959 1.2372798090726122577 -0.10203473920929410912 1.2385100865495559663 -0.45917852369662820466 -13.361177268259462281 -13.257054943139280567 2.3223561533967336246 -0.15934519675909117331;47.86937517107605089 -1.7930996727325214657 17.913057580571983607 17.114251038163583729 33.432851081660892589 -13.466293673343765036 -43.518160078296070026 13.857145377484673787 -11.673974373908100333 23.660119883599044499 -3.6298148888659076405 7.3026198251489837787 -22.530901976447861301 70.935739070554618024 -21.289319544732649092;-1.7182105791771855596 -2.5813241410595049175 -21.95163784563555609 7.8247558397062846325 -3.1842850531472666198 50.377796462060892679 -4.4208892898495895096 3.1355667001515179493 -2.4364642356069530216 13.988850961757099967 -3.0559885176407020069 0.12632712293746278798 -3.4968430118413951213 -6.7506912073103828575 -0.049533896458062831591;-6.872660079568767344 0.3567929541219447831 -17.422495903061992806 3.5146304486614408091 -1.5543386407212702593 6.3382554417624525911 0.19447061426380443661 1.4531140840571223904 20.020028859953111322 5.7796320797993754326 -1.5601209072462485228 -0.17417774155293075378 1.9258005486098808134 -0.29608240043778716322 0.46537035051804570562;-10.832694905634479454 3.0222421074145753295 -0.4769801908088436293 -5.8137259241368415275 -1.2714086476708046369 20.226526804570202245 0.9301837791255420651 -0.33759723165682414825 0.040030853018899657803 -22.56980847554538272 0.14645552165812722456 6.6532128868600857174 4.3927962945474261502 -1.3779668832298235248 -0.38468911283271450818];

% Layer 3
b3 = -0.36531412907827714998;
LW3_2 = [-0.20746943699693237706 0.014752725412768796093 -0.033568371547394330878 -0.059851636100220775849 -0.78972058857057736425];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0184350799221364;
y1_step1.xoffset = -0.488903820514679;

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
