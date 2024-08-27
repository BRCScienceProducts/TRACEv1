function [Y,Xf,Af] = TRACE_Preformed_Si_1_Atl_4(X,~,~)
%TRACE_PREFORMED_SI_1_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 21-Aug-2024 10:13:39.
% 
% [Y] = TRACE_Preformed_Si_1_Atl_4(X,~,~) takes these arguments:
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
b1 = [0.19579019037680589488;0.92173202347937899237;3.6356578205187384256;5.8779722783281895815;-0.57554032501642005126];
IW1_1 = [-0.046236509916038767698 0.030328418423579233537 0.092135664511798942566 -0.55588364265420764365 0.117968347276050739 -0.79726635366980058706;-0.030146707005643404165 -0.022802572493041472351 0.32994481482528775773 0.98155066221164377627 0.38404440539050405112 0.10306125812321917379;0.96908717209313399721 -0.61747111178849867219 2.3694156909316199666 0.61626212325755136412 3.6035538440698968721 3.9382591743819448382;-0.14917779202748504552 -0.098798518682828129212 1.0844874948950522509 5.7312796511332528837 1.1756731240226014901 0.42937137070933328875;0.2759909884317242601 -0.29135973582746760302 1.3579259688342735224 0.092953477071083659422 1.1994205039028600357 -1.5315456048323492499];

% Layer 2
b2 = [-15.312245755252893886;16.135221159302925997;-2.4591411730286982795;-6.1761346898546616302;-1.2858779641191442344;0.16007587028982966926;-5.206062838842143492;2.5444702307857256507;-7.091010418149087613;-6.4683659554057610919;-0.16984109524946844338;3.3382918655770832395;4.295756812827415061;2.4718475377087849587;-14.342525472335086789];
LW2_1 = [1.718408863089723182 1.784526970358638609 -0.52441808152192370684 11.444839371545048934 0.1874494375540191049;-15.162093809051278015 -3.2772589290020306407 4.4982867682159159983 -4.3902680537591702858 -1.4662762453767355098;11.896569533030580956 7.6775778824838996073 2.0320596104069834631 -7.4848865877351515863 -5.0820428908366013232;-2.5178988220133433629 1.6941719010136264689 -5.1113494785942048537 -3.320713168938664861 0.6201675373938315694;4.8448830176738866271 14.027315099163432777 1.0052407192959693649 -3.5390422022113798839 -3.0514307525435397928;-1.9230336336012858744 2.9308319908353812622 -0.33160511386679147883 0.36340308189967485131 0.44580034263805723427;3.6122586382256884896 -12.009936694969441717 2.4987918919865945711 -1.6323625956233860368 4.1573483828243951166;-4.0694363430215876321 -2.7994483907160305769 -1.7141903810472391978 -2.214845219557105338 4.1376975607292214576;-4.9220835726972920199 0.94648448728586320744 -5.2313339802648579635 6.8140449989815436993 9.1949949430215500712;-19.844052590587750728 21.620221716374736332 0.14210145938524135212 -9.2144035193685116525 -4.2718381722380396326;-0.80710416960358521976 -1.7245507637265782286 3.6068691408224751171 -2.3967241288117531184 2.9765907574854053586;-5.20675320223426219 -5.2844141123227306167 -0.28626218662771713008 1.6350449465919276726 0.48546278895129818531;-8.6217057260762395288 -3.2449785582073116252 0.19384196616350660158 -0.99895601869989647525 6.5707026693424470665;-12.252061363254176385 -7.9262304363231095294 -1.9720275520867915908 7.6653697982100004182 5.2343951302907827383;12.832318754532902716 10.702796948817518086 -0.24517857565030193312 -6.4019169790641674567 2.4039232602497428992];

% Layer 3
b3 = 5.3712552116465204222;
LW3_2 = [9.1932679069610525602 0.16857273281236956342 -5.2744170062008182143 -3.0820874741405179797 -0.41515199157719612222 -0.24418784088186865078 -0.047453741585597841723 -0.20400895583449971116 -0.12451505758058063755 0.44590556173911927162 -0.067881163075109704486 -0.23404230057426494982 -0.089787157058186559122 -4.9787730275863397722 -0.40383698149986291259];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0510398195720568;
y1_step1.xoffset = -0.473623931407928;

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
