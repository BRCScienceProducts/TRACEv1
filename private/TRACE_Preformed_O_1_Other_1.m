function [Y,Xf,Af] = TRACE_Preformed_O_1_Other_1(X,~,~)
%TRACE_PREFORMED_O_1_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 21-Aug-2024 10:13:37.
% 
% [Y] = TRACE_Preformed_O_1_Other_1(X,~,~) takes these arguments:
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
b1 = [-1.4245263646701140647;-0.53432884985264339495;0.90520960311672482579;-14.015653492163215432;8.3095987145943759344;1.2266128799764470791;-7.1727663105748185046;0.85714168819681524347;0.30682056362596388244;55.685674069916693441;0.24262851397347381943;-0.52432865897216873741;2.8523257044826015871;-1.126491887787474111;2.9006905963773879087;-1.2072189990611963406;-0.48150824196215558093;0.50735630028656530932;3.4118244651339511009;-5.4957039988129494645];
IW1_1 = [0.11607718411554744653 -0.023558106189151578613 -0.13662426377284286882 -1.465702406731530294 -1.6321013480066475232 -0.35076295887013320129;-0.040408657012607450409 0.047452014551024643274 2.4570867658481638784 -5.4633500734033617618 -8.2173383131363522125 2.2093565384292706;-2.5974383820677839907 -1.6556213003933000216 0.59003035930108416185 1.6737957706582256012 2.1509663711533888453 -3.1708875529669016125;2.9152222626678634398 -1.2287780499088083452 -11.643028704622176761 -0.86881223119817496681 16.070618963361276599 -4.299653388435368484;-3.8509009875015913771 1.7327435901414911523 -1.9065559147295279008 1.5847486131521950625 -8.018146099355107026 4.2628867849770673004;-0.093127855846950574614 0.019068509769589424285 -0.027392088475772444212 1.2917552582751528067 2.2312288184111932487 0.16202773615801160934;10.723121335335179083 -9.2196943448717174618 5.7928519058049054991 6.3716660861741933886 -3.0082950635464902334 0.78961829540423633667;-2.6096725433632110835 -1.6600455125666746525 0.58110509985401292976 1.6160186970751118096 2.3006084830955568776 -3.2051848143576084915;-0.027980481844617190179 0.029244529263720755313 2.0554368193500183182 -4.3187924430752575589 -7.8193772500733427222 2.8334145072990284397;0.08829868584327127401 -0.17824234807395139546 -0.13835813964989479574 53.733864468005400283 -0.59904485899589299258 0.57687389621372620052;0.051613467302812875148 0.01074057428049638406 0.6980678504613559987 -1.0088806940232006859 0.49985153280539085241 -0.25468395211044375648;-0.051292984366148222652 -0.012105549763953312239 -0.70166281170848188165 1.068092667130015716 -1.2265732207598811776 0.030141575926926449103;0.81745078136241278166 -0.039531227003071656156 -0.61010161781537375347 1.8913507775313525006 -2.6560167450538481049 1.8288679986820706169;0.066553066490943654787 -0.01062126095458450388 0.26773503244768154241 -1.1877373520368819104 -2.7291397695910473153 0.039405490419608238395;0.7942381221190264462 0.0252172614236992991 -0.5779401076805821269 1.8349035402976574538 -2.709489367180369257 1.9367934071409822483;0.0510071791877703129 -0.0047083430702746292482 0.44795096479463403227 -1.2636794962573449919 -3.0830320615904409642 0.20221176126021883634;-0.13025191588209672489 0.70502934329075694908 3.0196920612316837662 0.73355937462991205322 -2.7015954582053840483 3.9316853757093861432;0.11772613108352279754 -0.79485629940030200213 -3.111301757428925896 -0.72672379429722522293 2.7604084050639006165 -4.0420806924082386402;0.87814462876177235007 0.73725247870514731829 -0.48882164761596641478 1.4790379428226134717 1.7500904494159530422 2.4389163026201239326;-1.885987550667448609 1.1272813035472080134 -3.3046740179373523816 -1.9997905414151788683 1.3963796088725926037 2.051151822754980536];

% Layer 2
b2 = 2.9602036869649421824;
LW2_1 = [12.104024069594647273 0.61099053146115034352 2.1972845218350434138 0.034130773977363908223 -0.14032767817084809692 31.300915491570709293 0.026126087671249634964 -2.2022075257580970842 -0.9050371893567540349 -2.3044785004569741105 -4.8685892948221383136 -4.2513285613060114443 -2.1350178648322102326 36.135633640885963302 2.3120899887623687086 -15.85848660825604739 2.8685581982191661332 2.6945698018752293201 -0.22858051035259815964 -0.10507217287101645875];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00612908122568115;
y1_step1.xoffset = 89.6956371703064;

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
    a2 = repmat(b2,1,Q) + LW2_1*a1;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a2,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(2,0);

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
