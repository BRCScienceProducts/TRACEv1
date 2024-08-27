function [Y,Xf,Af] = TRACE_Preformed_P_1_Atl_3(X,~,~)
%TRACE_PREFORMED_P_1_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 21-Aug-2024 10:13:38.
% 
% [Y] = TRACE_Preformed_P_1_Atl_3(X,~,~) takes these arguments:
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
b1 = [-1.7820298053649974435;5.6987655359340880779;-1.3867712730150827394;-5.3338626288276094201;-1.5737540867766812713;0.32464291249229387804;0.062841050605324064171;-1.2050676537888178608;-0.83086550934299452376;-1.2361327584488845943;-5.6114220107155272643;1.5390104654065390744;1.958877329270688028;-0.91344278902490894101;-4.4980372802800525633];
IW1_1 = [1.2207992986684992331 0.96178572801749995236 -0.17957960629044056322 1.2525971391228285512 0.27570576399891816299 0.88740523818787298005;-2.2392668394897516926 -0.36630216571463070085 -1.242110254682889714 3.2308275446697143352 3.706935266445650079 1.0236894577719064525;0.2014257293017889805 -0.087355778717715465431 0.24802642834494020851 1.8633229924026009083 4.1660036517110103915 -4.1493705506177684939;0.65737659507458912067 -0.20492538867689405468 -0.62579187437441108344 0.71884141860947303737 -2.7830045687672755328 -6.8269918612260616086;2.2198620368832990835 -0.32570780418745137741 3.7303558249339086039 -1.5831321980881933875 0.76609283913841619285 1.5683874445516030605;0.87568498010719497682 1.0773011549021231925 1.7580199609391427007 -1.0270316302530235752 1.6254817425950538468 -2.6259285346815044448;-1.6447658875183233107 0.2127300122460420162 2.839961921044696247 1.2318862945398685316 0.46626069491149080282 -0.111340817407378459;-0.14698829243088465391 0.62411733949780623387 0.83138302429322419318 -0.61154724957320361067 -0.82511051129305956753 -1.6529882524944781164;1.1530688779834243896 -1.6957767695648797712 -2.6383522959532110796 -0.13312817292749182241 0.16200579783749463614 -2.7504096433580866865;0.94411597398509694568 -0.49567643274357686689 2.5086789653922849652 -1.373214955806396409 3.4993866488051361685 -3.522283633676166037;-0.48908415972057373011 0.49018456174802255143 0.48321238307496716669 -5.0333875113377137822 -0.77430410243635627854 -1.1858040720698268533;-0.2005380498250689536 0.70057326203722281033 0.38515722301147758966 0.52820945593543844332 0.46825219038835624863 1.4941302571410719491;0.68482139149908927145 1.0857110890787802226 0.4456657744854478076 0.63212965142207988567 0.93374930905644715207 0.68213415701837543104;-1.8364140487942282043 0.6978791138624331758 2.2828741224696957879 0.0029687224028381374412 1.5803203898491764434 -1.696088861692457872;-0.63597972992577167517 0.94344011904768743104 0.57851443506204081046 -3.2348105637180601946 0.15199700216774894579 -1.763832822696633551];

% Layer 2
b2 = [0.17190405321857848975;4.9885786135127503371;-1.4210605070317017251;-1.35541650310376105;7.3565435061304436104];
LW2_1 = [4.3234507790784899584 1.5245689967562459088 0.14940288966448622698 1.4443218581539285505 -1.740843710865260352 -3.3124934730923873971 -2.2676261541502742425 -3.2975741473984641061 -2.0106544405556383559 2.5447788226938534706 11.73079917067402711 -1.1730272801586558717 3.0587321190103295976 3.9454475752845734604 -9.9072159110244584213;1.0244007960327787821 1.3482028724775270145 0.20237989965414632243 1.2892293038962772123 -0.50171227844155863451 -1.2192277825028798866 -0.22721132992822304364 1.254569826173678404 -1.1069015422930923442 -0.53773648517987848727 -0.077093312730826382695 -1.8581664185041895099 1.3272811201127270309 -0.33263484142275029365 0.76881419831755426486;-0.17679519254344697177 -0.10528826610116254203 0.66832635234902360022 -0.12170309448435330524 0.025699706359389479227 -0.67058468982278429671 0.62317012377891389807 0.62283292606836693217 -0.23220378186948581467 1.497858248691352312 0.43173752019580463601 -0.74436276229311870711 0.10332935019424979717 -0.68262266446745700765 0.040384807716559306678;-2.6865723543056465594 -0.021835773247998985225 0.55542443165609589606 -0.04234149293663401481 0.59193556499762733925 1.1445893666345186013 1.7913688694355021092 2.9169321691205620795 0.68813564803110605084 -0.66793911285260498811 -5.851527446049527903 0.0032264221103146528691 -1.3682797573100409227 -2.3913873526041222739 4.2673567987030072146;-4.4238157197194247772 1.6312695495991493821 1.8862626388795005195 -1.3207523330146841989 0.37206283359422331891 0.4048453036141543282 9.54474104476827101 -0.20880261675435626834 0.50915091748981855524 -0.070018680055892359992 3.2979097035711197883 2.0476826704381485555 -0.36189088765847648999 -1.9483451098180146754 2.1175536716111049707];

% Layer 3
b3 = -0.039154658753706564345;
LW3_2 = [0.51470277265960495683 0.30563611014071356831 -0.24859737841532708602 0.5736408373797159399 -0.38963703733702864662];

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
