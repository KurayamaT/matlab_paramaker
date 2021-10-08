% This is the primary function named "runner.m".
% 引数：θ、θdot、γが引数。
% 出力：コマンドラインに”fixedpoint”を表示、一歩分の周期データをcsvで出力
% 必要な関数：ODE113, FSOLVE, INTERP1. 
% プログラム全体の流れ

function passivewalker_k(q,u,gamma)

tic

% motionデータ保存用のフォルダ確保
try rmdir MotionDataResults s;% フォルダが無くてもsustainさせるため
catch 
end
mkdir MotionDataResults;
% 
%     walker.M = 1000;
%     walker.m = 1.0;
%     walker.I = 0.00;
%     walker.l = .85;
%     walker.w = 0.0; 
%     walker.c = 1.0;
%     walker.r = 0.0;
%     walker.g = 1.0;

parfor k = 1:length(gamma)
    for i = 1:length(q)
        for j = 1:length(u)
            kij = append(num2str(k),'_',num2str(i),'_',num2str(j));
            disp(kij)
                
    
    %%%% Initial State %%%%%
    q1 =  q(i);
    u1 = -1*u(j);
    q2 =  2*q1; % φ　　 最初の股角度
    u2 =  -1*u1*(1-cos(q2)); % φdot 最初の股角速度
    z0 = [q1 u1 q2 u2];
    gam = gamma(k); 
    %%%
%     str_q1 = num2str(q1);
%     str_u1 = num2str(u1);
%     str_q2 = num2str(q2);
%     str_u2 = num2str(u2);
%     str_z0 = append('z0 =,  ',str_q1,',',str_u1,',',str_q2,',',str_u2);    

%%%%%%%%%%%%%%%%%%%%%%%%%
steps = 1; %number of steps to animate
% fps = 10; %Use low frames per second for low gravity

%%%% Root finding, Period one gait %%%%
options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
[zstar,~,exitflag] = fsolve(@fixedpt,z0,options,gam);
if exitflag ~= 1
    continue
else
    str_zstar = num2str(zstar);
    str_zstar = append('Fixed point = ' , str_zstar);
%   disp(str_z0)

%%% Stability, using eigenvalues of Poincare map %%%
J=partialder(@onestep,zstar,gam);
ramda = eig(J);
ramda_abs_max = max(abs(ramda));
if ramda_abs_max < 1
disp(str_zstar);
disp('Limitcycle is stable.')
disp('Motion data will be exported as a CSV file.')
disp(ramda)
%%%% Get data of leg motion. %%%
  csv_filename = filenamer(z0,gam);
  [z,~] = onestep(zstar,gam,steps);
  onestep_parameter = z;
  out = ('MotionDataResults');
  csvwrite(fullfile(out,csv_filename),onestep_parameter);   
end
end
        end
    end   
end
toc

%===================================================================
function poincaremap = fixedpt(z0,walker)
%===================================================================
poincaremap = onestep(z0,walker)-z0; 
%disp(zdiff')
%ここの右辺にある部分が、f(x)として認識され、fsolveで球解される。

%===================================================================
% データ保存用csvのファイル名を決定
function csv_filename = filenamer(z0,gam)
    str_q1 = num2str(z0(1));
    str_u1 = num2str(z0(2));
    str_gam = num2str(gam);
    csv_filename = append('onestep_parameter_',str_q1,'_',str_u1,'_',str_gam,'.csv');
%===================================================================

%===================================================================
function [z,t]=onestep(z0,gam,~)
%===================================================================

flag = 1;
if nargin<2
    error('need more inputs to onestep');
elseif nargin<3
    flag = 0; %send only last state, for root finder and jacobian
%     steps = 1;
end

q1 = z0(1);
u1 = z0(2);
q2 = z0(3);
u2 = z0(4);
z0 = [q1 u1 q2 u2];

t0 = 0; 
dt = 5; %計算間隔　might need to be changed based on time taken for one step
time_stamps = 20; %計算終了回数
t_ode = t0;
z_ode = z0;

%%% ODE solver used.

    options=odeset('abstol',1e-13,'reltol',1e-13,'events',@collision);
    %@collision関数で設定された条件（event）で計算停止を指示。
    tspan = linspace(t0,t0+dt,time_stamps);
    %y = linspace(x1,x2,n) は、x1 ～ x2 の間の等間隔の点を n 個含む行ベクトルを返します。
    [t_temp, z_temp] = ode113(@single_stance,tspan,z0,options,gam);
    % ode113⇒初期値を方程式に適用し、tspanの範囲で積分している
    % [t,y] = ode113(odefun,tspan,y0) は、tspan = [t0 tf] のときに、
    % 初期条件 y0 を使用して、微分方程式系 y′=f(t,y) を t0 から tf まで積分します。
    % 解の配列 y の各行は、列ベクトル t に返される値に対応します。
    
    % call heelstrike 
    zplus=heelstrike(t_temp(end),z_temp(end,:),gam); %　>>最終行を付加。
    
    %%%次の回転の初期値。一周しかしない場合は不要%%%
%     z0 = zplus;
%     t0 = t_temp(end);
    
    %%%%% Ignore time stamps for heelstrike and first integration point
    t_ode = [t_ode; t_temp(2:end)];
    z_ode = [z_ode; z_temp(2:end,:)];
    z = zplus(1:4);

    %%% Stability, using eigenvalues of Poincare map %%%
    
if flag==1
   z=z_ode;
   t=t_ode;
   
end

%%% 運動方程式の担当 %%%
%===================================================================
function zdot=single_stance(~,z,gam)  
%===================================================================

q1 = z(1);   u1 = z(2);                         
q2 = z(3);   u2 = z(4);                         

    M = 1000;
    m = 1.0;
    I = 0.00;
    l = .85;
    w = 0.0; 
    c = 1.0;
    r = 0.0;
    g = 1.0;
    
% 運動方程式の定義：I=w＝r=0、c=l（※l-aにてa＝0：池俣fig4.1）、で式を整理すると、池俣p20の、M11と完全に一致する！！！→20210809確認
%M11 = -1*(-2*w^2*m-2*I+2*m*l*c*cos(q2)+2*m*w*l*sin(q2)-2*m*c^2-2*m*l^2-M*l^2+2*m*l*c-2*m*r^2-M*r^2+2*m*r*c*cos(q1-q2)-2*m*r*w*sin(q1-q2)-2*M*r*l*cos(q1)-4*m*r*l*cos(q1)+2*m*r*c*cos(q1)-2*m*r*w*sin(q1)); 
M11 = -2*m*l*l*cos(q2)+2*m*l^2+M*l^2; 
M12 = -1*(w^2*m+I-m*l*c*cos(q2)-m*w*l*sin(q2)+m*c^2-m*r*c*cos(q1-q2)+m*r*w*sin(q1-q2)); 
M21 = -1*(m*w*l*sin(q2)+m*l*c*cos(q2)-m*r*w*sin(q1-q2)+m*r*c*cos(q1-q2)-m*c^2-w^2*m-I); 
M22 = -1*(w^2*m+m*c^2+I); 

% External hip torque = th.
% 以下は同じく池俣式のHとGに該当する
% RHS1 = -2*m*r*u1*u2*c*sin(q1-q2)-2*m*r*u1*u2*w*cos(q1-q2)+m*r*u1^2*w*cos(q1)+m*r*u1^2*c*sin(q1)-2*m*r*l*sin(q1)*u1^2+M*g*sin(gam)*r+2*m*g*sin(gam)*r+m*r*u2^2*w*cos(q1-q2)+m*r*u2^2*c*sin(q1-q2)+m*r*u1^2*w*cos(q1-q2)+m*r*u1^2*c*sin(q1-q2)-M*r*l*sin(q1)*u1^2+M*g*l*sin(gam-q1)+2*m*g*l*sin(gam-q1)-m*g*c*sin(gam-q1)+m*g*w*cos(gam-q1)-m*g*c*sin(gam-q1+q2)+m*g*w*cos(gam-q1+q2)-2*m*l*u1*u2*w*cos(q2)-m*l*u2^2*c*sin(q2)+2*m*l*u1*u2*c*sin(q2)+m*l*u2^2*w*cos(q2); 
RHS1 = M*g*l*sin(gam-q1)+2*m*g*l*sin(gam-q1)-m*g*c*sin(gam-q1)-m*g*c*sin(gam-q1+q2)-m*l*u2^2*c*sin(q2)+2*m*l*u1*u2*c*sin(q2); 
% RHS2 = -m*g*c*sin(gam-q1+q2)+m*g*w*cos(gam-q1+q2)-Th-m*l*u1^2*w*cos(q2)+m*l*u1^2*c*sin(q2); 
RHS2 = -m*g*c*sin(gam-q1+q2)+m*l*u1^2*c*sin(q2); 

MM = [M11 M12;                               
      M21 M22];                               

RHS = [RHS1; RHS2];                       

%　逆行列を使って、方程式noの形をθ''＝〜〜〜、とφ''＝〜〜〜にしている。「onestep.m」のode113にぶち込む。
%%% MM*X + RHS = 0; 池俣：2.1式。これより、「X＝」の式に書き直している。
X = -1*MM \ RHS;                                    

ud1 = X(1);  % θ''                                   
ud2 = X(2);  % φ''                                     

% DTE>>分解して少し書き出してみたがなんのことか良くわからん.今求めている計算過程では以下の要素は不要。axh、ayhも同様。
% DTE = -ud1*I*u2+2*ud1*m*u1*r^2+m*u1*r*u2^2*c*sin(q1-q2)+m*u1*r*u2^2*w*cos(q1-q2)-m*u2^2*l*u1*c*sin(q2)+u2*m*g*c*sin(gam-q1+q2)-u2*m*g*w*cos(gam-q1+q2)+2*ud1*I*u1+ud2*I*u2+m*u2^2*l*u1*w*cos(q2)+2*ud1*m*u1*c^2+ud2*m*u2*c^2-ud2*I*u1+ud1*m*u2*c*l*cos(q2)+ud1*m*u2*w*l*sin(q2)-2*ud1*m*u1*l*c*cos(q2)-2*ud1*m*l*u1*w*sin(q2)-m*u2*u1^2*w*l*cos(q2)+m*u2*u1^2*c*l*sin(q2)+2*ud1*m*u1*w^2+ud2*m*u2*w^2+ud2*m*u1*l*c*cos(q2)+ud1*M*l^2*u1-ud2*m*u1*w^2-ud1*m*u2*c^2-ud2*m*u1*c^2+2*ud1*m*l^2*u1-ud1*m*u2*w^2-2*ud1*m*l*u1*c+2*ud1*M*l*cos(q1)*u1*r+4*ud1*m*l*cos(q1)*u1*r-2*ud1*m*u1*r*c*cos(q1)+2*ud1*m*u1*r*w*sin(q1)-2*ud1*m*u1*r*c*cos(q1-q2)-2*m*u1^3*r*l*sin(q1)+m*u1^3*r*c*sin(q1)+m*u1^3*r*w*cos(q1)+m*u1^3*r*c*sin(q1-q2)+m*u1^3*r*w*cos(q1-q2)-2*m*u1^2*r*u2*c*sin(q1-q2)-2*m*u1^2*r*u2*w*cos(q1-q2)-M*u1^3*r*l*sin(q1)+2*u1*m*g*l*sin(gam-q1)-u1*m*g*c*sin(gam-q1)+u1*m*g*w*cos(gam-q1)+2*u1*m*g*sin(gam)*r+ud2*m*l*u1*w*sin(q2)+ud1*M*u1*r^2-u1*m*g*c*sin(gam-q1+q2)+u1*m*g*w*cos(gam-q1+q2)+u1*M*g*l*sin(gam-q1)+u1*M*g*sin(gam)*r+2*ud1*m*u1*r*w*sin(q1-q2)+ud1*m*u2*c*cos(q1-q2)*r-ud1*m*u2*w*sin(q1-q2)*r+ud2*m*u1*r*c*cos(q1-q2)-ud2*m*u1*r*w*sin(q1-q2); 
% axh = l*sin(q1)*u1^2+(-l*cos(q1)-r)*ud1; 
% ayh = -l*cos(q1)*u1^2-l*sin(q1)*ud1; 

%　ode113で使う、ud1、ud2を吐き出す。
%  zdot = [u1 ud1 u2 ud2 DTE vxh axh vyh ayh]';  %[θ’ θ'' φ’ φ''　 ]
zdot = [u1 ud1 u2 ud2]';  %[θ’ θ'' φ’ φ''　 ]

%===================================================================
function zplus=heelstrike(~,z,gam)      
%===================================================================

r1 = z(1);   v1 = z(2);                         
r2 = z(3);   v2 = z(4);                         

q1 = r1 - r2;                         
q2 = -r2;                                       

    M = 1000;
    m = 1.0;
    I = 0.00;
    l = .85;
    w = 0.0; 
    c = 1.0;
    r = 0.0;
    g = 1.0;
    

M11 = 2*m*l^2-2*m*l*c+2*m*c^2+2*m*w^2+2*m*r^2+4*m*r*l*cos(q1)-2*m*r*c*cos(q1)+2*m*w*sin(q1)*r-2*m*l*c*cos(q2)-2*m*l*w*sin(q2)-2*m*r*c*cos(q1-q2)+2*m*sin(q1-q2)*w*r+M*l^2+2*M*r*l*cos(q1)+M*r^2+2*I; 
M12 = m*l*c*cos(q2)+m*l*w*sin(q2)-m*c^2-m*w^2+m*r*c*cos(q1-q2)-m*sin(q1-q2)*w*r-I; 

M21 = -m*l*c*cos(q2)-m*l*w*sin(q2)+m*c^2+m*w^2-m*r*c*cos(q1-q2)+m*sin(q1-q2)*w*r+I; 
M22 = -m*w^2-m*c^2-I; 

RHS1 = m*l*v2*c-m*c^2*v2+M*v1*r^2-2*m*c*l*v1+M*v1*l^2*cos(r2)+2*m*l^2*v1*cos(r2)+2*I*v1-I*v2+2*m*v1*r^2-2*m*l*v1*c*cos(r2)+2*m*c^2*v1+2*m*w^2*v1-m*w^2*v2+2*m*r*v1*w*sin(r1)+M*r*v1*l*cos(r1)+M*l*cos(-r1+r2)*v1*r-2*m*r*v1*c*cos(r1)+2*m*l*cos(-r1+r2)*v1*r+2*m*r*v1*l*cos(r1)-2*m*r*v1*c*cos(-r1+r2)-2*m*r*v1*w*sin(-r1+r2)+m*r*v2*c*cos(-r1+r2)+m*r*v2*w*sin(-r1+r2); 
RHS2 = m*r*v1*w*sin(r1)-m*r*v1*c*cos(r1)+I*v1-I*v2+m*w^2*v1-m*c*l*v1+m*c^2*v1; 

MM = [M11 M12;                               
     M21 M22];    

RHS = [RHS1; RHS2];                      

X = MM \ RHS;                                    

u1 = X(1);                                       
u2 = X(2);                                      

zplus = [q1 u1 q2 u2 ];      


%===================================================================
function [gstop, isterminal,direction]=collision(~,z,~)
%===================================================================

q1 = z(1); q2 = z(3); 

gstop = -q2 + 2*q1;%毎回の計算におけるgstopの値を検査（gstop=0を検出。「イベントを記述する数式です。value(i) がゼロに等しくなると、イベントが発生します。」）
if (q2>-0.05) %allow legs to pass through for small hip angles (taken care in real walker using stepping stones)
    isterminal = 0;%股関節の角度が小さい時はカウントしない。
else%そうでない場合は
    isterminal=1; %ode should terminate is conveyed by 1, if you put 0 it goes till the final time u specify
    %　isterminal＝＝１で計算ストップ
end
direction=-1; % The t_final can be approached by any direction is indicated by the direction
% direction=-1 ではイベント関数(in this case, gstop.)が減少している零点のみが検出されます。


%===================================================================
function J=partialder(FUN,z,walker)
%===================================================================
pert=1e-5;%pertuabation 
n = length(z);
J = zeros(n,n);

%%%% Using forward difference, accuracy linear %%%
% y0=feval(FUN,z,walker); 
% for i=1:n
%     ztemp=z;
%     ztemp(i)=ztemp(i)+pert; 
%     J(:,i)=(feval(FUN,ztemp,walker)-y0) ;
% end
% J=(J/pert);

%%% Using central difference, accuracy quadratic %%%
for i=1:n
    ztemp1=z; ztemp2=z;
    ztemp1(i)=ztemp1(i)+pert; 
    ztemp2(i)=ztemp2(i)-pert; 
    J(:,i)=(feval(FUN,ztemp1,walker)-feval(FUN,ztemp2,walker)) ;
end
J=J/(2*pert);

