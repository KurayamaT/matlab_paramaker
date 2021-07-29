% make initial parameters [q1(sheta), u1(sheta-dot), g(gamma)] with the range of defined value.

function marged_runner_and_paramaker(q_start, q_end, q_points, u_start, u_end, u_points, gam_start, gam_end, gam_points);
%%%実行するときはここで定義しているfunction名はあまり関係なかったりする。ファイル名そのものに依存。

format long

delete q1.mat;
delete u1.mat;
delete gam.mat;

q  = (q_start:(q_end-q_start)/q_points:q_end);
u  = (u_start:(u_end-u_start)/u_points:u_end);
gam = (gam_start:(gam_end-gam_start)/gam_points:gam_end);

save q1.mat q;
save u1.mat u;
save gam.mat gam;

%run(runner);

% paramaker(0.017453293,0.13962634,10,0.017453293,0.13962634,10,5.81776E-05,0.000523599,10)
% marged_runner_and_paramaker(0.017453293,0.13962634,10,0.017453293,0.13962634,10,5.81776E-05,0.000523599,10)

% This is the primary function named "runner.m".
% 引数：θ、θdot、γが引数。
% 出力：コマンドラインに”fixedpoint”を表示、一歩分の周期データをcsvで出力
% 必要な関数：ODE113, FSOLVE, INTERP1. 

%function runner

delete *.csv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% parameters for tests. %%%

% theta_rad = 0.261799388;
% theta_dot = 0.34906585;
% gam    =.104719755;
% %%%%%%%%%%%%%%%%%%%%%
% 
% theta_rad = [0.017453293	0.034906585	0.052359878	0.06981317	0.087266463	0.104719755	0.122173048	0.13962634	0.157079633	0.174532925	0.191986218	0.20943951	0.226892803	0.244346095	0.261799388	0.27925268	0.296705973	0.314159265	0.331612558	0.34906585	0.366519143	0.383972435	0.401425728	0.41887902	0.436332313	0.453785606	0.471238898	0.488692191	0.506145483	0.523598776];
% theta_dot = [0.017453293	0.034906585	0.052359878	0.06981317	0.087266463	0.104719755	0.122173048	0.13962634	0.157079633	0.174532925	0.191986218	0.20943951	0.226892803	0.244346095	0.261799388	0.27925268	0.296705973	0.314159265	0.331612558	0.34906585	0.366519143	0.383972435	0.401425728	0.41887902	0.436332313	0.453785606	0.471238898	0.488692191	0.506145483	0.523598776];
% gam       = [5.81776E-05	0.000116355	0.000174533	0.000232711	0.000290888	0.000349066	0.000407243	0.000465421	0.000523599	0.000581776	0.000639954	0.000698132	0.000756309	0.000814487	0.000872665	0.000930842	0.00098902	0.001047198	0.001105375	0.001163553	0.00122173	0.001279908	0.001338086	0.001396263	0.001454441	0.001512619	0.001570796	0.001628974	0.001687152	0.001745329];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load q1.mat;  % q1
load u1.mat;  % u1
load gam.mat; % gam


for i = 1:length(q)
    for j = 1:length(u)
        for k = 1:length(gam)

    %% Garcia's simplest walker with roots for Validation
    %%% Dimensions %%
    %% c = COM on the leg from hip, w = COM fore-aft offset, r = radius of feet
    %% M = hip mass, m = leg mass, I = leg inertia, l = leg length
    %%%%% To get results close to Garcia's walker increase M %%%%%%
    walker.M = 1000;
    walker.m = 1.0;
    walker.I = 0.00;
    walker.l = .85;
    walker.w = 0.0; 
    walker.c = 1.0;
    walker.r = 0.0;
    walker.g = 1.0;
    walker.gam = gam(k);
    str_gam = num2str(walker.gam);
    str_gam = append('gam, ',str_gam);
    
    %%%% Initial State %%%%%
    q1 =  q(i);
    u1 = -1*u(j);
    q2 =  2*q1; % φ　　 最初の股角度
    u2 =  -1*u1*(1-cos(q2)); % φdot 最初の股角速度
    z0 = [q1 u1 q2 u2];
    
    %%%
    str_q1 = num2str(q1);
    str_u1 = num2str(u1);
    str_q2 = num2str(q2);
    str_u2 = num2str(u2);
    str_z0 = append('z0 =,  ',str_q1,',',str_u1,',',str_q2,',',str_u2);
    

%%%%%%%%%%%%%%%%%%%%%%%%%
steps = 1; %number of steps to animate
fps = 10; %Use low frames per second for low gravity

%%%% Root finding, Period one gait %%%%
options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
[zstar,fval,exitflag] = fsolve(@fixedpt,z0,options,walker);
if exitflag ~= 1
       % error('Root finder not converged, change guess or change system parameters')

    continue

else
    disp(str_z0);
    str_zstar = num2str(zstar);
    str_zstar = append('Fixed point =, ' , str_zstar);
    disp(str_zstar);
    disp (str_gam);
    disp('A Motion data was exported as a CSV file.')
    %disp(theta_dot);
end

%%% Stability, using eigenvalues of Poincare map %%%
J=partialder(@onestep,zstar,walker);
% disp('EigenValues for linearized map are');
% eig(J)
 
%%%% Get data for all the steps %%%
[z,t] = onestep(zstar,walker,steps);


        end
    end
end


              
