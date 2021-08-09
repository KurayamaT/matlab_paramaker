% This is the primary function named "runner.m".
% 引数：θ、θdot、γが引数。
% 出力：コマンドラインに”fixedpoint”を表示、一歩分の周期データをcsvで出力
% 必要な関数：ODE113, FSOLVE, INTERP1. 

function runner(q,u,gam)

delete *.csv 

for i = 1:length(q)
    for j = 1:length(u)
        for k = 1:length(gam)

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
[zstar,fval,exitflag] = fsolve(fx,z0,options,walker);
if exitflag ~= 1
       % error('Root finder not converged, change guess or change system parameters')

    continue

else
    %disp(str_z0);
    str_zstar = num2str(zstar);
    str_zstar = append('Fixed point = ' , str_zstar);
    disp(str_zstar);
    %disp (str_gam);
    disp('Motion data were exported as a CSV file.')
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


              
