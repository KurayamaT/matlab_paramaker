% make initial parameters [q1(theta), u1(theta-dot), g(gamma)] with the range of defined value.
% θの開始と終了値、point数、θ'の開始と終了値、point数、γの開始と終了値、point数を入力。

function [q,u,gamma,tau] = paramaker_k_tau(q_start, q_end, q_points, u_start, u_end, u_points, gam_start, gam_end, gam_points, tau_start, tau_end, tau_point)
%%%実行するときはここで定義しているfunction名はあまり関係なかったりする。ファイル名そのものに依存。

format long;

% delete q1.mat;
% delete u1.mat;
% delete gam.mat;
if q_start == q_end
q = q_start;
else
q  = q_start:(q_end-q_start)/q_points:q_end;
end

if u_start == u_end
    u = u_start;
else
u  = u_start:(u_end-u_start)/u_points:u_end;
end

if gam_start == gam_end
    gamma = gam_start;
else
gamma = gam_start:(gam_end-gam_start)/gam_points:gam_end;
end

if tau_start == tau_end
    tau = tau_start;
else
tau = tau_start:(tau_end-tau_start)/tau_point:tau_end;
end

total_calctimes = (q_points+1) * (u_points+1) * (gam_points+1) * (tau_point+1);
disp(total_calctimes); 
passivewalker_k_tau(q,u,gamma,tau);
%  save q1.mat q;
%  save u1.mat u;
%  save gam.mat gamma;
%  save th.mat tau;

%%% ↓お試し用：入力数値コマンド
%paramaker_k_tau(0.01745, 0.5236, 3, 0.256, 2.56, 3, 0.001, 1.57, 3, 0.001, 0.01, 5)

%% Pranav setting
% paramaker_k_tau(0.3604,0.3605,1, -0.3736, -0.3737, 1, 0.115, 0.116, 1, 0.001, 0.1, 10)

%20211019PM1
%paramaker_k_tau(0.3,0.4,10, -0.4, -0.3, 10, 0.087, 0.122, 10, 0.01, 0.1, 10)
%20211019PM2
%paramaker_k_tau(0.36,0.37,1, -0.37, -0.38, 1, 0.115, 0.116, 1, 0.01, 0.1, 10)
%20211019last
%paramaker_k_tau(0.34,0.47,13, -3.0, -2.0, 10, 0.1, 0.3, 20, 0.01, 0.3, 30)





