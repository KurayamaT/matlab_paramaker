% make initial parameters [q1(theta), u1(theta-dot), g(gamma)] with the range of defined value.
% θの開始と終了値、point数、θ'の開始と終了値、point数、γの開始と終了値、point数を入力。

function [q,u,gamma,tau] = paramaker_k_tau(q_start, q_end, q_points, u_start, u_end, u_points, gam_start, gam_end, gam_points, tau_start, tau_end, tau_point)
%%%実行するときはここで定義しているfunction名はあまり関係なかったりする。ファイル名そのものに依存。

format long;

% delete q1.mat;
% delete u1.mat;
% delete gam.mat;

q  = q_start:(q_end-q_start)/q_points:q_end;
u  = u_start:(u_end-u_start)/u_points:u_end;
gamma = gam_start:(gam_end-gam_start)/gam_points:gam_end;
tau = tau_start:(tau_end-tau_start)/tau_point:tau_end;
passivewalker_k_tau(q,u,gamma,tau);
%  save q1.mat q;
%  save u1.mat u;
%  save gam.mat gamma;
%  save th.mat tau;

%%% ↓お試し用：入力数値コマンド
%paramaker_k_tau(0.01745, 0.5236, 3, 0.256, 2.56, 3, 0.001, 1.57, 3, 0.001, 0.01, 5)




