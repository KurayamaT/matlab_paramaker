% make initial parameters [q1(theta), u1(theta-dot), g(gamma)] with the range of defined value.
% θの開始と終了値、point数、θ'の開始と終了値、point数、γの開始と終了値、point数を入力。

function [q,u,gamma] = paramaker_k(q_start, q_end, q_points, u_start, u_end, u_points, gam_start, gam_end, gam_points);
%%%実行するときはここで定義しているfunction名はあまり関係なかったりする。ファイル名そのものに依存。

format long;

% delete q1.mat;
% delete u1.mat;
% delete gam.mat;

q  = q_start:(q_end-q_start)/q_points:q_end;
u  = u_start:(u_end-u_start)/u_points:u_end;
gamma = gam_start:(gam_end-gam_start)/gam_points:gam_end;

passivewalker_k(q,u,gamma);
% save q1.mat q;
% save u1.mat u;
% save gam.mat gam;

%%% ↓お試し用：入力数値コマンド
% paramaker(0.01745,0.5236,10,0.256,2.56,10,0.001,1.57,10)
%%%　2021/10/8：入力数値
% paramaker(0.2618,0.5236,30,0.256,2.56,30,0.0175,0.5236,100)




