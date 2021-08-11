% make initial parameters [q1(theta), u1(theta-dot), g(gamma)] with the range of defined value.
% θの開始と終了値、point数、θ'の開始と終了値、point数、γの開始と終了値、point数を入力。

function [q,u,gamma] = paramaker(q_start, q_end, q_points, u_start, u_end, u_points, gam_start, gam_end, gam_points);
%%%実行するときはここで定義しているfunction名はあまり関係なかったりする。ファイル名そのものに依存。

format long;

% delete q1.mat;
% delete u1.mat;
% delete gam.mat;

q  = q_start:(q_end-q_start)/q_points:q_end;
u  = u_start:(u_end-u_start)/u_points:u_end;
gamma = gam_start:(gam_end-gam_start)/gam_points:gam_end;

passivewalker_k(q,u,gamma)
% save q1.mat q;
% save u1.mat u;
% save gam.mat gam;

%%% ↓お試し用：入力数値コマンド
% paramaker(0.017453293,0.13962634,3,0.017453293,0.13962634,3,5.81776E-05,0.000523599,3)
% marged_runner_and_paramaker(0.017453293,0.13962634,3,0.017453293,0.13962634,3,5.81776E-05,0.000523599,3)




