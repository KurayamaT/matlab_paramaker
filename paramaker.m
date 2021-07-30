% make initial parameters [q1(sheta), u1(sheta-dot), g(gamma)] with the range of defined value.

function [q,u,gam] = paramaker(q_start, q_end, q_points, u_start, u_end, u_points, gam_start, gam_end, gam_points);
%%%実行するときはここで定義しているfunction名はあまり関係なかったりする。ファイル名そのものに依存。

format long;

% delete q1.mat;
% delete u1.mat;
% delete gam.mat;

q  = (q_start:(q_end-q_start)/q_points:q_end);
u  = (u_start:(u_end-u_start)/u_points:u_end);
gam = (gam_start:(gam_end-gam_start)/gam_points:gam_end);

runner(q,u,gam)
% save q1.mat q;
% save u1.mat u;
% save gam.mat gam;


%%% ↓
% paramaker(0.017453293,0.13962634,10,0.017453293,0.13962634,10,5.81776E-05,0.000523599,10)
% marged_runner_and_paramaker(0.017453293,0.13962634,3,0.017453293,0.13962634,3,5.81776E-05,0.000523599,3)




