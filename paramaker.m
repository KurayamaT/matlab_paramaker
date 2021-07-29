% make initial parameters [q1(sheta), u1(sheta-dot), g(gamma)] with the
% range of defined value.

function [q1, u1, gam] = paramaker(q_start, q_end, q_points, u_start, u_end, u_points, gam_start, gam_end, gam_points);

delete q1.mat;
delete u1.mat;
delete gam.mat;

q1 = (q_start:(q_end-q_start)/q_points:q_end)';
u1 = (u_start:(u_end-u_start)/u_points:u_end)';
gam = (gam_start:(gam_end-gam_start)/gam_points:gam_end)';

save q1.mat q1;
save u1.mat u1;
save gam.mat gam;



