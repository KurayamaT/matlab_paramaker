% make initial parameters [q1(sheta), u1(sheta-dot), g(gamma)] with the
% range of defined value.



function [q1, u1, g] = paramaker(q_start, q_end, q_points, u_start, u_end, u_points, gam_start, gam_end, gam_points)

q1 = q_start:(q_end-q_start)/q_points:q_end
u1 = u_start:(u_end-u_start)/u_points:u_end
q1 = q_start:(q_end-q_start)/q_points:q_end



