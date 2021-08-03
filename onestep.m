

%===================================================================
function [z,t]=onestep(z0,walker,steps)
%===================================================================

M = walker.M;  m = walker.m; I = walker.I;   
l = walker.l;  c = walker.c; w = walker.w;   
r = walker.r;  g = walker.g; gam = walker.gam;

%%% written by TK %%%
    % for save onstep data.
    str_q1 = num2str(z0(1));
    str_u1 = num2str(z0(2));
    str_gam = num2str(gam);
    fname2 = append('onestep_parameter_',str_q1,'_',str_u1,'_',str_gam,'.csv');
%%%%%%%%%%%%%%%%%%%%%%%

flag = 1;
if nargin<2
    error('need more inputs to onestep');
elseif nargin<3
    flag = 0; %send only last state, for root finder and jacobian
    steps = 1;
end

q1 = z0(1);
u1 = z0(2);
q2 = z0(3);
u2 = z0(4);

    %%%% Derived variables %%%%
    TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1; 
    xp1 = 0;
    xh = -l*sin(q1) - r*q1 + xp1;
    vxh = (-l*cos(q1)-r)*u1; 
    yh =  l*cos(q1) + r;
    vyh = -l*sin(q1)*u1; 
    
z0 = [q1 u1 q2 u2 TE xh vxh yh vyh];

t0 = 0; 
dt = 5; %might need to be changed based on time taken for one step
time_stamps = 20;
t_ode = t0;
z_ode = z0;

for i=1:steps
    options=odeset('abstol',1e-13,'reltol',1e-13,'events',@collision);
    tspan = linspace(t0,t0+dt,time_stamps);
    [t_temp, z_temp] = ode113(@single_stance,tspan,z0,options,walker);
    
    zplus=heelstrike(t_temp(end),z_temp(end,:),walker); 
    
    z0 = zplus;
    t0 = t_temp(end);
    
    %%%%% Ignore time stamps for heelstrike and first integration point
    t_ode = [t_ode; t_temp(2:end)];
    z_ode = [z_ode; z_temp(2:end,:)];
    onestep_parameter = z_temp(2:end,:);
    
end


z = zplus(1:4);

if flag==1
   z=z_ode;
   t=t_ode;
   
   %%% written by TK %%%
    csvwrite(fname2,onestep_parameter);  
   %%%%%%%%%%%%%%%%%%%%%
   
end
