
%===================================================================
function [gstop, isterminal,direction]=collision(t,z,walker)
%===================================================================

M = walker.M;  m = walker.m; I = walker.I;   
l = walker.l;  c = walker.c; w = walker.w;   
r = walker.r;  g = walker.g; gam = walker.gam; 

q1 = z(1); q2 = z(3); 

gstop = -q2 + 2*q1;%毎回の計算におけるgstopの値を検査（gstop=0を検出。「イベントを記述する数式です。value(i) がゼロに等しくなると、イベントが発生します。」）
if (q2>-0.05) %allow legs to pass through for small hip angles (taken care in real walker using stepping stones)
    isterminal = 0;%股関節の角度が小さい時はカウントしない。
else%そうでない場合は
    isterminal=1; %ode should terminate is conveyed by 1, if you put 0 it goes till the final time u specify
    %　isterminal＝＝１で計算ストップ
end
direction=-1; % The t_final can be approached by any direction is indicated by the direction
% direction=-1 ではイベント関数(in this case, gstop.)が減少している零点のみが検出されます。
