
%===================================================================
function poincaremap = fixedpt(z0,walker)
%===================================================================
poincaremap = onestep(z0,walker)-z0; 
%disp(zdiff')