%% calculate the value of K
function K=computeK(E,A,v0,mstar)

K=mstar*v0^2*E*A/(E*A-mstar*v0^2);

end
