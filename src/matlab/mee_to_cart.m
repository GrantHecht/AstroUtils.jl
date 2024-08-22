clear; close all; clc

syms p f g h k L real;
syms rx ry rz vx vy vz real;
syms mu real positive;

kk = k*k;
hh = h*h;
tkh = 2.0*k*h;
s2 = 1.0 + hh + kk;
cL = cos(L);
sL = sin(L);
w = 1.0 + f*cL + g*sL;
r = p / w;
smp = sqrt(mu/p);
fhat = [1.0 - kk + hh; tkh; -2.0*k];
ghat = [tkh; 1.0 + kk - hh; 2.0*h];
fhat = fhat / s2;
ghat = ghat / s2;
x = r*cL;
y = r*sL;
xdot = -smp*(g + sL);
ydot = smp*(f + cL);

r = x*fhat + y*ghat;
v = xdot*fhat + ydot*ghat;

mee = [p;f;g;h;k;L];
x = [r; v];

J = sym('J', [6 6]);
for i = 1:6
    for j = 1:6
        J(i,j) = simplify(diff(x(i),mee(j)));
    end
end

ccode(J, 'file', './mee_to_cart_partials.c');