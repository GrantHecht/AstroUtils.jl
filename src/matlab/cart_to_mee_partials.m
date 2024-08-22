clear; close all; clc

syms p f g h k L real;
syms rx ry rz vx vy vz real;
syms mu real positive;

r = [rx;ry;rz];
v = [vx;vy;vz];
x = [r; v];

rdv = dot(r,v);
rmag = norm(r);
rhat = r / rmag;
hvec = cross(r,v);
hmag = norm(hvec);
hhat = hvec / hmag;
vhat = (rmag*v - rdv*rhat) / hmag;
p = hmag*hmag / mu;
k = hhat(1)/(1.0 + hhat(3));
h = -hhat(2)/(1.0 + hhat(3));
kk = k*k;
hh = h*h;
s2 = 1.0 + hh + kk;
tkh = 2.0*k*h;
ecc = cross(v,hvec)/mu - rhat;
fhat(1) = 1.0 - kk + hh;
fhat(2) = tkh;
fhat(3) = -2.0*k;
ghat(1) = tkh;
ghat(2) = 1.0 + kk - hh;
ghat(3) = 2.0*h;
fhat = fhat/s2;
ghat = ghat/s2;
f = dot(ecc,fhat);
g = dot(ecc,ghat);
L = atan2(rhat(2)-vhat(1),rhat(1)+vhat(2));

mee = simplify([p,f,g,h,k,L]);

J = sym('J', [6 6]);
for i = 1:6
    for j = 1:6
        J(i,j) = simplify(diff(mee(i),x(j)));
    end
end

ccode(J, 'file', './cart_to_mee_partials.c');