addpath('./KPP');

tau0 = 0.1;
H = 10;
h = 5;
z = - linspace(0, 1, 1000) * H;
d = -z;

kpp = KPPConstants();
[ u_star, L_star ] = calMOSTscales(kpp, tau0, 0.1 * 0.0001);

phi_s_conv = calPhi_s(d, L_star);
phi_s_stab = calPhi_s(d, -L_star);

% This subplot should look like LMD94 Figure B1. scalar
subplot(111);
title(sprintf('u_* = %f, L_* = %f', u_star, L_star));
hold on;
plot(d / L_star, phi_s_conv);
plot(-d / L_star, phi_s_stab);
xlim([-2 .2]);
hold off;
