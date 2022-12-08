clear;

tau0 = 0.1;
H = 10;
h = 5;
z = - linspace(0, 1, 1000) * H;
d = -z;

kpp = KPP();

[ u_star, L_star ] = kpp.calMOSTscales(tau0, 0.1 * 0.0001);

phi_s_conv = kpp.calPhi_s(d, L_star);
phi_s_stab = kpp.calPhi_s(d, -L_star);

phi_m_conv = kpp.calPhi_m(d, L_star);
phi_m_stab = kpp.calPhi_m(d, -L_star);

% This subplot should look like LMD94 Figure B1. scalar
subplot(111);
title(sprintf('u_* = %f, L_* = %f', u_star, L_star));
hold on;
plot(d / L_star, phi_s_conv, 'k--');
plot(-d / L_star, phi_s_stab, 'k--');

plot(d / L_star, phi_m_conv, 'k-');
plot(-d / L_star, phi_m_stab, 'k-');
xlim([-2 .2]);
hold off;

