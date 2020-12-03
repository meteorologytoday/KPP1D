clear;

Kv_iso = 1e-3;
Nz = 100;
H  = 25; % m
f = 1e-4;
dt = 10*60;

m = Model(H, Nz, Kv_iso, dt, f);
m.showModelInfo()

u = m.state.u * 0;
v = m.state.v * 0;
b = m.state.b * 0;

dbdz = 0.01;
N = sqrt(dbdz);
b = dbdz * m.grid.z_T;

shear_bnd = [-5 -15];
x = (m.grid.z_T - shear_bnd(1)) / (shear_bnd(2) - shear_bnd(1));
U = 1 - 3 * x.^2 + 2 * x.^3;
U(m.grid.z_T > shear_bnd(1)) = 1;
U(m.grid.z_T <= shear_bnd(2)) = 0;

angle = 30 * pi / 180;
u = U * cos(angle);
v = U * sin(angle);

x = (m.grid.z_W - shear_bnd(1)) / (shear_bnd(2) - shear_bnd(1));
k = 1 / (shear_bnd(2) - shear_bnd(1));
Ri_g_true = N^2 ./ ( k * (-6 * x + 6 * x.^2) ).^2;
Ri_g_true(m.grid.z_W > shear_bnd(1))  = inf;
Ri_g_true(m.grid.z_W <= shear_bnd(2)) = inf;


Ri_g = m.kpp.calGradientRichardsonNumber(m.grid, b, u, v);
K_s = m.kpp.calInteriorK_sh(m.grid, b, u, v);

subplot(1, 5, 1);
plot(b, m.grid.z_T);
title("b");

subplot(1, 5, 2);
plot(u, m.grid.z_T);
title("u");

subplot(1, 5, 3);
plot(v, m.grid.z_T);
title("v");

subplot(1, 5, 4);
hold on;
plot(Ri_g_true, m.grid.z_W, 'r-');
plot(Ri_g,      m.grid.z_W, 'k--');
title("Ri_g");
xlim([0 10])
ylim([-H 0])
legend(gca, {'Analytical', 'Numerical'}, 'Location', 'south');
hold off;

subplot(1, 5, 5);
plot(K_s, m.grid.z_W);
title("\nu^{shear}_s")
