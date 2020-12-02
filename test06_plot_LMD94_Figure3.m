clear;

Kv_iso = 1e-3;
Nz = 100;
H  = 25; % m
f = 1e-4;
dt = 10*60;

m = Model(H, Nz, Kv_iso, dt, f);

%%%%% plot the shape function in LMD94 figure 3 %%%%%
x = linspace(-1, 2, 300);
shape = m.kpp.shapeInterior(x);
plot(x, shape);
xlim([0 0.8]);

