clear;

Kv_iso = 1e-3;
Nz = 100;
H  = 25; % m

dt = 10*60;

m = Model(H, Nz, Kv_iso, dt);

%%%%% plot the shape function in LMD94 figure 3 %%%%%
x = linspace(-1, 2, 300);
shape = m.kpp.shapeInterior(x);
plot(x, shape);
xlim([0 0.8]);

