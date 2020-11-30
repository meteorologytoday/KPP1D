addpath('./KPP');

tau0 = 0.1;
H = 10;
h = 5;
z = - linspace(0, 1, 1000) * h;
d = -z;
sig = d / h;
