function [] = showModelInfo(m)
    fprintf("===== Model Info =====\n");
    fprintf("Nz= %d\n", m.Nz);
    fprintf("H = %f\n", m.H);
    fprintf("K = %f\n", m.Kv);
    fprintf("dt = %f\n", m.dt);
end