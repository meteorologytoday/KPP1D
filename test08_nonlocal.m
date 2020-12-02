clear;

Kv_iso = 1e-3;
Nz = 100;
H  = 25;    % m
dt = 10*60;
f = 1e-4;
ws_sfc = 1;
h_k = floor(Nz / 2);

m = Model(H, Nz, Kv_iso, dt, f);
m.showModelInfo();

wb_sfc_arr = [-0.1 0.05 0 1 0.05 0.1];
legend_str = {};
hold on;
for i=1:length(wb_sfc_arr)
    wb_sfc = wb_sfc_arr(i);
    nloc_flux = m.kpp.calNonLocalFlux_s(m.grid, h_k, wb_sfc, ws_sfc);
    plot(nloc_flux, m.grid.z_W);
    legend_str{end+1} = sprintf('wb_{sfc} = %f', wb_sfc);
end
legend(gca, legend_str, 'Location', 'south');
hold off;


