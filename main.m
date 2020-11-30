
Kv = 1e-3;
dt = 60;

total_steps = 86400 / dt;

Nz = 10;
H  = 25; % m


kpp = KPPConstants();

model = makeModel(H, Nz, Kv, dt);
showModelInfo(model);

T = 3+(1:Nz)';
S = (1:Nz)';

hold off;
figure;

for step=1:total_steps
    if mod(step, 60) == 1
        plot(S, model.grid.zt);
        hold on;

        fprintf("Integrate S: %f\n", sum( model.grid.dzt .* S) )
    end
    
    [T, S] = stepModel(model, T, S);

    
end



