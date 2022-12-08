clear;

sf = SurfaceFlux();


% Examine C_D
U10_arr = linspace(1, 20, 100);
[C_D, C_T, C_q] = calCoefficients(sf, U10_arr, 1);

figure;
subplot(3,1,1);
plot(U10_arr, C_D);
title('C_D');

subplot(3,1,2);
plot(U10_arr, C_T);
title('C_T');

subplot(3,1,3);
plot(U10_arr, C_q);
title('C_q');

xlabel('U10 [ m/s ]');


% Plot: vary with U10
U10_arr = linspace(1, 20, 100);
angle = 30;
o = getEmptyObj();
for i = 1:length(U10_arr)
    U10 = U10_arr(i);
    
    u10 = U10 * cos(angle * pi / 180);
    v10 = U10 * sin(angle * pi / 180);
    T_a = 25;
    T_o = 20;
    q_a = 20 / 1e3;
    [ o.wu(end+1), o.wv(end+1), o.wT_sen(end+1), o.wT_R(end+1), o.wq(end+1) ] = sf.calSurfaceFluxes(u10, v10, T_a, q_a, T_o);
end
o.x = U10_arr;
o.xlabel='U10 [m/s]';
o.title='Vary with U10';
do_plot(o);

% Plot: vary with T_a
T_a_arr = linspace(1, 30, 100);
o = getEmptyObj();
for i = 1:length(T_a_arr)
    T_a = T_a_arr(i);
    
    u10 = 5;
    v10 = 0;
    T_o = 20;
    q_a = 20 / 1e3;
    [ o.wu(end+1), o.wv(end+1), o.wT_sen(end+1), o.wT_R(end+1), o.wq(end+1) ] = sf.calSurfaceFluxes(u10, v10, T_a, q_a, T_o);
end
o.x = T_a_arr;
o.xlabel='T_a [degC]';
o.title='Vary with T_a';
do_plot(o);


% Plot: vary with q_a
q_a_arr = linspace(0, 30, 100) / 1e3;
o = getEmptyObj();
for i = 1:length(q_a_arr)
    q_a = q_a_arr(i);
    
    u10 = 5;
    v10 = 0;
    T_a = 19;
    T_o = 20;
    
    [ o.wu(end+1), o.wv(end+1), o.wT_sen(end+1), o.wT_R(end+1), o.wq(end+1) ] = sf.calSurfaceFluxes(u10, v10, T_a, q_a, T_o);
end
o.x = q_a_arr;
o.xlabel='q_a [ g / m^3 ]';
o.title='Vary with q_a';
do_plot(o);


function o = getEmptyObj()
    o.wu = [];
    o.wv = [];
    o.wT_sen = [];
    o.wT_R = [];
    o.wq = [];
end


function [ ] = do_plot(o)
    figure;
    subplot(5,1,1);
    plot(o.x, o.wu);
    title('wu');

    subplot(5,1,2);
    plot(o.x, o.wv);
    title('wv');

    subplot(5,1,3);
    plot(o.x, o.wT_sen);
    title('wT_{sen}');

    subplot(5,1,4);
    plot(o.x, o.wT_R);
    title('wT_R');

    subplot(5,1,5);
    plot(o.x, o.wq);
    title('wq');
    
    xlabel(o.xlabel);
    suptitle(o.title);
end

