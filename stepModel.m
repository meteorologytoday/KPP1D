function [ T_out, S_out ] = stepModel( m, T, S )

    % x_t+1 - x_t = dt * OP * x_t+1
    % (I - dt * OP) x_t+1 = x_t
    % x_t+1 = (I - dt*OP) \ x_t

    S_out = (m.grid.T_I_T - m.dt * m.op.diffz) \ S;
    T_out = T;
end