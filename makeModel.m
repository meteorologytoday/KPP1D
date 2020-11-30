function [ m ] = makeModel(H, Nz, Kv, dt)

    grid = makeGrid(H, Nz);
    sop = makeSpatialOperators(grid);
    
    op.diffz = Kv * sop.T_ddz2_T;
    
    m.grid = grid;
    m.sop = sop;
    m.H = H;
    m.Nz = Nz;
    m.Kv = Kv;
    m.dt = dt;
    m.op = op;
    

end