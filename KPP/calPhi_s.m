% LMD1994 Appendex B for scalars
function [ phi_s ] = calPhi_s(d, L_star)
    zeta = d / L_star;
    if (L_star >= 0)
        phi_s = 1 + 5 * zeta;
    else
        % 98.96 = c_s in KPPconstants
        phi_s = (1 - 16*zeta).^(-1/2) .* (zeta > -1.0) + ...
            (-28.86 - 98.96 * zeta).^(-1/3) .* (zeta <= -1.0);
    end
end