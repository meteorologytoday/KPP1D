

const Ri_c  = 0.3  # Critical Richardson number
const κ     = 0.4  # von Karman constant
const ϵ     = 0.1   # epsilon in LMD1994
const c_s   = 98.96   # LMD94 equation (B1) and (B2)
const c_m   = 8.38   # LMD94 equation (B1) and (B2)



const ρ_0  = ρ_sw   # seawater density
const β_T  = -0.2 
const C_v  = 1.6 

# Mininum unresolved shear.
# However, although it is seen in the Cvmix project source code but
# I cannot find the actual value anywhere in either LMD94 or Van Roekel et al. (2018)
const min_Vtr_sqr = 1e-4  

const C_star = 10                                # LMD94 equation (20)
const C_s    = C_star * κ * (c_s * κ * ϵ)^(1/3)  # LMD94 equation (20)
