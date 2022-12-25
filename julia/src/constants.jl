const g               = 9.80616     # m / s^2      copied from models/csm_share/shr/shr_const_mod.F90
const Lv_w            = 2.257e6     # Heat of evaporation of water in J/kg
const ρ_fw            = 1000.0 
const ρ_sw            = 1026.0
const cp_sw           = 3996.0
const ρ_a             = 1.22
const cp_a            = 1004.0
const σ_sb            = 5.670374419e-8 # Stefan-Boltzmann constant
const Cel_Kel_offset  = 273.15

const cp_ρ_sw = ρ_sw * cp_sw
const cp_ρ_a  = ρ_a  * cp_a

const ρ0 = ρ_sw

const BUOYANCY_TYPE = :LINEAR  # Possible options: [ :LINEAR, :NONLINEAR ]
