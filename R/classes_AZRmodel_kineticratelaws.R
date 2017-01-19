###############################################################################
###############################################################################
# This file contains kinetic rate law functions
###############################################################################
###############################################################################

###############################################################################
# Return all names of the implemented kinetic rate laws
# Note that this function needs to be maintained manually if kinetic rate laws
# are added deleted or renamed. Main use of info: namespace adjustments in deSolve
# simulation functions.
###############################################################################
getAllKineticRateLaws <- function() {
  allRateLawsNames <- c('kin_allosteric_inihib_empirical_rev',
                        'kin_allosteric_inihib_mwc_irr',
                        'kin_catalytic_activation_irr',
                        'kin_catalytic_activation_rev',
                        'kin_comp_inihib_irr',
                        'kin_comp_inihib_rev',
                        'kin_constantflux',
                        'kin_degradation',
                        'kin_hill_1_modifier_rev',
                        'kin_hill_2_modifiers_rev',
                        'kin_hill_cooperativity_irr',
                        'kin_hill_rev',
                        'kin_hyperbolic_modifier_irr',
                        'kin_hyperbolic_modifier_rev',
                        'kin_iso_uni_uni_rev',
                        'kin_mass_action_irr',
                        'kin_mass_action_rev',
                        'kin_michaelis_menten_irr',
                        'kin_michaelis_menten_rev',
                        'kin_mixed_activation_irr',
                        'kin_mixed_activation_rev',
                        'kin_mixed_inihib_irr',
                        'kin_mixed_inihib_rev',
                        'kin_noncomp_inihib_irr',
                        'kin_noncomp_inihib_rev',
                        'kin_ordered_bi_bi_rev',
                        'kin_ordered_bi_uni_rev',
                        'kin_ordered_uni_bi_rev',
                        'kin_ping_pong_bi_bi_rev',
                        'kin_specific_activation_irr',
                        'kin_specific_activation_rev',
                        'kin_substrate_activation_irr',
                        'kin_substrate_inihib_irr',
                        'kin_substrate_inihib_rev',
                        'kin_uncomp_inihib_irr',
                        'kin_uncomp_inihib_rev',
                        'kin_uni_uni_rev')
  return(allRateLawsNames)
}

###############################################################################
# kin_mixed_activation_irr
###############################################################################
# Kinetic Rate Law: Mixed activation irreversible
#
# Needed for AZRmodels
#
# R = V*substrate*activator / ( Kms*(Kas+activator) + substrate*(Kac+activator) )
kin_mixed_activation_irr <- function(V,substrate,activator,Kms,Kas,Kac) {
  R <- V*substrate*activator / ( Kms*(Kas+activator) + substrate*(Kac+activator) )
  return(R)
}

###############################################################################
# kin_michaelis_menten_rev
###############################################################################
# Kinetic Rate Law: Michaelis Menten (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = (Vf*substrate/Kms+Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp )
kin_michaelis_menten_rev <- function(Vf,substrate,Kms,Vr,product,Kmp) {
  R <- (Vf*substrate/Kms+Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp )
  return(R)
}

###############################################################################
# kin_michaelis_menten_irr
###############################################################################
# Kinetic Rate Law: Michaelis Menten (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#    - Only irreversible
#    - exactly 1 substrate
#
# R = V*substrate / ( Km + substrate )
kin_michaelis_menten_irr <- function(V,substrate,Km) {
  R <- V*substrate / ( Km + substrate )
  return(R)
}

###############################################################################
# kin_mass_action_rev
###############################################################################
# Kinetic Rate Law: Mass action (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate (several substrates can be realized by
#     substrate=substrate1*substrate2*...
#   - exactly 1 product (several products can be realized by
#     product=product1*product2*...
#
# R = k1*substrate-k2*product
kin_mass_action_rev <- function(k1,substrate,k2,product) {
  R <- k1*substrate-k2*product
  return(R)
}

###############################################################################
# kin_mass_action_irr
###############################################################################
# Kinetic Rate Law: Mass action (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate (several substrates can be realized by
#     substrate=substrate1*substrate2*...
#
# R = k*substrate
kin_mass_action_irr <- function(k,substrate) {
  R <- k*substrate
  return(R)
}

###############################################################################
# kin_allosteric_inihib_empirical_rev
###############################################################################
# Kinetic Rate Law: Allosteric inhibition (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = (Vf*substrate/Kms - Vr*product/Kmp) / (1 + substrate/Kms + product/Kmp + (inhibitor/Ki)^n)
kin_allosteric_inihib_empirical_rev <- function(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Ki,n) {
  R <- (Vf*substrate/Kms - Vr*product/Kmp) / (1 + substrate/Kms + product/Kmp + (inhibitor/Ki)^n)
  return(R)
}

###############################################################################
# kin_allosteric_inihib_mwc_irr
###############################################################################
# Kinetic Rate Law: Allosteric inhibition (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*substrate*(Ks+substrate)^(n-1) / ( L*(Ks*(1+inhibitor/Ki))^n + (Ks+substrate)^n )
kin_allosteric_inihib_mwc_irr <- function(V,substrate,Ks,n,L,inhibitor,Ki) {
  R <- V*substrate*(Ks+substrate)^(n-1) / ( L*(Ks*(1+inhibitor/Ki))^n + (Ks+substrate)^n )
  return(R)
}

###############################################################################
# kin_catalytic_activation_irr
###############################################################################
# Kinetic Rate Law: Catalytic activation (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*substrate*activator / ( (Kms + substrate)*(Ka+activator) )
kin_catalytic_activation_irr <- function(V,substrate,activator,Kms,Ka) {
  R <- V*substrate*activator / ( (Kms + substrate)*(Ka+activator) )
  return(R)
}

###############################################################################
# kin_catalytic_activation_rev
###############################################################################
# Kinetic Rate Law: Catalytic activation (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = ( Vf*substrate/Kms - Vr*product/Kmp )*activator / ( (1+substrate/Kms+product/Kmp)*(Ka+activator) )
kin_catalytic_activation_rev <- function(Vf,substrate,Kms,Vr,product,Kmp,activator,Ka) {
  R <- ( Vf*substrate/Kms - Vr*product/Kmp )*activator / ( (1+substrate/Kms+product/Kmp)*(Ka+activator) )
  return(R)
}

###############################################################################
# kin_comp_inihib_irr
###############################################################################
# Kinetic Rate Law: Competitive inhibition (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*substrate / ( Km*(1+inhibitor/Ki) + substrate )
kin_comp_inihib_irr <- function(V,substrate,Km,inhibitor,Ki) {
  R <- V*substrate / ( Km*(1+inhibitor/Ki) + substrate )
  return(R)
}

###############################################################################
# kin_comp_inihib_rev
###############################################################################
# Kinetic Rate Law: Competitive inhibition (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = V*substrate / ( Km*(1+inhibitor/Ki) + substrate )
kin_comp_inihib_rev <- function(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Ki) {
  R <- (Vf*substrate/Kms - Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp+inhibitor/Ki )
  return(R)
}

###############################################################################
# kin_constantflux
###############################################################################
# Kinetic Rate Law: Constant flux rate
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Be careful that concentrations do not become negative!
#
# R = v
kin_constantflux <- function(v) {
  R <- v
  return(R)
}

###############################################################################
# kin_degradation
###############################################################################
# Kinetic Rate Law: Linear degradation kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = kdeg*substrate
kin_degradation <- function(kdeg,substrate) {
  R <- kdeg*substrate
  return(R)
}

###############################################################################
# kin_hill_1_modifier_rev
###############################################################################
# Kinetic Rate Law: Reversible Hill type kinetics with one modifier
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = Vf*substrate/Shalve*(1-product/(substrate*Keq))*(substrate/Shalve+product/Phalve)^(h-1) / ( (1+(modifier/Mhalve)^h)/(1+alpha*(modifier/Mhalve)^h)+(substrate/Shalve + product/Phalve)^h )
kin_hill_1_modifier_rev <- function(Vf,substrate,Shalve,product,Keq,Phalve,h,modifier,Mhalve,alpha) {
  R <- Vf*substrate/Shalve*(1-product/(substrate*Keq))*(substrate/Shalve+product/Phalve)^(h-1) / ( (1+(modifier/Mhalve)^h)/(1+alpha*(modifier/Mhalve)^h)+(substrate/Shalve + product/Phalve)^h )
  return(R)
}

###############################################################################
# kin_hill_2_modifiers_rev
###############################################################################
# Kinetic Rate Law: Reversible Hill type kinetics with two modifiers
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = Vf*substrate/Shalve*(1-product/(substrate*Keq))*(substrate/Shalve+product/Phalve)^(h-1) / ( (1+(modifierA/MAhalve)^h + 1+(modifierB/MBhalve)^h) / ( 1+alphaA*(modifierA/MAhalve)^h+alphaB*(modifierB/MBhalve)^h+alphaA*alphaB*alphaAB*(modifierA/MAhalve)^h*(modifierB/MBhalve)^h ) + (substrate/Shalve + product/Phalve)^h )
kin_hill_2_modifiers_rev <- function(Vf,substrate,Shalve,product,Keq,Phalve,h,modifierA,MAhalve,modifierB,MBhalve,alphaA,alphaB,alphaAB) {
  R <- Vf*substrate/Shalve*(1-product/(substrate*Keq))*(substrate/Shalve+product/Phalve)^(h-1) / ( (1+(modifierA/MAhalve)^h + 1+(modifierB/MBhalve)^h) / ( 1+alphaA*(modifierA/MAhalve)^h+alphaB*(modifierB/MBhalve)^h+alphaA*alphaB*alphaAB*(modifierA/MAhalve)^h*(modifierB/MBhalve)^h ) + (substrate/Shalve + product/Phalve)^h )
  return(R)
}

###############################################################################
# kin_hill_cooperativity_irr
###############################################################################
# Kinetic Rate Law: Hill type (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*substrate^h / ( Shalve^h + substrate^h )
kin_hill_cooperativity_irr <- function(V,substrate,h,Shalve) {
  R <- V*substrate^h / ( Shalve^h + substrate^h )
  return(R)
}

###############################################################################
# kin_hill_rev
###############################################################################
# Kinetic Rate Law: Hill type (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = Vf*substrate/Shalve*(1-product/(substrate*Keq))*(substrate/Shalve+product/Phalve)^(h-1) / ( 1+(substrate/Shalve + product/Phalve)^h )
kin_hill_rev <- function(Vf,substrate,Shalve,product,Keq,Phalve,h) {
  R <- Vf*substrate/Shalve*(1-product/(substrate*Keq))*(substrate/Shalve+product/Phalve)^(h-1) / ( 1+(substrate/Shalve + product/Phalve)^h )
  return(R)
}

###############################################################################
# kin_hyperbolic_modifier_irr
###############################################################################
# Kinetic Rate Law: Hyperbolic modifier (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*substrate*(1+b*modifier/(a*Kd)) / ( Km*(1+modifier/Kd) + substrate*(1+modifier/(a*Kd)) )
kin_hyperbolic_modifier_irr <- function(V,substrate,b,modifier,a,Kd,Km) {
  R <- V*substrate*(1+b*modifier/(a*Kd)) / ( Km*(1+modifier/Kd) + substrate*(1+modifier/(a*Kd)) )
  return(R)
}

###############################################################################
# kin_hyperbolic_modifier_rev
###############################################################################
# Kinetic Rate Law: Hyperbolic modifier (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = (Vf*substrate/Kms - Vr*product/Kmp)*(1+b*modifier/(a*Kd)) / ( 1+modifier/Kd+(substrate/Kms+product/Kmp)*(1+modifier/(a*Kd)) )
kin_hyperbolic_modifier_rev <- function(Vf,substrate,Kms,Vr,product,Kmp,b,modifier,a,Kd) {
  R <- (Vf*substrate/Kms - Vr*product/Kmp)*(1+b*modifier/(a*Kd)) / ( 1+modifier/Kd+(substrate/Kms+product/Kmp)*(1+modifier/(a*Kd)) )
  return(R)
}

###############################################################################
# kin_iso_uni_uni_rev
###############################################################################
# Kinetic Rate Law: enzyme isomerization product inhibition
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = Vf*(substrate-product/Keq) / ( substrate*(1+product/Kii) + Kms*(1+product/Kmp) )
kin_iso_uni_uni_rev <- function(Vf,substrate,product,Keq,Kii,Kms,Kmp) {
  R <- Vf*(substrate-product/Keq) / ( substrate*(1+product/Kii) + Kms*(1+product/Kmp) )
  return(R)
}

###############################################################################
# kin_mixed_activation_rev
###############################################################################
# Kinetic Rate Law: Mixed activation reversible
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = (Vf*substrate/Kms - Vr*product/Kmp)*activator / ( Kas+activator+(substrate/Kms+product/Kmp)*(Kac+activator) )
kin_mixed_activation_rev <- function(Vf,substrate,Kms,Vr,product,Kmp,activator,Kas,Kac) {
  R <- (Vf*substrate/Kms - Vr*product/Kmp)*activator / ( Kas+activator+(substrate/Kms+product/Kmp)*(Kac+activator) )
  return(R)
}

###############################################################################
# kin_mixed_inihib_irr
###############################################################################
# Kinetic Rate Law: Mixed inhibition (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*substrate / ( Km*(1+inhibitor/Kis) + substrate*(1+inhibitor/Kic) )
kin_mixed_inihib_irr <- function(V,substrate,Km,inhibitor,Kis,Kic) {
  R <- V*substrate / ( Km*(1+inhibitor/Kis) + substrate*(1+inhibitor/Kic) )
  return(R)
}

###############################################################################
# kin_mixed_inihib_rev
###############################################################################
# Kinetic Rate Law: Mixed inhibition (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = (Vf*substrate/Kms - Vr*product/Kmp) / ( 1+inhibitor/Kis+(substrate/Kms+product/Kmp)*(1+inhibitor/Kic) )
kin_mixed_inihib_rev <- function(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Kis,Kic) {
  R <- (Vf*substrate/Kms - Vr*product/Kmp) / ( 1+inhibitor/Kis+(substrate/Kms+product/Kmp)*(1+inhibitor/Kic) )
  return(R)
}

###############################################################################
# kin_noncomp_inihib_irr
###############################################################################
# Kinetic Rate Law: Noncompetitive inhibition (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*substrate / ( (Km+substrate)*(1+inhibitor/Ki) )
kin_noncomp_inihib_irr <- function(V,substrate,Km,inhibitor,Ki) {
  R <- V*substrate / ( (Km+substrate)*(1+inhibitor/Ki) )
  return(R)
}

###############################################################################
# kin_noncomp_inihib_rev
###############################################################################
# Kinetic Rate Law: Noncompetitive inhibition (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = (Vf*substrate/Kms-Vr*product/Kmp) / ( (1+substrate/Kms+product/Kmp)*(1+inhibitor/Ki) )
kin_noncomp_inihib_rev <- function(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Ki) {
  R <- (Vf*substrate/Kms-Vr*product/Kmp) / ( (1+substrate/Kms+product/Kmp)*(1+inhibitor/Ki) )
  return(R)
}

###############################################################################
# kin_ordered_bi_bi_rev
###############################################################################
# Kinetic Rate Law: ordered bi-bi reversible
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 2 substrates
#   - exactly 2 products
#
# R = Vf*(substratea*substrateb-productp*productq/Keq) / (substratea*substrateb*(1+productp/Kip) + Kma*substrateb + Kmb*(substratea+Kia)+Vf/(Vr*Keq)*(Kmq*productp*(1+substratea/Kia) + productq*(Kmp*(1+Kia*substrateb/(Kma*Kmb))+productp*(1+substrateb/Kib))) )
kin_ordered_bi_bi_rev <- function(Vf,substratea,substrateb,productp,productq,Keq,Kip,Kma,Kmb,Kia,Vr,Kmq,Kmp,Kib) {
  R <- Vf*(substratea*substrateb-productp*productq/Keq) / (substratea*substrateb*(1+productp/Kip) + Kma*substrateb + Kmb*(substratea+Kia)+Vf/(Vr*Keq)*(Kmq*productp*(1+substratea/Kia) + productq*(Kmp*(1+Kia*substrateb/(Kma*Kmb))+productp*(1+substrateb/Kib))) )
  return(R)
}

###############################################################################
# kin_ordered_bi_uni_rev
###############################################################################
# Kinetic Rate Law: ordered bi-uni reversible
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 2 substrates
#   - exactly 1 product
#
# R = Vf*(substratea*substrateb-product/Keq) / ( substratea*substrateb+Kma*substrateb+Kmb*substratea+Vf/(Vr*Keq)*(Kmp+product*(1+substratea/Kia)) )
kin_ordered_bi_uni_rev <- function(Vf,substratea,substrateb,product,Keq,Kma,Kmb,Vr,Kmp,Kia) {
  R <- Vf*(substratea*substrateb-product/Keq) / ( substratea*substrateb+Kma*substrateb+Kmb*substratea+Vf/(Vr*Keq)*(Kmp+product*(1+substratea/Kia)) )
  return(R)
}

###############################################################################
# kin_ordered_uni_bi_rev
###############################################################################
# Kinetic Rate Law: ordered uni-bi reversible
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 2 products
#
# R = Vf*(substrate-productp*productq/Keq) / ( Kms+substrate*(1+productp/Kip)+Vf/(Vr*Keq)*(Kmq*productp+Kmp*productq+productp*productq) )
kin_ordered_uni_bi_rev <- function(Vf,substrate,productp,productq,Keq,Kms,Kip,Vr,Kmq,Kmp) {
  R <- Vf*(substrate-productp*productq/Keq) / ( Kms+substrate*(1+productp/Kip)+Vf/(Vr*Keq)*(Kmq*productp+Kmp*productq+productp*productq) )
  return(R)
}

###############################################################################
# kin_ping_pong_bi_bi_rev
###############################################################################
# Kinetic Rate Law: Ping pong bi bi kinetics (reversible)
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 2 substrates
#   - exactly 2 products
#
# R = Vf*( substratea*substrateb-productp*productq/Keq ) / (substratea*substrateb*(1+productq/Kiq)+Kma*substrateb+Kmb*substratea+Vf/(Vr*Keq)*(Kmq*productp*(1+substratea/Kia)+productq*(Kmp+productp)))
kin_ping_pong_bi_bi_rev <- function(Vf,substratea,substrateb,productp,productq,Keq,Kiq,Kma,Kmb,Vr,Kmq,Kia,Kmp) {
  R <- Vf*( substratea*substrateb-productp*productq/Keq ) / (substratea*substrateb*(1+productq/Kiq)+Kma*substrateb+Kmb*substratea+Vf/(Vr*Keq)*(Kmq*productp*(1+substratea/Kia)+productq*(Kmp+productp)))
  return(R)
}

###############################################################################
# kin_specific_activation_irr
###############################################################################
# Kinetic Rate Law: Specific activation (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*substrate*activator/( Kms*Ka+(Kms+substrate)*activator )
kin_specific_activation_irr <- function(V,substrate,activator,Kms,Ka) {
  R <- V*substrate*activator/( Kms*Ka+(Kms+substrate)*activator )
  return(R)
}

###############################################################################
# kin_specific_activation_rev
###############################################################################
# Kinetic Rate Law: Specific activation (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = ( Vf*substrate/Kms-Vr*product/Kmp )*activator / ( Ka+(1+substrate/Kms+product/Kmp)*activator )
kin_specific_activation_rev <- function(Vf,substrate,Kms,Vr,product,Kmp,activator,Ka) {
  R <- ( Vf*substrate/Kms-Vr*product/Kmp )*activator / ( Ka+(1+substrate/Kms+product/Kmp)*activator )
  return(R)
}

###############################################################################
# kin_substrate_activation_irr
###############################################################################
# Kinetic Rate Law: Substrate activation (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*(substrate/Ksa)^2 / ( 1+substrate/Ksc+substrate/Ksa+(substrate/Ksa)^2 )
kin_substrate_activation_irr <- function(V,substrate,Ksa,Ksc) {
  R <- V*(substrate/Ksa)^2 / ( 1+substrate/Ksc+substrate/Ksa+(substrate/Ksa)^2 )
  return(R)
}

###############################################################################
# kin_substrate_inihib_irr
###############################################################################
# Kinetic Rate Law: Substrate inhibition (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*substrate / ( Km + substrate + Km*(substrate/Ki)^2 )
kin_substrate_inihib_irr <- function(V,substrate,Km,Ki) {
  R <- V*substrate / ( Km + substrate + Km*(substrate/Ki)^2 )
  return(R)
}

###############################################################################
# kin_substrate_inihib_rev
###############################################################################
# Kinetic Rate Law: Substrate inhibition (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = (Vf*substrate/Kms-Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp+(substrate/Ki)^2 )
kin_substrate_inihib_rev <- function(Vf,substrate,Kms,Vr,product,Kmp,Ki) {
  R <- (Vf*substrate/Kms-Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp+(substrate/Ki)^2 )
  return(R)
}

###############################################################################
# kin_uncomp_inihib_irr
###############################################################################
# Kinetic Rate Law: Uncompetitive inhibition (irreversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only irreversible
#   - exactly 1 substrate
#
# R = V*substrate / ( Km + substrate*(1+inhibitor/Ki) )
kin_uncomp_inihib_irr <- function(V,substrate,Km,inhibitor,Ki) {
  R <- V*substrate / ( Km + substrate*(1+inhibitor/Ki) )
  return(R)
}

###############################################################################
# kin_uncomp_inihib_rev
###############################################################################
# Kinetic Rate Law: Uncompetitive inhibition (reversible) kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = ( Vf*substrate/Kms-Vr*product/Kmp ) / ( 1+(substrate/Kms+product/Kmp)*(1+inhibitor/Ki) )
kin_uncomp_inihib_rev <- function(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Ki) {
  R <- ( Vf*substrate/Kms-Vr*product/Kmp ) / ( 1+(substrate/Kms+product/Kmp)*(1+inhibitor/Ki) )
  return(R)
}

###############################################################################
# kin_uni_uni_rev
###############################################################################
# Kinetic Rate Law: uni uni reversible kinetics
#
# Needed for AZRmodels
#
# Application restrictions:
# =========================
#   - Only reversible
#   - exactly 1 substrate
#   - exactly 1 product
#
# R = Vf*( substrate-product/Keq ) / ( substrate+Kms*(1+product/Kmp) )
kin_uni_uni_rev <- function(Vf,substrate,product,Keq,Kms,Kmp) {
  R <- Vf*( substrate-product/Keq ) / ( substrate+Kms*(1+product/Kmp) )
  return(R)
}








