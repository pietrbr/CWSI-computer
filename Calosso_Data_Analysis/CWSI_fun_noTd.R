# mettere la pressione in [hPa]

CWSI_fun_noTd <- function(Ta, Tc, UV, p, Ur, u2, h){
  source('functions_noTd.R')
  
  Ra_ <- Ra(h, u2)
  print(paste('Ra: ', Ra_))
  
  #Ur_ <- Ur(Ta, Td)           # cambio Ta -> Tc
  #Ur_ <- Ur(Tc, Td) 
  #print(paste('Ur: ', Ur_))
  #pVsat_ <- pVsat(Ta)        # cambio Ta -> Tc
  
  pVsat_ <- pVsat(Tc)
  print(paste('pVsat: ', pVsat_))
  
  pv_ <- pv(pVsat_, Ur)
  print(paste('pv: ', pv_))
  
  x_ <- x(p, pv_)
  print(paste('x: ', x_))
  
  cp_ <- cp(x_)
  print(paste('cp: ', cp_))
  
  gamma_ <- gamma(p, cp_)
  print(paste('gamma: ', gamma_))
  
  pas_ <- pas(p, pv_)
  print(paste('pas: ', pas_))
  
  #rho_ <- rho(pas_, pv_, Ta)        # cambio Ta -> Tc
  rho_ <- rho(pas_, pv_, Tc) 
  print(paste('rho: ', rho_))
  
  #delta_ <- delta(Ta)              # cambio Ta -> Tc 
  delta_ <- delta(Tc) 
  print(paste('delta: ', delta_))
  
  Em_ <- Em(pv_)
  print(paste('Em: ', Em_))
  
  
  Rn_ <- Rn(UV, Em_, Tc, Ta)
  print(paste('Rn: ', Rn_))
  
  Rc_Ra_ <- rC_rA(gamma_, Ra_, Rn_, rho_, cp_, Tc, Ta, delta_, pVsat_, pv_)
  print(paste('Rc_Ra: ', Rc_Ra_))
  
  CWSI_ <- CWSI(gamma_, Rc_Ra_, delta_)
  return(CWSI_)
}