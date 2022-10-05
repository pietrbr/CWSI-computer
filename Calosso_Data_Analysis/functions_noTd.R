CWSI <- function(gamma, rC_rA, delta, gamma_star=0.66){
  num <- gamma*(1 + rC_rA) - gamma_star
  den <- delta + gamma*(1 + rC_rA)
  return(num/den)
}


rC_rA <- function(gamma, Ra, Rn, rho, cp, Tc, Ta, delta, pVsat, pV){
  num <- (gamma*Ra*Rn)/(rho*cp) - (Tc - Ta)*(delta + gamma) - (pVsat - pV)
  den <- gamma*((Tc - Ta) - (Ra*Rn)/(rho*cp))
  return(num/den)
}


gamma <- function(p, cp, lambdaV=2.62e6, mW=0.622){
  num <- cp*p
  den <- lambdaV*mW
  return(num/den)
}


# pVsat <- function(Ta){
#   Ta <- Ta - 273.15
#   return(6.108*exp((17.271*Ta)/(Ta + 237.3)))
# } # [hPa]
pVsat <- function(Tc){
  return(6.108*exp((17.271*Tc)/(Tc + 237.3)))
} # [hPa]


pv <- function(pVsat, Ur){
  Ur <- Ur/100
  return(pVsat*Ur)
} # [hPa]


#[Ur <- function(Ta, Td){
#  a <- 17.271
#  b <- 237.7
#  Ta <- Ta - 273.15
#  Td <- Td - 273.15
#  s1 <- Td/(b + Td)
#  s2 <- Ta/(b + Ta)
#  return(exp(a*(s1 - s2)))
#}


pas <- function(p, pV){
  return(p - pV)
}


rho <- function(pas, pV, Ta, Ras=287.058, Rv=461.495){
  Ta <- Ta + 273.15
  pas <- pas*100
  pV <- pV*100
  sum1 <- pas/(Ras*Ta)
  sum2 <- pV/(Rv*Ta)
  return(sum1 + sum2)
}


x <- function(p, pV, mw=0.622){
  return(mw*pV/(p - pV))
}


cp <- function(x, cpas=1006, cph20=1860){
  return(cpas + x*cph20)
}


delta <- function(Ta){
  p1 <- (105.49*237.3)/((Ta + 237.3)^2)
  p2 <- exp((17.271*Ta)/(Ta + 237.3))
  return(p1*p2)
}


Ra <- function(h, u2, z2=2){
  d <- 0.63*h
  z0 <- 0.13*h
  num <- 4.72*(log((z2 - d)/z0))^2
  den <- 1 + 0.54*u2
  return(num/den)
}

Em <- function(pv){
  return(0.52+0.065*sqrt(pv))
}


Rn <- function(UV, Em, Tc, Ta, r=0.19, Csis=0.98, sigB=5.67*1e-8){
  Tc <- Tc + 273.15
  return( UV*(1 - r) + Em*sigB*(Ta^4) - Csis*sigB*(Tc^4) )
}

#Rn <- function(UV, IR, Tc, r=0.19, Csis=0.98, sigB=5.67*1e-8){
#  IR <- IR/100
#  Tc <- Tc + 273.15
#  return( UV*(1 - r) + IR - Csis*sigB*(Tc^4) )
#}















