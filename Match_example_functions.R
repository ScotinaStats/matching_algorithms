match5.results <- function(){
  rownames(data) <- 1:nrow(data)
  data$id <- 1:nrow(data)
  data$both.1 <- data$id %in% match12$index.treated &
    data$id %in% match13$index.treated &
    data$id %in% match14$index.treated & 
    data$id %in% match15$index.treated
  temp <- data[data$both.1 == "TRUE", ]
  m12 <- cbind(match12$index.treated, match12$index.control)
  m13 <- cbind(match13$index.treated, 
               match13$index.control + sum(data$treat == "Treatment 2"))
  m14 <- cbind(match14$index.treated, 
               match14$index.control + sum(data$treat == "Treatment 2") + 
                 sum(data$treat == "Treatment 3"))
  m15 <- cbind(match15$index.treated, 
               match15$index.control + sum(data$treat == "Treatment 2") + 
                 sum(data$treat == "Treatment 3") + 
                 sum(data$treat == "Treatment 4"))
  m12 <- m12[m12[,1] %in% rownames(temp), ]
  m13 <- m13[m13[,1] %in% rownames(temp), ]
  m14 <- m14[m14[,1] %in% rownames(temp), ]
  m15 <- m15[m15[,1] %in% rownames(temp), ]
  quintets <- cbind(m12[order(m12[,1]), ], m13[order(m13[,1]), ], 
                    m14[order(m14[,1]), ], m15[order(m15[,1]), ])
  quintets <- as.matrix(quintets[,c(1, 2, 4, 6, 8)])
  n.quint <- nrow(quintets)
  
  data.quintets <- rbind(data[as.vector(t(quintets)), ])
  data.matched <- rbind(data[quintets[,1], ], 
                        data[quintets[,2], ], 
                        data[quintets[,3], ], 
                        data[quintets[,4], ], 
                        data[quintets[,5], ])
  
  if (P == 5){
    data.matched <- with(data.matched, 
                         data.frame(id, X1, X2, X3, X4, X5,
                                    p1, p2, p3, p4, p5, treat, T1, T2, T3, T4, T5))
  }
  if (P == 10){
    data.matched <- with(data.matched, 
                         data.frame(id, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10,
                                    p1, p2, p3, p4, p5, treat, T1, T2, T3, T4, T5))
  }
  if (P == 20){
    data.matched <- with(data.matched, 
                         data.frame(id, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, 
                                    X11, X12, X13, X14, X15, X16, X17, X18, X19, X20,
                                    p1, p2, p3, p4, p5, treat, T1, T2, T3, T4, T5))
  }
  data.unique <- unique(data.matched)
  
  n1m <- sum(data.unique$treat == "Treatment 1")
  n2m <- sum(data.unique$treat == "Treatment 2")
  n3m <- sum(data.unique$treat == "Treatment 3")
  n4m <- sum(data.unique$treat == "Treatment 4")
  n5m <- sum(data.unique$treat == "Treatment 5")
  Nm <- n1m + n2m + n3m + n4m + n5m
  
  # Matching performance
  data.unique$psi <- c()
  for (i in 1:Nm){
    data.unique$psi[i] <- sum(as.vector(quintets) == data.unique$id[i])
  }
  
  if (P == 5){
    Xbar.11 <- (1/n.quint)*with(data.unique, sum(X1*T1*psi))
    Xbar.21 <- (1/n.quint)*with(data.unique, sum(X2*T1*psi))
    Xbar.31 <- (1/n.quint)*with(data.unique, sum(X3*T1*psi))
    Xbar.41 <- (1/n.quint)*with(data.unique, sum(X4*T1*psi))
    Xbar.51 <- (1/n.quint)*with(data.unique, sum(X5*T1*psi))
    
    Xbar.12 <- (1/n.quint)*with(data.unique, sum(X1*T2*psi))
    Xbar.22 <- (1/n.quint)*with(data.unique, sum(X2*T2*psi))
    Xbar.32 <- (1/n.quint)*with(data.unique, sum(X3*T2*psi))
    Xbar.42 <- (1/n.quint)*with(data.unique, sum(X4*T2*psi))
    Xbar.52 <- (1/n.quint)*with(data.unique, sum(X5*T2*psi))
    
    Xbar.13 <- (1/n.quint)*with(data.unique, sum(X1*T3*psi))
    Xbar.23 <- (1/n.quint)*with(data.unique, sum(X2*T3*psi))
    Xbar.33 <- (1/n.quint)*with(data.unique, sum(X3*T3*psi))
    Xbar.43 <- (1/n.quint)*with(data.unique, sum(X4*T3*psi))
    Xbar.53 <- (1/n.quint)*with(data.unique, sum(X5*T3*psi))
    
    Xbar.14 <- (1/n.quint)*with(data.unique, sum(X1*T4*psi))
    Xbar.24 <- (1/n.quint)*with(data.unique, sum(X2*T4*psi))
    Xbar.34 <- (1/n.quint)*with(data.unique, sum(X3*T4*psi))
    Xbar.44 <- (1/n.quint)*with(data.unique, sum(X4*T4*psi))
    Xbar.54 <- (1/n.quint)*with(data.unique, sum(X5*T4*psi))
    
    Xbar.15 <- (1/n.quint)*with(data.unique, sum(X1*T5*psi))
    Xbar.25 <- (1/n.quint)*with(data.unique, sum(X2*T5*psi))
    Xbar.35 <- (1/n.quint)*with(data.unique, sum(X3*T5*psi))
    Xbar.45 <- (1/n.quint)*with(data.unique, sum(X4*T5*psi))
    Xbar.55 <- (1/n.quint)*with(data.unique, sum(X5*T5*psi))
    
    sb112 <- with(data.unique, (Xbar.11 - Xbar.12)/sd(data$X1[data$treat == "Treatment 1"]))
    sb212 <- with(data.unique, (Xbar.21 - Xbar.22)/sd(data$X2[data$treat == "Treatment 1"]))
    sb312 <- with(data.unique, (Xbar.31 - Xbar.32)/sd(data$X3[data$treat == "Treatment 1"]))
    sb412 <- with(data.unique, (Xbar.41 - Xbar.42)/sd(data$X4[data$treat == "Treatment 1"]))
    sb512 <- with(data.unique, (Xbar.51 - Xbar.52)/sd(data$X5[data$treat == "Treatment 1"]))
    
    sb113 <- with(data.unique, (Xbar.11 - Xbar.13)/sd(data$X1[data$treat == "Treatment 1"]))
    sb213 <- with(data.unique, (Xbar.21 - Xbar.23)/sd(data$X2[data$treat == "Treatment 1"]))
    sb313 <- with(data.unique, (Xbar.31 - Xbar.33)/sd(data$X3[data$treat == "Treatment 1"]))
    sb413 <- with(data.unique, (Xbar.41 - Xbar.43)/sd(data$X4[data$treat == "Treatment 1"]))
    sb513 <- with(data.unique, (Xbar.51 - Xbar.53)/sd(data$X5[data$treat == "Treatment 1"]))
    
    sb114 <- with(data.unique, (Xbar.11 - Xbar.14)/sd(data$X1[data$treat == "Treatment 1"]))
    sb214 <- with(data.unique, (Xbar.21 - Xbar.24)/sd(data$X2[data$treat == "Treatment 1"]))
    sb314 <- with(data.unique, (Xbar.31 - Xbar.34)/sd(data$X3[data$treat == "Treatment 1"]))
    sb414 <- with(data.unique, (Xbar.41 - Xbar.44)/sd(data$X4[data$treat == "Treatment 1"]))
    sb514 <- with(data.unique, (Xbar.51 - Xbar.54)/sd(data$X5[data$treat == "Treatment 1"]))
    
    sb115 <- with(data.unique, (Xbar.11 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb215 <- with(data.unique, (Xbar.21 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb315 <- with(data.unique, (Xbar.31 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb415 <- with(data.unique, (Xbar.41 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb515 <- with(data.unique, (Xbar.51 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    
    sb123 <- with(data.unique, (Xbar.12 - Xbar.13)/sd(data$X1[data$treat == "Treatment 1"]))
    sb223 <- with(data.unique, (Xbar.22 - Xbar.23)/sd(data$X2[data$treat == "Treatment 1"]))
    sb323 <- with(data.unique, (Xbar.32 - Xbar.33)/sd(data$X3[data$treat == "Treatment 1"]))
    sb423 <- with(data.unique, (Xbar.42 - Xbar.43)/sd(data$X4[data$treat == "Treatment 1"]))
    sb523 <- with(data.unique, (Xbar.52 - Xbar.53)/sd(data$X5[data$treat == "Treatment 1"]))
    
    sb124 <- with(data.unique, (Xbar.12 - Xbar.14)/sd(data$X1[data$treat == "Treatment 1"]))
    sb224 <- with(data.unique, (Xbar.22 - Xbar.24)/sd(data$X2[data$treat == "Treatment 1"]))
    sb324 <- with(data.unique, (Xbar.32 - Xbar.34)/sd(data$X3[data$treat == "Treatment 1"]))
    sb424 <- with(data.unique, (Xbar.42 - Xbar.44)/sd(data$X4[data$treat == "Treatment 1"]))
    sb524 <- with(data.unique, (Xbar.52 - Xbar.54)/sd(data$X5[data$treat == "Treatment 1"]))
    
    sb125 <- with(data.unique, (Xbar.12 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb225 <- with(data.unique, (Xbar.22 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb325 <- with(data.unique, (Xbar.32 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb425 <- with(data.unique, (Xbar.42 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb525 <- with(data.unique, (Xbar.52 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    
    sb134 <- with(data.unique, (Xbar.13 - Xbar.14)/sd(data$X1[data$treat == "Treatment 1"]))
    sb234 <- with(data.unique, (Xbar.23 - Xbar.24)/sd(data$X2[data$treat == "Treatment 1"]))
    sb334 <- with(data.unique, (Xbar.33 - Xbar.34)/sd(data$X3[data$treat == "Treatment 1"]))
    sb434 <- with(data.unique, (Xbar.43 - Xbar.44)/sd(data$X4[data$treat == "Treatment 1"]))
    sb534 <- with(data.unique, (Xbar.53 - Xbar.54)/sd(data$X5[data$treat == "Treatment 1"]))
    
    sb135 <- with(data.unique, (Xbar.13 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb235 <- with(data.unique, (Xbar.23 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb335 <- with(data.unique, (Xbar.33 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb435 <- with(data.unique, (Xbar.43 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb535 <- with(data.unique, (Xbar.53 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    
    sb145 <- with(data.unique, (Xbar.14 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb245 <- with(data.unique, (Xbar.24 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb345 <- with(data.unique, (Xbar.34 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb445 <- with(data.unique, (Xbar.44 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb545 <- with(data.unique, (Xbar.54 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    
    max2sb1 <- max(abs(sb112), abs(sb113), abs(sb114), abs(sb115),
                   abs(sb123), abs(sb124), abs(sb125), 
                   abs(sb134), abs(sb135), abs(sb145))
    max2sb2 <- max(abs(sb212), abs(sb213), abs(sb214), abs(sb215),
                   abs(sb223), abs(sb224), abs(sb225), 
                   abs(sb234), abs(sb235), abs(sb245))
    max2sb3 <- max(abs(sb312), abs(sb313), abs(sb314), abs(sb315),
                   abs(sb323), abs(sb324), abs(sb325), 
                   abs(sb334), abs(sb335), abs(sb345))
    max2sb4 <- max(abs(sb412), abs(sb413), abs(sb414), abs(sb415),
                   abs(sb423), abs(sb424), abs(sb425), 
                   abs(sb434), abs(sb435), abs(sb445))
    max2sb5 <- max(abs(sb512), abs(sb513), abs(sb514), abs(sb515),
                   abs(sb523), abs(sb524), abs(sb525), 
                   abs(sb534), abs(sb535), abs(sb545))
    
    max2sb <- mean(c(max2sb1, max2sb2, max2sb3, max2sb4, max2sb5))
    maxmax2sb <- max(c(max2sb1, max2sb2, max2sb3, max2sb4, max2sb5))
  }
  if (P == 10){
    Xbar.11 <- (1/n.quint)*with(data.unique, sum(X1*T1*psi))
    Xbar.21 <- (1/n.quint)*with(data.unique, sum(X2*T1*psi))
    Xbar.31 <- (1/n.quint)*with(data.unique, sum(X3*T1*psi))
    Xbar.41 <- (1/n.quint)*with(data.unique, sum(X4*T1*psi))
    Xbar.51 <- (1/n.quint)*with(data.unique, sum(X5*T1*psi))
    Xbar.61 <- (1/n.quint)*with(data.unique, sum(X6*T1*psi))
    Xbar.71 <- (1/n.quint)*with(data.unique, sum(X7*T1*psi))
    Xbar.81 <- (1/n.quint)*with(data.unique, sum(X8*T1*psi))
    Xbar.91 <- (1/n.quint)*with(data.unique, sum(X9*T1*psi))
    Xbar.101 <- (1/n.quint)*with(data.unique, sum(X10*T1*psi))
    
    Xbar.12 <- (1/n.quint)*with(data.unique, sum(X1*T2*psi))
    Xbar.22 <- (1/n.quint)*with(data.unique, sum(X2*T2*psi))
    Xbar.32 <- (1/n.quint)*with(data.unique, sum(X3*T2*psi))
    Xbar.42 <- (1/n.quint)*with(data.unique, sum(X4*T2*psi))
    Xbar.52 <- (1/n.quint)*with(data.unique, sum(X5*T2*psi))
    Xbar.62 <- (1/n.quint)*with(data.unique, sum(X6*T2*psi))
    Xbar.72 <- (1/n.quint)*with(data.unique, sum(X7*T2*psi))
    Xbar.82 <- (1/n.quint)*with(data.unique, sum(X8*T2*psi))
    Xbar.92 <- (1/n.quint)*with(data.unique, sum(X9*T2*psi))
    Xbar.102 <- (1/n.quint)*with(data.unique, sum(X10*T2*psi))
    
    Xbar.13 <- (1/n.quint)*with(data.unique, sum(X1*T3*psi))
    Xbar.23 <- (1/n.quint)*with(data.unique, sum(X2*T3*psi))
    Xbar.33 <- (1/n.quint)*with(data.unique, sum(X3*T3*psi))
    Xbar.43 <- (1/n.quint)*with(data.unique, sum(X4*T3*psi))
    Xbar.53 <- (1/n.quint)*with(data.unique, sum(X5*T3*psi))
    Xbar.63 <- (1/n.quint)*with(data.unique, sum(X6*T3*psi))
    Xbar.73 <- (1/n.quint)*with(data.unique, sum(X7*T3*psi))
    Xbar.83 <- (1/n.quint)*with(data.unique, sum(X8*T3*psi))
    Xbar.93 <- (1/n.quint)*with(data.unique, sum(X9*T3*psi))
    Xbar.103 <- (1/n.quint)*with(data.unique, sum(X10*T3*psi))
    
    Xbar.14 <- (1/n.quint)*with(data.unique, sum(X1*T4*psi))
    Xbar.24 <- (1/n.quint)*with(data.unique, sum(X2*T4*psi))
    Xbar.34 <- (1/n.quint)*with(data.unique, sum(X3*T4*psi))
    Xbar.44 <- (1/n.quint)*with(data.unique, sum(X4*T4*psi))
    Xbar.54 <- (1/n.quint)*with(data.unique, sum(X5*T4*psi))
    Xbar.64 <- (1/n.quint)*with(data.unique, sum(X6*T4*psi))
    Xbar.74 <- (1/n.quint)*with(data.unique, sum(X7*T4*psi))
    Xbar.84 <- (1/n.quint)*with(data.unique, sum(X8*T4*psi))
    Xbar.94 <- (1/n.quint)*with(data.unique, sum(X9*T4*psi))
    Xbar.104 <- (1/n.quint)*with(data.unique, sum(X10*T4*psi))
    
    Xbar.15 <- (1/n.quint)*with(data.unique, sum(X1*T5*psi))
    Xbar.25 <- (1/n.quint)*with(data.unique, sum(X2*T5*psi))
    Xbar.35 <- (1/n.quint)*with(data.unique, sum(X3*T5*psi))
    Xbar.45 <- (1/n.quint)*with(data.unique, sum(X4*T5*psi))
    Xbar.55 <- (1/n.quint)*with(data.unique, sum(X5*T5*psi))
    Xbar.65 <- (1/n.quint)*with(data.unique, sum(X6*T5*psi))
    Xbar.75 <- (1/n.quint)*with(data.unique, sum(X7*T5*psi))
    Xbar.85 <- (1/n.quint)*with(data.unique, sum(X8*T5*psi))
    Xbar.95 <- (1/n.quint)*with(data.unique, sum(X9*T5*psi))
    Xbar.105 <- (1/n.quint)*with(data.unique, sum(X10*T5*psi))
    
    sb112 <- with(data.unique, (Xbar.11 - Xbar.12)/sd(data$X1[data$treat == "Treatment 1"]))
    sb212 <- with(data.unique, (Xbar.21 - Xbar.22)/sd(data$X2[data$treat == "Treatment 1"]))
    sb312 <- with(data.unique, (Xbar.31 - Xbar.32)/sd(data$X3[data$treat == "Treatment 1"]))
    sb412 <- with(data.unique, (Xbar.41 - Xbar.42)/sd(data$X4[data$treat == "Treatment 1"]))
    sb512 <- with(data.unique, (Xbar.51 - Xbar.52)/sd(data$X5[data$treat == "Treatment 1"]))
    sb612 <- with(data.unique, (Xbar.61 - Xbar.62)/sd(data$X6[data$treat == "Treatment 1"]))
    sb712 <- with(data.unique, (Xbar.71 - Xbar.72)/sd(data$X7[data$treat == "Treatment 1"]))
    sb812 <- with(data.unique, (Xbar.81 - Xbar.82)/sd(data$X8[data$treat == "Treatment 1"]))
    sb912 <- with(data.unique, (Xbar.91 - Xbar.92)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1012 <- with(data.unique, (Xbar.101 - Xbar.102)/sd(data$X10[data$treat == "Treatment 1"]))
    
    sb113 <- with(data.unique, (Xbar.11 - Xbar.13)/sd(data$X1[data$treat == "Treatment 1"]))
    sb213 <- with(data.unique, (Xbar.21 - Xbar.23)/sd(data$X2[data$treat == "Treatment 1"]))
    sb313 <- with(data.unique, (Xbar.31 - Xbar.33)/sd(data$X3[data$treat == "Treatment 1"]))
    sb413 <- with(data.unique, (Xbar.41 - Xbar.43)/sd(data$X4[data$treat == "Treatment 1"]))
    sb513 <- with(data.unique, (Xbar.51 - Xbar.53)/sd(data$X5[data$treat == "Treatment 1"]))
    sb613 <- with(data.unique, (Xbar.61 - Xbar.63)/sd(data$X6[data$treat == "Treatment 1"]))
    sb713 <- with(data.unique, (Xbar.71 - Xbar.73)/sd(data$X7[data$treat == "Treatment 1"]))
    sb813 <- with(data.unique, (Xbar.81 - Xbar.83)/sd(data$X8[data$treat == "Treatment 1"]))
    sb913 <- with(data.unique, (Xbar.91 - Xbar.93)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1013 <- with(data.unique, (Xbar.101 - Xbar.103)/sd(data$X10[data$treat == "Treatment 1"]))
    
    sb114 <- with(data.unique, (Xbar.11 - Xbar.14)/sd(data$X1[data$treat == "Treatment 1"]))
    sb214 <- with(data.unique, (Xbar.21 - Xbar.24)/sd(data$X2[data$treat == "Treatment 1"]))
    sb314 <- with(data.unique, (Xbar.31 - Xbar.34)/sd(data$X3[data$treat == "Treatment 1"]))
    sb414 <- with(data.unique, (Xbar.41 - Xbar.44)/sd(data$X4[data$treat == "Treatment 1"]))
    sb514 <- with(data.unique, (Xbar.51 - Xbar.54)/sd(data$X5[data$treat == "Treatment 1"]))
    sb614 <- with(data.unique, (Xbar.61 - Xbar.64)/sd(data$X6[data$treat == "Treatment 1"]))
    sb714 <- with(data.unique, (Xbar.71 - Xbar.74)/sd(data$X7[data$treat == "Treatment 1"]))
    sb814 <- with(data.unique, (Xbar.81 - Xbar.84)/sd(data$X8[data$treat == "Treatment 1"]))
    sb914 <- with(data.unique, (Xbar.91 - Xbar.94)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1014 <- with(data.unique, (Xbar.101 - Xbar.104)/sd(data$X10[data$treat == "Treatment 1"]))
    
    sb115 <- with(data.unique, (Xbar.11 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb215 <- with(data.unique, (Xbar.21 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb315 <- with(data.unique, (Xbar.31 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb415 <- with(data.unique, (Xbar.41 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb515 <- with(data.unique, (Xbar.51 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    sb615 <- with(data.unique, (Xbar.61 - Xbar.65)/sd(data$X6[data$treat == "Treatment 1"]))
    sb715 <- with(data.unique, (Xbar.71 - Xbar.75)/sd(data$X7[data$treat == "Treatment 1"]))
    sb815 <- with(data.unique, (Xbar.81 - Xbar.85)/sd(data$X8[data$treat == "Treatment 1"]))
    sb915 <- with(data.unique, (Xbar.91 - Xbar.95)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1015 <- with(data.unique, (Xbar.101 - Xbar.105)/sd(data$X10[data$treat == "Treatment 1"]))
    
    sb123 <- with(data.unique, (Xbar.12 - Xbar.13)/sd(data$X1[data$treat == "Treatment 1"]))
    sb223 <- with(data.unique, (Xbar.22 - Xbar.23)/sd(data$X2[data$treat == "Treatment 1"]))
    sb323 <- with(data.unique, (Xbar.32 - Xbar.33)/sd(data$X3[data$treat == "Treatment 1"]))
    sb423 <- with(data.unique, (Xbar.42 - Xbar.43)/sd(data$X4[data$treat == "Treatment 1"]))
    sb523 <- with(data.unique, (Xbar.52 - Xbar.53)/sd(data$X5[data$treat == "Treatment 1"]))
    sb623 <- with(data.unique, (Xbar.62 - Xbar.63)/sd(data$X6[data$treat == "Treatment 1"]))
    sb723 <- with(data.unique, (Xbar.72 - Xbar.73)/sd(data$X7[data$treat == "Treatment 1"]))
    sb823 <- with(data.unique, (Xbar.82 - Xbar.83)/sd(data$X8[data$treat == "Treatment 1"]))
    sb923 <- with(data.unique, (Xbar.92 - Xbar.93)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1023 <- with(data.unique, (Xbar.102 - Xbar.103)/sd(data$X10[data$treat == "Treatment 1"]))
    
    sb124 <- with(data.unique, (Xbar.12 - Xbar.14)/sd(data$X1[data$treat == "Treatment 1"]))
    sb224 <- with(data.unique, (Xbar.22 - Xbar.24)/sd(data$X2[data$treat == "Treatment 1"]))
    sb324 <- with(data.unique, (Xbar.32 - Xbar.34)/sd(data$X3[data$treat == "Treatment 1"]))
    sb424 <- with(data.unique, (Xbar.42 - Xbar.44)/sd(data$X4[data$treat == "Treatment 1"]))
    sb524 <- with(data.unique, (Xbar.52 - Xbar.54)/sd(data$X5[data$treat == "Treatment 1"]))
    sb624 <- with(data.unique, (Xbar.62 - Xbar.64)/sd(data$X6[data$treat == "Treatment 1"]))
    sb724 <- with(data.unique, (Xbar.72 - Xbar.74)/sd(data$X7[data$treat == "Treatment 1"]))
    sb824 <- with(data.unique, (Xbar.82 - Xbar.84)/sd(data$X8[data$treat == "Treatment 1"]))
    sb924 <- with(data.unique, (Xbar.92 - Xbar.94)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1024 <- with(data.unique, (Xbar.102 - Xbar.104)/sd(data$X10[data$treat == "Treatment 1"]))
    
    sb125 <- with(data.unique, (Xbar.12 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb225 <- with(data.unique, (Xbar.22 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb325 <- with(data.unique, (Xbar.32 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb425 <- with(data.unique, (Xbar.42 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb525 <- with(data.unique, (Xbar.52 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    sb625 <- with(data.unique, (Xbar.62 - Xbar.65)/sd(data$X6[data$treat == "Treatment 1"]))
    sb725 <- with(data.unique, (Xbar.72 - Xbar.75)/sd(data$X7[data$treat == "Treatment 1"]))
    sb825 <- with(data.unique, (Xbar.82 - Xbar.85)/sd(data$X8[data$treat == "Treatment 1"]))
    sb925 <- with(data.unique, (Xbar.92 - Xbar.95)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1025 <- with(data.unique, (Xbar.102 - Xbar.105)/sd(data$X10[data$treat == "Treatment 1"]))
    
    sb134 <- with(data.unique, (Xbar.13 - Xbar.14)/sd(data$X1[data$treat == "Treatment 1"]))
    sb234 <- with(data.unique, (Xbar.23 - Xbar.24)/sd(data$X2[data$treat == "Treatment 1"]))
    sb334 <- with(data.unique, (Xbar.33 - Xbar.34)/sd(data$X3[data$treat == "Treatment 1"]))
    sb434 <- with(data.unique, (Xbar.43 - Xbar.44)/sd(data$X4[data$treat == "Treatment 1"]))
    sb534 <- with(data.unique, (Xbar.53 - Xbar.54)/sd(data$X5[data$treat == "Treatment 1"]))
    sb634 <- with(data.unique, (Xbar.63 - Xbar.64)/sd(data$X6[data$treat == "Treatment 1"]))
    sb734 <- with(data.unique, (Xbar.73 - Xbar.74)/sd(data$X7[data$treat == "Treatment 1"]))
    sb834 <- with(data.unique, (Xbar.83 - Xbar.84)/sd(data$X8[data$treat == "Treatment 1"]))
    sb934 <- with(data.unique, (Xbar.93 - Xbar.94)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1034 <- with(data.unique, (Xbar.103 - Xbar.104)/sd(data$X10[data$treat == "Treatment 1"]))
    
    sb135 <- with(data.unique, (Xbar.13 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb235 <- with(data.unique, (Xbar.23 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb335 <- with(data.unique, (Xbar.33 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb435 <- with(data.unique, (Xbar.43 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb535 <- with(data.unique, (Xbar.53 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    sb635 <- with(data.unique, (Xbar.63 - Xbar.65)/sd(data$X6[data$treat == "Treatment 1"]))
    sb735 <- with(data.unique, (Xbar.73 - Xbar.75)/sd(data$X7[data$treat == "Treatment 1"]))
    sb835 <- with(data.unique, (Xbar.83 - Xbar.85)/sd(data$X8[data$treat == "Treatment 1"]))
    sb935 <- with(data.unique, (Xbar.93 - Xbar.95)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1035 <- with(data.unique, (Xbar.103 - Xbar.105)/sd(data$X10[data$treat == "Treatment 1"]))
    
    sb145 <- with(data.unique, (Xbar.14 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb245 <- with(data.unique, (Xbar.24 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb345 <- with(data.unique, (Xbar.34 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb445 <- with(data.unique, (Xbar.44 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb545 <- with(data.unique, (Xbar.54 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    sb645 <- with(data.unique, (Xbar.64 - Xbar.65)/sd(data$X6[data$treat == "Treatment 1"]))
    sb745 <- with(data.unique, (Xbar.74 - Xbar.75)/sd(data$X7[data$treat == "Treatment 1"]))
    sb845 <- with(data.unique, (Xbar.84 - Xbar.85)/sd(data$X8[data$treat == "Treatment 1"]))
    sb945 <- with(data.unique, (Xbar.94 - Xbar.95)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1045 <- with(data.unique, (Xbar.104 - Xbar.105)/sd(data$X10[data$treat == "Treatment 1"]))
    
    max2sb1 <- max(abs(sb112), abs(sb113), abs(sb114), abs(sb115),
                   abs(sb123), abs(sb124), abs(sb125), 
                   abs(sb134), abs(sb135), abs(sb145))
    max2sb2 <- max(abs(sb212), abs(sb213), abs(sb214), abs(sb215),
                   abs(sb223), abs(sb224), abs(sb225), 
                   abs(sb234), abs(sb235), abs(sb245))
    max2sb3 <- max(abs(sb312), abs(sb313), abs(sb314), abs(sb315),
                   abs(sb323), abs(sb324), abs(sb325), 
                   abs(sb334), abs(sb335), abs(sb345))
    max2sb4 <- max(abs(sb412), abs(sb413), abs(sb414), abs(sb415),
                   abs(sb423), abs(sb424), abs(sb425), 
                   abs(sb434), abs(sb435), abs(sb445))
    max2sb5 <- max(abs(sb512), abs(sb513), abs(sb514), abs(sb515),
                   abs(sb523), abs(sb524), abs(sb525), 
                   abs(sb534), abs(sb535), abs(sb545))
    max2sb6 <- max(abs(sb612), abs(sb613), abs(sb614), abs(sb615),
                   abs(sb623), abs(sb624), abs(sb625), 
                   abs(sb634), abs(sb635), abs(sb645))
    max2sb7 <- max(abs(sb712), abs(sb713), abs(sb714), abs(sb715),
                   abs(sb723), abs(sb724), abs(sb725), 
                   abs(sb734), abs(sb735), abs(sb745))
    max2sb8 <- max(abs(sb812), abs(sb813), abs(sb814), abs(sb815),
                   abs(sb823), abs(sb824), abs(sb825), 
                   abs(sb834), abs(sb835), abs(sb845))
    max2sb9 <- max(abs(sb912), abs(sb913), abs(sb914), abs(sb915),
                   abs(sb923), abs(sb924), abs(sb925), 
                   abs(sb934), abs(sb935), abs(sb945))
    max2sb10 <- max(abs(sb1012), abs(sb1013), abs(sb1014), abs(sb1015),
                    abs(sb1023), abs(sb1024), abs(sb1025), 
                    abs(sb1034), abs(sb1035), abs(sb1045))
    
    max2sb <- mean(c(max2sb1, max2sb2, max2sb3, max2sb4, max2sb5,
                     max2sb6, max2sb7, max2sb8, max2sb9, max2sb10))
    maxmax2sb <- max(c(max2sb1, max2sb2, max2sb3, max2sb4, max2sb5,
                       max2sb6, max2sb7, max2sb8, max2sb9, max2sb10))
  }
  if (P == 20){
    Xbar.11 <- (1/n.quint)*with(data.unique, sum(X1*T1*psi))
    Xbar.21 <- (1/n.quint)*with(data.unique, sum(X2*T1*psi))
    Xbar.31 <- (1/n.quint)*with(data.unique, sum(X3*T1*psi))
    Xbar.41 <- (1/n.quint)*with(data.unique, sum(X4*T1*psi))
    Xbar.51 <- (1/n.quint)*with(data.unique, sum(X5*T1*psi))
    Xbar.61 <- (1/n.quint)*with(data.unique, sum(X6*T1*psi))
    Xbar.71 <- (1/n.quint)*with(data.unique, sum(X7*T1*psi))
    Xbar.81 <- (1/n.quint)*with(data.unique, sum(X8*T1*psi))
    Xbar.91 <- (1/n.quint)*with(data.unique, sum(X9*T1*psi))
    Xbar.101 <- (1/n.quint)*with(data.unique, sum(X10*T1*psi))
    Xbar.111 <- (1/n.quint)*with(data.unique, sum(X11*T1*psi))
    Xbar.121 <- (1/n.quint)*with(data.unique, sum(X12*T1*psi))
    Xbar.131 <- (1/n.quint)*with(data.unique, sum(X13*T1*psi))
    Xbar.141 <- (1/n.quint)*with(data.unique, sum(X14*T1*psi))
    Xbar.151 <- (1/n.quint)*with(data.unique, sum(X15*T1*psi))
    Xbar.161 <- (1/n.quint)*with(data.unique, sum(X16*T1*psi))
    Xbar.171 <- (1/n.quint)*with(data.unique, sum(X17*T1*psi))
    Xbar.181 <- (1/n.quint)*with(data.unique, sum(X18*T1*psi))
    Xbar.191 <- (1/n.quint)*with(data.unique, sum(X19*T1*psi))
    Xbar.201 <- (1/n.quint)*with(data.unique, sum(X20*T1*psi))
    
    Xbar.12 <- (1/n.quint)*with(data.unique, sum(X1*T2*psi))
    Xbar.22 <- (1/n.quint)*with(data.unique, sum(X2*T2*psi))
    Xbar.32 <- (1/n.quint)*with(data.unique, sum(X3*T2*psi))
    Xbar.42 <- (1/n.quint)*with(data.unique, sum(X4*T2*psi))
    Xbar.52 <- (1/n.quint)*with(data.unique, sum(X5*T2*psi))
    Xbar.62 <- (1/n.quint)*with(data.unique, sum(X6*T2*psi))
    Xbar.72 <- (1/n.quint)*with(data.unique, sum(X7*T2*psi))
    Xbar.82 <- (1/n.quint)*with(data.unique, sum(X8*T2*psi))
    Xbar.92 <- (1/n.quint)*with(data.unique, sum(X9*T2*psi))
    Xbar.102 <- (1/n.quint)*with(data.unique, sum(X10*T2*psi))
    Xbar.112 <- (1/n.quint)*with(data.unique, sum(X11*T2*psi))
    Xbar.122 <- (1/n.quint)*with(data.unique, sum(X12*T2*psi))
    Xbar.132 <- (1/n.quint)*with(data.unique, sum(X13*T2*psi))
    Xbar.142 <- (1/n.quint)*with(data.unique, sum(X14*T2*psi))
    Xbar.152 <- (1/n.quint)*with(data.unique, sum(X15*T2*psi))
    Xbar.162 <- (1/n.quint)*with(data.unique, sum(X16*T2*psi))
    Xbar.172 <- (1/n.quint)*with(data.unique, sum(X17*T2*psi))
    Xbar.182 <- (1/n.quint)*with(data.unique, sum(X18*T2*psi))
    Xbar.192 <- (1/n.quint)*with(data.unique, sum(X19*T2*psi))
    Xbar.202 <- (1/n.quint)*with(data.unique, sum(X20*T2*psi))
    
    Xbar.13 <- (1/n.quint)*with(data.unique, sum(X1*T3*psi))
    Xbar.23 <- (1/n.quint)*with(data.unique, sum(X2*T3*psi))
    Xbar.33 <- (1/n.quint)*with(data.unique, sum(X3*T3*psi))
    Xbar.43 <- (1/n.quint)*with(data.unique, sum(X4*T3*psi))
    Xbar.53 <- (1/n.quint)*with(data.unique, sum(X5*T3*psi))
    Xbar.63 <- (1/n.quint)*with(data.unique, sum(X6*T3*psi))
    Xbar.73 <- (1/n.quint)*with(data.unique, sum(X7*T3*psi))
    Xbar.83 <- (1/n.quint)*with(data.unique, sum(X8*T3*psi))
    Xbar.93 <- (1/n.quint)*with(data.unique, sum(X9*T3*psi))
    Xbar.103 <- (1/n.quint)*with(data.unique, sum(X10*T3*psi))
    Xbar.113 <- (1/n.quint)*with(data.unique, sum(X11*T3*psi))
    Xbar.123 <- (1/n.quint)*with(data.unique, sum(X12*T3*psi))
    Xbar.133 <- (1/n.quint)*with(data.unique, sum(X13*T3*psi))
    Xbar.143 <- (1/n.quint)*with(data.unique, sum(X14*T3*psi))
    Xbar.153 <- (1/n.quint)*with(data.unique, sum(X15*T3*psi))
    Xbar.163 <- (1/n.quint)*with(data.unique, sum(X16*T3*psi))
    Xbar.173 <- (1/n.quint)*with(data.unique, sum(X17*T3*psi))
    Xbar.183 <- (1/n.quint)*with(data.unique, sum(X18*T3*psi))
    Xbar.193 <- (1/n.quint)*with(data.unique, sum(X19*T3*psi))
    Xbar.203 <- (1/n.quint)*with(data.unique, sum(X20*T3*psi))
    
    Xbar.14 <- (1/n.quint)*with(data.unique, sum(X1*T4*psi))
    Xbar.24 <- (1/n.quint)*with(data.unique, sum(X2*T4*psi))
    Xbar.34 <- (1/n.quint)*with(data.unique, sum(X3*T4*psi))
    Xbar.44 <- (1/n.quint)*with(data.unique, sum(X4*T4*psi))
    Xbar.54 <- (1/n.quint)*with(data.unique, sum(X5*T4*psi))
    Xbar.64 <- (1/n.quint)*with(data.unique, sum(X6*T4*psi))
    Xbar.74 <- (1/n.quint)*with(data.unique, sum(X7*T4*psi))
    Xbar.84 <- (1/n.quint)*with(data.unique, sum(X8*T4*psi))
    Xbar.94 <- (1/n.quint)*with(data.unique, sum(X9*T4*psi))
    Xbar.104 <- (1/n.quint)*with(data.unique, sum(X10*T4*psi))
    Xbar.114 <- (1/n.quint)*with(data.unique, sum(X11*T4*psi))
    Xbar.124 <- (1/n.quint)*with(data.unique, sum(X12*T4*psi))
    Xbar.134 <- (1/n.quint)*with(data.unique, sum(X13*T4*psi))
    Xbar.144 <- (1/n.quint)*with(data.unique, sum(X14*T4*psi))
    Xbar.154 <- (1/n.quint)*with(data.unique, sum(X15*T4*psi))
    Xbar.164 <- (1/n.quint)*with(data.unique, sum(X16*T4*psi))
    Xbar.174 <- (1/n.quint)*with(data.unique, sum(X17*T4*psi))
    Xbar.184 <- (1/n.quint)*with(data.unique, sum(X18*T4*psi))
    Xbar.194 <- (1/n.quint)*with(data.unique, sum(X19*T4*psi))
    Xbar.204 <- (1/n.quint)*with(data.unique, sum(X20*T4*psi))
    
    Xbar.15 <- (1/n.quint)*with(data.unique, sum(X1*T5*psi))
    Xbar.25 <- (1/n.quint)*with(data.unique, sum(X2*T5*psi))
    Xbar.35 <- (1/n.quint)*with(data.unique, sum(X3*T5*psi))
    Xbar.45 <- (1/n.quint)*with(data.unique, sum(X4*T5*psi))
    Xbar.55 <- (1/n.quint)*with(data.unique, sum(X5*T5*psi))
    Xbar.65 <- (1/n.quint)*with(data.unique, sum(X6*T5*psi))
    Xbar.75 <- (1/n.quint)*with(data.unique, sum(X7*T5*psi))
    Xbar.85 <- (1/n.quint)*with(data.unique, sum(X8*T5*psi))
    Xbar.95 <- (1/n.quint)*with(data.unique, sum(X9*T5*psi))
    Xbar.105 <- (1/n.quint)*with(data.unique, sum(X10*T5*psi))
    Xbar.115 <- (1/n.quint)*with(data.unique, sum(X11*T5*psi))
    Xbar.125 <- (1/n.quint)*with(data.unique, sum(X12*T5*psi))
    Xbar.135 <- (1/n.quint)*with(data.unique, sum(X13*T5*psi))
    Xbar.145 <- (1/n.quint)*with(data.unique, sum(X14*T5*psi))
    Xbar.155 <- (1/n.quint)*with(data.unique, sum(X15*T5*psi))
    Xbar.165 <- (1/n.quint)*with(data.unique, sum(X16*T5*psi))
    Xbar.175 <- (1/n.quint)*with(data.unique, sum(X17*T5*psi))
    Xbar.185 <- (1/n.quint)*with(data.unique, sum(X18*T5*psi))
    Xbar.195 <- (1/n.quint)*with(data.unique, sum(X19*T5*psi))
    Xbar.205 <- (1/n.quint)*with(data.unique, sum(X20*T5*psi))
    
    sb112 <- with(data.unique, (Xbar.11 - Xbar.12)/sd(data$X1[data$treat == "Treatment 1"]))
    sb212 <- with(data.unique, (Xbar.21 - Xbar.22)/sd(data$X2[data$treat == "Treatment 1"]))
    sb312 <- with(data.unique, (Xbar.31 - Xbar.32)/sd(data$X3[data$treat == "Treatment 1"]))
    sb412 <- with(data.unique, (Xbar.41 - Xbar.42)/sd(data$X4[data$treat == "Treatment 1"]))
    sb512 <- with(data.unique, (Xbar.51 - Xbar.52)/sd(data$X5[data$treat == "Treatment 1"]))
    sb612 <- with(data.unique, (Xbar.61 - Xbar.62)/sd(data$X6[data$treat == "Treatment 1"]))
    sb712 <- with(data.unique, (Xbar.71 - Xbar.72)/sd(data$X7[data$treat == "Treatment 1"]))
    sb812 <- with(data.unique, (Xbar.81 - Xbar.82)/sd(data$X8[data$treat == "Treatment 1"]))
    sb912 <- with(data.unique, (Xbar.91 - Xbar.92)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1012 <- with(data.unique, (Xbar.101 - Xbar.102)/sd(data$X10[data$treat == "Treatment 1"]))
    sb1112 <- with(data.unique, (Xbar.111 - Xbar.112)/sd(data$X11[data$treat == "Treatment 1"]))
    sb1212 <- with(data.unique, (Xbar.121 - Xbar.122)/sd(data$X12[data$treat == "Treatment 1"]))
    sb1312 <- with(data.unique, (Xbar.131 - Xbar.132)/sd(data$X13[data$treat == "Treatment 1"]))
    sb1412 <- with(data.unique, (Xbar.141 - Xbar.142)/sd(data$X14[data$treat == "Treatment 1"]))
    sb1512 <- with(data.unique, (Xbar.151 - Xbar.152)/sd(data$X15[data$treat == "Treatment 1"]))
    sb1612 <- with(data.unique, (Xbar.161 - Xbar.162)/sd(data$X16[data$treat == "Treatment 1"]))
    sb1712 <- with(data.unique, (Xbar.171 - Xbar.172)/sd(data$X17[data$treat == "Treatment 1"]))
    sb1812 <- with(data.unique, (Xbar.181 - Xbar.182)/sd(data$X18[data$treat == "Treatment 1"]))
    sb1912 <- with(data.unique, (Xbar.191 - Xbar.192)/sd(data$X19[data$treat == "Treatment 1"]))
    sb2012 <- with(data.unique, (Xbar.201 - Xbar.202)/sd(data$X20[data$treat == "Treatment 1"]))
    
    sb113 <- with(data.unique, (Xbar.11 - Xbar.13)/sd(data$X1[data$treat == "Treatment 1"]))
    sb213 <- with(data.unique, (Xbar.21 - Xbar.23)/sd(data$X2[data$treat == "Treatment 1"]))
    sb313 <- with(data.unique, (Xbar.31 - Xbar.33)/sd(data$X3[data$treat == "Treatment 1"]))
    sb413 <- with(data.unique, (Xbar.41 - Xbar.43)/sd(data$X4[data$treat == "Treatment 1"]))
    sb513 <- with(data.unique, (Xbar.51 - Xbar.53)/sd(data$X5[data$treat == "Treatment 1"]))
    sb613 <- with(data.unique, (Xbar.61 - Xbar.63)/sd(data$X6[data$treat == "Treatment 1"]))
    sb713 <- with(data.unique, (Xbar.71 - Xbar.73)/sd(data$X7[data$treat == "Treatment 1"]))
    sb813 <- with(data.unique, (Xbar.81 - Xbar.83)/sd(data$X8[data$treat == "Treatment 1"]))
    sb913 <- with(data.unique, (Xbar.91 - Xbar.93)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1013 <- with(data.unique, (Xbar.101 - Xbar.103)/sd(data$X10[data$treat == "Treatment 1"]))
    sb1113 <- with(data.unique, (Xbar.111 - Xbar.113)/sd(data$X11[data$treat == "Treatment 1"]))
    sb1213 <- with(data.unique, (Xbar.121 - Xbar.123)/sd(data$X12[data$treat == "Treatment 1"]))
    sb1313 <- with(data.unique, (Xbar.131 - Xbar.133)/sd(data$X13[data$treat == "Treatment 1"]))
    sb1413 <- with(data.unique, (Xbar.141 - Xbar.143)/sd(data$X14[data$treat == "Treatment 1"]))
    sb1513 <- with(data.unique, (Xbar.151 - Xbar.153)/sd(data$X15[data$treat == "Treatment 1"]))
    sb1613 <- with(data.unique, (Xbar.161 - Xbar.163)/sd(data$X16[data$treat == "Treatment 1"]))
    sb1713 <- with(data.unique, (Xbar.171 - Xbar.173)/sd(data$X17[data$treat == "Treatment 1"]))
    sb1813 <- with(data.unique, (Xbar.181 - Xbar.183)/sd(data$X18[data$treat == "Treatment 1"]))
    sb1913 <- with(data.unique, (Xbar.191 - Xbar.193)/sd(data$X19[data$treat == "Treatment 1"]))
    sb2013 <- with(data.unique, (Xbar.201 - Xbar.203)/sd(data$X20[data$treat == "Treatment 1"]))
    
    sb114 <- with(data.unique, (Xbar.11 - Xbar.14)/sd(data$X1[data$treat == "Treatment 1"]))
    sb214 <- with(data.unique, (Xbar.21 - Xbar.24)/sd(data$X2[data$treat == "Treatment 1"]))
    sb314 <- with(data.unique, (Xbar.31 - Xbar.34)/sd(data$X3[data$treat == "Treatment 1"]))
    sb414 <- with(data.unique, (Xbar.41 - Xbar.44)/sd(data$X4[data$treat == "Treatment 1"]))
    sb514 <- with(data.unique, (Xbar.51 - Xbar.54)/sd(data$X5[data$treat == "Treatment 1"]))
    sb614 <- with(data.unique, (Xbar.61 - Xbar.64)/sd(data$X6[data$treat == "Treatment 1"]))
    sb714 <- with(data.unique, (Xbar.71 - Xbar.74)/sd(data$X7[data$treat == "Treatment 1"]))
    sb814 <- with(data.unique, (Xbar.81 - Xbar.84)/sd(data$X8[data$treat == "Treatment 1"]))
    sb914 <- with(data.unique, (Xbar.91 - Xbar.94)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1014 <- with(data.unique, (Xbar.101 - Xbar.104)/sd(data$X10[data$treat == "Treatment 1"]))
    sb1114 <- with(data.unique, (Xbar.111 - Xbar.114)/sd(data$X11[data$treat == "Treatment 1"]))
    sb1214 <- with(data.unique, (Xbar.121 - Xbar.124)/sd(data$X12[data$treat == "Treatment 1"]))
    sb1314 <- with(data.unique, (Xbar.131 - Xbar.134)/sd(data$X13[data$treat == "Treatment 1"]))
    sb1414 <- with(data.unique, (Xbar.141 - Xbar.144)/sd(data$X14[data$treat == "Treatment 1"]))
    sb1514 <- with(data.unique, (Xbar.151 - Xbar.154)/sd(data$X15[data$treat == "Treatment 1"]))
    sb1614 <- with(data.unique, (Xbar.161 - Xbar.164)/sd(data$X16[data$treat == "Treatment 1"]))
    sb1714 <- with(data.unique, (Xbar.171 - Xbar.174)/sd(data$X17[data$treat == "Treatment 1"]))
    sb1814 <- with(data.unique, (Xbar.181 - Xbar.184)/sd(data$X18[data$treat == "Treatment 1"]))
    sb1914 <- with(data.unique, (Xbar.191 - Xbar.194)/sd(data$X19[data$treat == "Treatment 1"]))
    sb2014 <- with(data.unique, (Xbar.201 - Xbar.204)/sd(data$X20[data$treat == "Treatment 1"]))
    
    sb115 <- with(data.unique, (Xbar.11 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb215 <- with(data.unique, (Xbar.21 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb315 <- with(data.unique, (Xbar.31 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb415 <- with(data.unique, (Xbar.41 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb515 <- with(data.unique, (Xbar.51 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    sb615 <- with(data.unique, (Xbar.61 - Xbar.65)/sd(data$X6[data$treat == "Treatment 1"]))
    sb715 <- with(data.unique, (Xbar.71 - Xbar.75)/sd(data$X7[data$treat == "Treatment 1"]))
    sb815 <- with(data.unique, (Xbar.81 - Xbar.85)/sd(data$X8[data$treat == "Treatment 1"]))
    sb915 <- with(data.unique, (Xbar.91 - Xbar.95)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1015 <- with(data.unique, (Xbar.101 - Xbar.105)/sd(data$X10[data$treat == "Treatment 1"]))
    sb1115 <- with(data.unique, (Xbar.111 - Xbar.115)/sd(data$X11[data$treat == "Treatment 1"]))
    sb1215 <- with(data.unique, (Xbar.121 - Xbar.125)/sd(data$X12[data$treat == "Treatment 1"]))
    sb1315 <- with(data.unique, (Xbar.131 - Xbar.135)/sd(data$X13[data$treat == "Treatment 1"]))
    sb1415 <- with(data.unique, (Xbar.141 - Xbar.145)/sd(data$X14[data$treat == "Treatment 1"]))
    sb1515 <- with(data.unique, (Xbar.151 - Xbar.155)/sd(data$X15[data$treat == "Treatment 1"]))
    sb1615 <- with(data.unique, (Xbar.161 - Xbar.165)/sd(data$X16[data$treat == "Treatment 1"]))
    sb1715 <- with(data.unique, (Xbar.171 - Xbar.175)/sd(data$X17[data$treat == "Treatment 1"]))
    sb1815 <- with(data.unique, (Xbar.181 - Xbar.185)/sd(data$X18[data$treat == "Treatment 1"]))
    sb1915 <- with(data.unique, (Xbar.191 - Xbar.195)/sd(data$X19[data$treat == "Treatment 1"]))
    sb2015 <- with(data.unique, (Xbar.201 - Xbar.205)/sd(data$X20[data$treat == "Treatment 1"]))
    
    sb123 <- with(data.unique, (Xbar.12 - Xbar.13)/sd(data$X1[data$treat == "Treatment 1"]))
    sb223 <- with(data.unique, (Xbar.22 - Xbar.23)/sd(data$X2[data$treat == "Treatment 1"]))
    sb323 <- with(data.unique, (Xbar.32 - Xbar.33)/sd(data$X3[data$treat == "Treatment 1"]))
    sb423 <- with(data.unique, (Xbar.42 - Xbar.43)/sd(data$X4[data$treat == "Treatment 1"]))
    sb523 <- with(data.unique, (Xbar.52 - Xbar.53)/sd(data$X5[data$treat == "Treatment 1"]))
    sb623 <- with(data.unique, (Xbar.62 - Xbar.63)/sd(data$X6[data$treat == "Treatment 1"]))
    sb723 <- with(data.unique, (Xbar.72 - Xbar.73)/sd(data$X7[data$treat == "Treatment 1"]))
    sb823 <- with(data.unique, (Xbar.82 - Xbar.83)/sd(data$X8[data$treat == "Treatment 1"]))
    sb923 <- with(data.unique, (Xbar.92 - Xbar.93)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1023 <- with(data.unique, (Xbar.102 - Xbar.103)/sd(data$X10[data$treat == "Treatment 1"]))
    sb1123 <- with(data.unique, (Xbar.112 - Xbar.113)/sd(data$X11[data$treat == "Treatment 1"]))
    sb1223 <- with(data.unique, (Xbar.122 - Xbar.123)/sd(data$X12[data$treat == "Treatment 1"]))
    sb1323 <- with(data.unique, (Xbar.132 - Xbar.133)/sd(data$X13[data$treat == "Treatment 1"]))
    sb1423 <- with(data.unique, (Xbar.142 - Xbar.143)/sd(data$X14[data$treat == "Treatment 1"]))
    sb1523 <- with(data.unique, (Xbar.152 - Xbar.153)/sd(data$X15[data$treat == "Treatment 1"]))
    sb1623 <- with(data.unique, (Xbar.162 - Xbar.163)/sd(data$X16[data$treat == "Treatment 1"]))
    sb1723 <- with(data.unique, (Xbar.172 - Xbar.173)/sd(data$X17[data$treat == "Treatment 1"]))
    sb1823 <- with(data.unique, (Xbar.182 - Xbar.183)/sd(data$X18[data$treat == "Treatment 1"]))
    sb1923 <- with(data.unique, (Xbar.192 - Xbar.193)/sd(data$X19[data$treat == "Treatment 1"]))
    sb2023 <- with(data.unique, (Xbar.202 - Xbar.203)/sd(data$X20[data$treat == "Treatment 1"]))
    
    sb124 <- with(data.unique, (Xbar.12 - Xbar.14)/sd(data$X1[data$treat == "Treatment 1"]))
    sb224 <- with(data.unique, (Xbar.22 - Xbar.24)/sd(data$X2[data$treat == "Treatment 1"]))
    sb324 <- with(data.unique, (Xbar.32 - Xbar.34)/sd(data$X3[data$treat == "Treatment 1"]))
    sb424 <- with(data.unique, (Xbar.42 - Xbar.44)/sd(data$X4[data$treat == "Treatment 1"]))
    sb524 <- with(data.unique, (Xbar.52 - Xbar.54)/sd(data$X5[data$treat == "Treatment 1"]))
    sb624 <- with(data.unique, (Xbar.62 - Xbar.64)/sd(data$X6[data$treat == "Treatment 1"]))
    sb724 <- with(data.unique, (Xbar.72 - Xbar.74)/sd(data$X7[data$treat == "Treatment 1"]))
    sb824 <- with(data.unique, (Xbar.82 - Xbar.84)/sd(data$X8[data$treat == "Treatment 1"]))
    sb924 <- with(data.unique, (Xbar.92 - Xbar.94)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1024 <- with(data.unique, (Xbar.102 - Xbar.104)/sd(data$X10[data$treat == "Treatment 1"]))
    sb1124 <- with(data.unique, (Xbar.112 - Xbar.114)/sd(data$X11[data$treat == "Treatment 1"]))
    sb1224 <- with(data.unique, (Xbar.122 - Xbar.124)/sd(data$X12[data$treat == "Treatment 1"]))
    sb1324 <- with(data.unique, (Xbar.132 - Xbar.134)/sd(data$X13[data$treat == "Treatment 1"]))
    sb1424 <- with(data.unique, (Xbar.142 - Xbar.144)/sd(data$X14[data$treat == "Treatment 1"]))
    sb1524 <- with(data.unique, (Xbar.152 - Xbar.154)/sd(data$X15[data$treat == "Treatment 1"]))
    sb1624 <- with(data.unique, (Xbar.162 - Xbar.164)/sd(data$X16[data$treat == "Treatment 1"]))
    sb1724 <- with(data.unique, (Xbar.172 - Xbar.174)/sd(data$X17[data$treat == "Treatment 1"]))
    sb1824 <- with(data.unique, (Xbar.182 - Xbar.184)/sd(data$X18[data$treat == "Treatment 1"]))
    sb1924 <- with(data.unique, (Xbar.192 - Xbar.194)/sd(data$X19[data$treat == "Treatment 1"]))
    sb2024 <- with(data.unique, (Xbar.202 - Xbar.204)/sd(data$X20[data$treat == "Treatment 1"]))
    
    sb125 <- with(data.unique, (Xbar.12 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb225 <- with(data.unique, (Xbar.22 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb325 <- with(data.unique, (Xbar.32 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb425 <- with(data.unique, (Xbar.42 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb525 <- with(data.unique, (Xbar.52 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    sb625 <- with(data.unique, (Xbar.62 - Xbar.65)/sd(data$X6[data$treat == "Treatment 1"]))
    sb725 <- with(data.unique, (Xbar.72 - Xbar.75)/sd(data$X7[data$treat == "Treatment 1"]))
    sb825 <- with(data.unique, (Xbar.82 - Xbar.85)/sd(data$X8[data$treat == "Treatment 1"]))
    sb925 <- with(data.unique, (Xbar.92 - Xbar.95)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1025 <- with(data.unique, (Xbar.102 - Xbar.105)/sd(data$X10[data$treat == "Treatment 1"]))
    sb1125 <- with(data.unique, (Xbar.112 - Xbar.115)/sd(data$X11[data$treat == "Treatment 1"]))
    sb1225 <- with(data.unique, (Xbar.122 - Xbar.125)/sd(data$X12[data$treat == "Treatment 1"]))
    sb1325 <- with(data.unique, (Xbar.132 - Xbar.135)/sd(data$X13[data$treat == "Treatment 1"]))
    sb1425 <- with(data.unique, (Xbar.142 - Xbar.145)/sd(data$X14[data$treat == "Treatment 1"]))
    sb1525 <- with(data.unique, (Xbar.152 - Xbar.155)/sd(data$X15[data$treat == "Treatment 1"]))
    sb1625 <- with(data.unique, (Xbar.162 - Xbar.165)/sd(data$X16[data$treat == "Treatment 1"]))
    sb1725 <- with(data.unique, (Xbar.172 - Xbar.175)/sd(data$X17[data$treat == "Treatment 1"]))
    sb1825 <- with(data.unique, (Xbar.182 - Xbar.185)/sd(data$X18[data$treat == "Treatment 1"]))
    sb1925 <- with(data.unique, (Xbar.192 - Xbar.195)/sd(data$X19[data$treat == "Treatment 1"]))
    sb2025 <- with(data.unique, (Xbar.202 - Xbar.205)/sd(data$X20[data$treat == "Treatment 1"]))
    
    sb134 <- with(data.unique, (Xbar.13 - Xbar.14)/sd(data$X1[data$treat == "Treatment 1"]))
    sb234 <- with(data.unique, (Xbar.23 - Xbar.24)/sd(data$X2[data$treat == "Treatment 1"]))
    sb334 <- with(data.unique, (Xbar.33 - Xbar.34)/sd(data$X3[data$treat == "Treatment 1"]))
    sb434 <- with(data.unique, (Xbar.43 - Xbar.44)/sd(data$X4[data$treat == "Treatment 1"]))
    sb534 <- with(data.unique, (Xbar.53 - Xbar.54)/sd(data$X5[data$treat == "Treatment 1"]))
    sb634 <- with(data.unique, (Xbar.63 - Xbar.64)/sd(data$X6[data$treat == "Treatment 1"]))
    sb734 <- with(data.unique, (Xbar.73 - Xbar.74)/sd(data$X7[data$treat == "Treatment 1"]))
    sb834 <- with(data.unique, (Xbar.83 - Xbar.84)/sd(data$X8[data$treat == "Treatment 1"]))
    sb934 <- with(data.unique, (Xbar.93 - Xbar.94)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1034 <- with(data.unique, (Xbar.103 - Xbar.104)/sd(data$X10[data$treat == "Treatment 1"]))
    sb1134 <- with(data.unique, (Xbar.113 - Xbar.114)/sd(data$X11[data$treat == "Treatment 1"]))
    sb1234 <- with(data.unique, (Xbar.123 - Xbar.124)/sd(data$X12[data$treat == "Treatment 1"]))
    sb1334 <- with(data.unique, (Xbar.133 - Xbar.134)/sd(data$X13[data$treat == "Treatment 1"]))
    sb1434 <- with(data.unique, (Xbar.143 - Xbar.144)/sd(data$X14[data$treat == "Treatment 1"]))
    sb1534 <- with(data.unique, (Xbar.153 - Xbar.154)/sd(data$X15[data$treat == "Treatment 1"]))
    sb1634 <- with(data.unique, (Xbar.163 - Xbar.164)/sd(data$X16[data$treat == "Treatment 1"]))
    sb1734 <- with(data.unique, (Xbar.173 - Xbar.174)/sd(data$X17[data$treat == "Treatment 1"]))
    sb1834 <- with(data.unique, (Xbar.183 - Xbar.184)/sd(data$X18[data$treat == "Treatment 1"]))
    sb1934 <- with(data.unique, (Xbar.193 - Xbar.194)/sd(data$X19[data$treat == "Treatment 1"]))
    sb2034 <- with(data.unique, (Xbar.203 - Xbar.204)/sd(data$X20[data$treat == "Treatment 1"]))
    
    sb135 <- with(data.unique, (Xbar.13 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb235 <- with(data.unique, (Xbar.23 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb335 <- with(data.unique, (Xbar.33 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb435 <- with(data.unique, (Xbar.43 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb535 <- with(data.unique, (Xbar.53 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    sb635 <- with(data.unique, (Xbar.63 - Xbar.65)/sd(data$X6[data$treat == "Treatment 1"]))
    sb735 <- with(data.unique, (Xbar.73 - Xbar.75)/sd(data$X7[data$treat == "Treatment 1"]))
    sb835 <- with(data.unique, (Xbar.83 - Xbar.85)/sd(data$X8[data$treat == "Treatment 1"]))
    sb935 <- with(data.unique, (Xbar.93 - Xbar.95)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1035 <- with(data.unique, (Xbar.103 - Xbar.105)/sd(data$X10[data$treat == "Treatment 1"]))
    sb1135 <- with(data.unique, (Xbar.113 - Xbar.115)/sd(data$X11[data$treat == "Treatment 1"]))
    sb1235 <- with(data.unique, (Xbar.123 - Xbar.125)/sd(data$X12[data$treat == "Treatment 1"]))
    sb1335 <- with(data.unique, (Xbar.133 - Xbar.135)/sd(data$X13[data$treat == "Treatment 1"]))
    sb1435 <- with(data.unique, (Xbar.143 - Xbar.145)/sd(data$X14[data$treat == "Treatment 1"]))
    sb1535 <- with(data.unique, (Xbar.153 - Xbar.155)/sd(data$X15[data$treat == "Treatment 1"]))
    sb1635 <- with(data.unique, (Xbar.163 - Xbar.165)/sd(data$X16[data$treat == "Treatment 1"]))
    sb1735 <- with(data.unique, (Xbar.173 - Xbar.175)/sd(data$X17[data$treat == "Treatment 1"]))
    sb1835 <- with(data.unique, (Xbar.183 - Xbar.185)/sd(data$X18[data$treat == "Treatment 1"]))
    sb1935 <- with(data.unique, (Xbar.193 - Xbar.195)/sd(data$X19[data$treat == "Treatment 1"]))
    sb2035 <- with(data.unique, (Xbar.203 - Xbar.205)/sd(data$X20[data$treat == "Treatment 1"]))
    
    sb145 <- with(data.unique, (Xbar.14 - Xbar.15)/sd(data$X1[data$treat == "Treatment 1"]))
    sb245 <- with(data.unique, (Xbar.24 - Xbar.25)/sd(data$X2[data$treat == "Treatment 1"]))
    sb345 <- with(data.unique, (Xbar.34 - Xbar.35)/sd(data$X3[data$treat == "Treatment 1"]))
    sb445 <- with(data.unique, (Xbar.44 - Xbar.45)/sd(data$X4[data$treat == "Treatment 1"]))
    sb545 <- with(data.unique, (Xbar.54 - Xbar.55)/sd(data$X5[data$treat == "Treatment 1"]))
    sb645 <- with(data.unique, (Xbar.64 - Xbar.65)/sd(data$X6[data$treat == "Treatment 1"]))
    sb745 <- with(data.unique, (Xbar.74 - Xbar.75)/sd(data$X7[data$treat == "Treatment 1"]))
    sb845 <- with(data.unique, (Xbar.84 - Xbar.85)/sd(data$X8[data$treat == "Treatment 1"]))
    sb945 <- with(data.unique, (Xbar.94 - Xbar.95)/sd(data$X9[data$treat == "Treatment 1"]))
    sb1045 <- with(data.unique, (Xbar.104 - Xbar.105)/sd(data$X10[data$treat == "Treatment 1"]))
    sb1145 <- with(data.unique, (Xbar.114 - Xbar.115)/sd(data$X11[data$treat == "Treatment 1"]))
    sb1245 <- with(data.unique, (Xbar.124 - Xbar.125)/sd(data$X12[data$treat == "Treatment 1"]))
    sb1345 <- with(data.unique, (Xbar.134 - Xbar.135)/sd(data$X13[data$treat == "Treatment 1"]))
    sb1445 <- with(data.unique, (Xbar.144 - Xbar.145)/sd(data$X14[data$treat == "Treatment 1"]))
    sb1545 <- with(data.unique, (Xbar.154 - Xbar.155)/sd(data$X15[data$treat == "Treatment 1"]))
    sb1645 <- with(data.unique, (Xbar.164 - Xbar.165)/sd(data$X16[data$treat == "Treatment 1"]))
    sb1745 <- with(data.unique, (Xbar.174 - Xbar.175)/sd(data$X17[data$treat == "Treatment 1"]))
    sb1845 <- with(data.unique, (Xbar.184 - Xbar.185)/sd(data$X18[data$treat == "Treatment 1"]))
    sb1945 <- with(data.unique, (Xbar.194 - Xbar.195)/sd(data$X19[data$treat == "Treatment 1"]))
    sb2045 <- with(data.unique, (Xbar.204 - Xbar.205)/sd(data$X20[data$treat == "Treatment 1"]))
    
    max2sb1 <- max(abs(sb112), abs(sb113), abs(sb114), abs(sb115),
                   abs(sb123), abs(sb124), abs(sb125), 
                   abs(sb134), abs(sb135), abs(sb145))
    max2sb2 <- max(abs(sb212), abs(sb213), abs(sb214), abs(sb215),
                   abs(sb223), abs(sb224), abs(sb225), 
                   abs(sb234), abs(sb235), abs(sb245))
    max2sb3 <- max(abs(sb312), abs(sb313), abs(sb314), abs(sb315),
                   abs(sb323), abs(sb324), abs(sb325), 
                   abs(sb334), abs(sb335), abs(sb345))
    max2sb4 <- max(abs(sb412), abs(sb413), abs(sb414), abs(sb415),
                   abs(sb423), abs(sb424), abs(sb425), 
                   abs(sb434), abs(sb435), abs(sb445))
    max2sb5 <- max(abs(sb512), abs(sb513), abs(sb514), abs(sb515),
                   abs(sb523), abs(sb524), abs(sb525), 
                   abs(sb534), abs(sb535), abs(sb545))
    max2sb6 <- max(abs(sb612), abs(sb613), abs(sb614), abs(sb615),
                   abs(sb623), abs(sb624), abs(sb625), 
                   abs(sb634), abs(sb635), abs(sb645))
    max2sb7 <- max(abs(sb712), abs(sb713), abs(sb714), abs(sb715),
                   abs(sb723), abs(sb724), abs(sb725), 
                   abs(sb734), abs(sb735), abs(sb745))
    max2sb8 <- max(abs(sb812), abs(sb813), abs(sb814), abs(sb815),
                   abs(sb823), abs(sb824), abs(sb825), 
                   abs(sb834), abs(sb835), abs(sb845))
    max2sb9 <- max(abs(sb912), abs(sb913), abs(sb914), abs(sb915),
                   abs(sb923), abs(sb924), abs(sb925), 
                   abs(sb934), abs(sb935), abs(sb945))
    max2sb10 <- max(abs(sb1012), abs(sb1013), abs(sb1014), abs(sb1015),
                    abs(sb1023), abs(sb1024), abs(sb1025), 
                    abs(sb1034), abs(sb1035), abs(sb1045))
    max2sb11 <- max(abs(sb1112), abs(sb1113), abs(sb1114), abs(sb1115),
                    abs(sb1123), abs(sb1124), abs(sb1125), 
                    abs(sb1134), abs(sb1135), abs(sb1145))
    max2sb12 <- max(abs(sb1212), abs(sb1213), abs(sb1214), abs(sb1215),
                    abs(sb1223), abs(sb1224), abs(sb1225), 
                    abs(sb1234), abs(sb1235), abs(sb1245))
    max2sb13 <- max(abs(sb1312), abs(sb1313), abs(sb1314), abs(sb1315),
                    abs(sb1323), abs(sb1324), abs(sb1325), 
                    abs(sb1334), abs(sb1335), abs(sb1345))
    max2sb14 <- max(abs(sb1412), abs(sb1413), abs(sb1414), abs(sb1415),
                    abs(sb1423), abs(sb1424), abs(sb1425), 
                    abs(sb1434), abs(sb1435), abs(sb1445))
    max2sb15 <- max(abs(sb1512), abs(sb1513), abs(sb1514), abs(sb1515),
                    abs(sb1523), abs(sb1524), abs(sb1525), 
                    abs(sb1534), abs(sb1535), abs(sb1545))
    max2sb16 <- max(abs(sb1612), abs(sb1613), abs(sb1614), abs(sb1615),
                    abs(sb1623), abs(sb1624), abs(sb1625), 
                    abs(sb1634), abs(sb1635), abs(sb1645))
    max2sb17 <- max(abs(sb1712), abs(sb1713), abs(sb1714), abs(sb1715),
                    abs(sb1723), abs(sb1724), abs(sb1725), 
                    abs(sb1734), abs(sb1735), abs(sb1745))
    max2sb18 <- max(abs(sb1812), abs(sb1813), abs(sb1814), abs(sb1815),
                    abs(sb1823), abs(sb1824), abs(sb1825), 
                    abs(sb1834), abs(sb1835), abs(sb1845))
    max2sb19 <- max(abs(sb1912), abs(sb1913), abs(sb1914), abs(sb1915),
                    abs(sb1923), abs(sb1924), abs(sb1925), 
                    abs(sb1934), abs(sb1935), abs(sb1945))
    max2sb20 <- max(abs(sb2012), abs(sb2013), abs(sb2014), abs(sb2015),
                    abs(sb2023), abs(sb2024), abs(sb2025), 
                    abs(sb2034), abs(sb2035), abs(sb2045))
    
    max2sb <- mean(c(max2sb1, max2sb2, max2sb3, max2sb4, max2sb5, 
                     max2sb6, max2sb7, max2sb8, max2sb9, max2sb10, 
                     max2sb11, max2sb12, max2sb13, max2sb14, max2sb15,
                     max2sb16, max2sb17, max2sb18, max2sb19, max2sb20))
    maxmax2sb <- max(c(max2sb1, max2sb2, max2sb3, max2sb4, max2sb5, 
                       max2sb6, max2sb7, max2sb8, max2sb9, max2sb10, 
                       max2sb11, max2sb12, max2sb13, max2sb14, max2sb15,
                       max2sb16, max2sb17, max2sb18, max2sb19, max2sb20))
  }
  percent.matched <- n.quint/sum(data$T1)
  
  c(percent.matched, maxmax2sb)
}



Matchby.fuzzy <- function (Y = NULL, Tr, X, by, estimand = "ATT", M = 1, ties = FALSE, 
                           replace = TRUE, exact = NULL, caliper = NULL, AI = FALSE, 
                           Var.calc = 0, Weight = 1, Weight.matrix = NULL, distance.tolerance = 1e-05, 
                           tolerance = sqrt(.Machine$double.eps), print.level = 1, version = "Matchby", 
                           ...) 
{
  require(IRanges)
  Tr <- as.double(Tr)
  nobs <- length(Tr)
  orig.treated.nobs <- sum(Tr == 1)
  if (is.null(Y)) {
    Y <- rep(0, nobs)
  }
  else {
    Y <- as.double(Y)
    if (nobs != length(Y)) {
      stop("length(Tr) != length(Y)")
    }
  }
  X <- as.matrix(X)
  if (nobs != nrow(X)) {
    stop("length(Tr) != nrow(X)")
  }
  index.nobs <- 1:nobs
  t.index.nobs <- multisplit(index.nobs, by)
  nindx <- length(t.index.nobs)
  weights <- NULL
  index.treated <- NULL
  index.control <- NULL
  if (Var.calc < 0) {
    warning("User set 'Var.calc' to less than 0.  Resetting to the default which is 0.")
    Var.calc <- 0
  }
  if (Var.calc > 0) 
    AI <- TRUE
  if (AI & estimand != "ATT") {
    warning("This function can only calculate Abadie-Imbens SEs for 'ATT'.  For AI SEs and other estimands, please use 'Match()'. Setting 'AI=FALSE'")
    AI <- FALSE
    Var.calc <- 0
  }
  if (AI & ties != TRUE) {
    ties = TRUE
    warning("Abadie-Imbens SEs have been requested. Setting 'ties=TRUE'")
  }
  if (AI & replace != TRUE) {
    replace = TRUE
    warning("Abadie-Imbens SEs have been requested. Setting 'replace=TRUE'")
  }
  if (AI) 
    version <- "MatchbyAI"
  if (replace != FALSE & replace != TRUE) {
    warning("'replace' must be TRUE or FALSE.  Setting to TRUE")
    replace <- TRUE
  }
  if (replace == FALSE) {
    ties <- FALSE
  }
  if (ties != FALSE & ties != TRUE) {
    warning("'ties' must be TRUE or FALSE.  Setting to TRUE")
    ties <- TRUE
  }
  if (AI) {
    Kcount <- as.matrix(rep(0, nobs))
    KKcount <- as.matrix(rep(0, nobs))
    YCAUS <- matrix(0, nrow = nobs, ncol = 1)
    if (Var.calc > 0) 
      Sigs <- matrix(0, nrow = nobs, ncol = 1)
  }
  for (i in 1:nindx) {
    if (print.level > 0) 
      cat(i, "of", nindx, "groups\n")
    tmp.index.nobs <- t.index.nobs[[i]]
    f.Tr <- Tr[tmp.index.nobs]
    if (length(f.Tr) < 2) {
      next
    }
    if (var(f.Tr) == 0) {
      next
    }
    t1 <- Match(Y = Y[tmp.index.nobs], Tr = f.Tr, X = X[tmp.index.nobs, 
                                                        ], estimand = estimand, M = M, Var.calc = Var.calc, 
                exact = exact, caliper = caliper, replace = replace, 
                ties = ties, Weight = Weight, Weight.matrix = Weight.matrix, 
                tolerance = tolerance, distance.tolerance = distance.tolerance, 
                version = version, ...)
    if (is.na(t1[1])) {
      if (!is.null(exact) | !is.null(caliper)) {
        warning("no matches found in group ", i, " (probably because of the exact or caliper option) continuing")
      }
      else {
        warning("no matches found in group ", i, " continuing")
      }
      next
    }
    weights <- c(weights, t1$weights)
    index.treated <- c(index.treated, (tmp.index.nobs)[t1$index.treated])
    index.control <- c(index.control, (tmp.index.nobs)[t1$index.control])
    if (AI) {
      YCAUS[tmp.index.nobs] <- t1$YCAUS
      Kcount[tmp.index.nobs] <- t1$Kcount
      KKcount[tmp.index.nobs] <- t1$KKcount
      if (Var.calc > 0) {
        if (is.null(t1$Sigs)) {
          pfoo <- paste("Var.calc=", Var.calc, " cannot be calculated (group ", 
                        i, ").  Var.calc is probably set to a number larger than the possible number of matches in a subgroup", 
                        sep = "")
          stop(pfoo)
        }
        else {
          Sigs[tmp.index.nobs] <- t1$Sigs
        }
      }
    }
  }
  if (is.null(index.treated)) {
    warning("no valid matches were found")
    z <- NA
    class(z) <- "Matchby"
    return(z)
  }
  Yt <- Y[index.treated]
  Yc <- Y[index.control]
  sum.tw <- sum(weights)
  est <- sum(Yt * weights)/sum.tw - sum(Yc * weights)/sum.tw
  Tau <- Yt - Yc
  varest <- sum(((Tau - est)^2) * weights)/(sum.tw * sum.tw)
  se.standard <- sqrt(varest)
  ret <- list(est = est, se = NULL, se.standard = se.standard, 
              index.treated = index.treated, index.control = index.control, 
              weights = weights)
  ret$orig.nobs <- nobs
  ret$orig.wnobs <- nobs
  ret$nobs <- length(Yt)
  ret$wnobs <- sum.tw
  ret$orig.treated.nobs <- orig.treated.nobs
  ret$ndrops <- orig.treated.nobs - ret$wnobs
  ret$estimand <- estimand
  ret$version <- version
  class(ret) <- "Matchby"
  if (AI) {
    if (Var.calc == 0) {
      eps <- Tau - est
      eps.sq <- eps * eps
      Sigs <- 0.5 * matrix(1, nobs, 1) %*% (t(eps.sq) %*% 
                                              weights)/sum(weights)
    }
    SN <- orig.treated.nobs
    var.pop = sum((Sigs * ((1 - Tr) * Kcount * Kcount - (1 - 
                                                           Tr) * KKcount))/(SN * SN))
    dvar.pop <- sum(Tr * (YCAUS - est) * (YCAUS - est))/(SN * 
                                                           SN)
    var <- var.pop + dvar.pop
    se <- sqrt(var)
    ret$se <- se
  }
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  
  dist <- c()
  for (i in 1:length(ret$index.treated)){
    if (ncol(X) > 1){
      dist[i] <- euc.dist(X[ret$index.treated[i], ], X[ret$index.control[i], ])
    }
    if (ncol(X) == 1){
      dist[i] <- euc.dist(X[ret$index.treated[i]], X[ret$index.control[i]])
    }
  }
  n.occur <- data.frame(table(ret$index.treated))
  n.occur <- filter(n.occur, Freq > 1)
  
  stuff <- c()
  for (i in 1:nrow(n.occur)){
    treated.index <- which(ret$index.treated == n.occur[i,1])
    stuff <- c(stuff, treated.index[-which.min(dist[treated.index])])
  }
  ret$index.treated <- ret$index.treated[-stuff]
  ret$index.control <- ret$index.control[-stuff]
  invisible(return(ret))
}