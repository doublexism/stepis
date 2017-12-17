

# Implementations of the likelihood ratio test of Hardy-weinberg equilibrium for uncertain genotype in sibsips.
# Qiong Li, Helene Massam and Xin Gao
# June 27, 2013

# Arguments
# num: a matrix of the numbers of sibships, with genptype configuration (AA), (Aa) and (aa) in rows and columns.
# u: the probability of being incorrectly coded as a "a" allele when "A" is a true allele.
# y: the probability of being incorrectly coded as a "A" allele when "a" is a true allele.

# Examples
# y <- 0.01   #epsilon1
# u <- 0.01   #epsilon2
# c <- c(2,3,0,1,22,10,0,8,54) # under HWE
# c <- c(10,5,1,4,39,13,1,11,16) # without HWE
# num <- matrix(c,nrow =3, ncol=3) # num[1,2] is the number of sibships with genotype (AA, Aa)
# HWEtest.error(num, y, u)

HWEtest.error <- function(num, y, u)
{ 
  
  
  q <- matrix(0, nrow=6, ncol=6) # Penetrance matrix
  
  q[1,1] <- (1-y)^2*(1-y)^2
  q[1,2] <- 1/4*(1-y)^2*(1-y)^2+2*1/4*(1-y)^2*u*(1-y)+1/4*u*(1-y)*u*(1-y)
  q[1,3] <- u*(1-y)*u*(1-y)
  q[1,4] <- 1/16*(1-y)^2*(1-y)^2+2*1/8*(1-y)^2*u*(1-y)+2*1/16*(1-y)^2*u^2+1/4*u*(1-y)*u*(1-y)+2*1/8*u*(1-y)*u^2+1/16*u^2*u^2
  q[1,5] <- 1/4*u*(1-y)*u*(1-y)+2*1/4*u*(1-y)*u^2+1/4*u^2*u^2
  q[1,6] <- u^2*u^2
  
  q[2,1] <- (1-y)^2*2*y*(1-y)
  q[2,2] <- 1/4*(1-y)^2*2*y*(1-y)+1/4*(1-y)^2*(y*u+(1-y)*(1-u))+1/4*u*(1-y)*2*y*(1-y)+1/4*u*(1-y)*(y*u+(1-y)*(1-u))
  q[2,3] <- u*(1-y)*(y*u+(1-y)*(1-u))
  q[2,4] <- 1/16*(1-y)^2*2*y*(1-y)+1/8*(1-y)^2*(y*u+(1-y)*(1-u))+1/8*u*(1-y)*2*y*(1-y)+1/16*(1-y)^2*2*u*(1-u)+
    1/16*u^2*2*y*(1-y)+1/4*u*(1-y)*(y*u+(1-y)*(1-u))+1/8*u*(1-y)*2*u*(1-u)+1/8*u^2*(y*u+(1-y)*(1-u))+1/16*u^2*2*u*(1-u)
  q[2,5] <- 1/4*u*(1-y)*(y*u+(1-y)*(1-u))+1/4*u*(1-y)*2*u*(1-u)+1/4*u^2*(y*u+(1-y)*(1-u))+1/4*u^2*2*u*(1-u)
  q[2,6] <- u^2*2*u*(1-u)
  
  q[3,1] <- (1-y)^2*y^2
  q[3,2] <- 1/4*(1-y)^2*y^2+1/4*(1-y)^2*y*(1-u)+1/4*u*(1-y)*y^2+1/4*u*(1-y)*y*(1-u)
  q[3,3] <- u*(1-y)*y*(1-u)
  q[3,4] <- 1/16*(1-y)^2*y^2+1/8*(1-y)^2*y*(1-u)+1/8*u*(1-y)*y^2+1/16*(1-y)^2*(1-u)^2+1/16*u^2*y^2+1/4*u*(1-y)*y*(1-u)+
    1/8*u*(1-y)*(1-u)^2+1/8*u^2*y*(1-u)+1/16*u^2*(1-u)^2
  q[3,5] <- 1/4*u*(1-y)*y*(1-u)+1/4*u*(1-y)*(1-u)^2+1/4*u^2*y*(1-u)+1/4*u^2*(1-u)^2
  q[3,6] <- u^2*(1-u)^2
  
  q[4,1] <- 2*y*(1-y)*2*y*(1-y)
  q[4,2] <- 1/4*2*y*(1-y)*2*y*(1-y)+2*1/4*2*y*(1-y)*(y*u+(1-y)*(1-u))+1/4*(y*u+(1-y)*(1-u))^2
  q[4,3] <- (y*u+(1-y)*(1-u))^2
  q[4,4] <- 1/16*2*y*(1-y)*2*y*(1-y)+2*1/8*2*y*(1-y)*(y*u+(1-y)*(1-u))+2*1/16*2*y*(1-y)*2*u*(1-u)+1/4*(y*u+(1-y)*(1-u))^2+
    2*1/8*(y*u+(1-y)*(1-u))*2*u*(1-u)+1/16*2*u*(1-u)*2*u*(1-u)
  q[4,5] <- 1/4*(y*u+(1-y)*(1-u))^2+2*1/4*(y*u+(1-y)*(1-u))*2*u*(1-u)+1/4*2*u*(1-u)*2*u*(1-u)
  q[4,6] <- 2*u*(1-u)*2*u*(1-u)
  
  q[5,1] <- 2*y*(1-y)*y^2
  q[5,2] <- 1/4*2*y*(1-y)*y^2+1/4*2*y*(1-y)*y*(1-u)+1/4*(y*u+(1-y)*(1-u))*y^2+1/4*(y*u+(1-y)*(1-u))*y*(1-u)
  q[5,3] <- (y*u+(1-y)*(1-u))*y*(1-u)
  q[5,4] <- 1/16*2*y*(1-y)*y^2+1/8*2*y*(1-y)*y*(1-u)+1/8*(y*u+(1-y)*(1-u))*y^2+1/16*2*y*(1-y)*(1-u)^2+
    1/16*2*u*(1-u)*y^2+1/4*(y*u+(1-y)*(1-u))*y*(1-u)+1/8*(y*u+(1-y)*(1-u))*(1-u)^2+1/8*y*(1-u)*2*u*(1-u)+
    1/16*2*u*(1-u)*(1-u)^2
  q[5,5] <- 1/4*(y*u+(1-y)*(1-u))*y*(1-u)+1/4*(y*u+(1-y)*(1-u))*(1-u)^2+1/4*y*(1-u)*2*u*(1-u)+
    1/4*2*u*(1-u)*(1-u)^2
  q[5,6] <- 2*u*(1-u)*(1-u)^2
  
  q[6,1] <- y^2*y^2
  q[6,2] <- 1/4*y^2*y^2+2*1/4*y^2*y*(1-u)+1/4*y*(1-u)*y*(1-u)
  q[6,3] <- y*(1-u)*y*(1-u)
  q[6,4] <- 1/16*y^2*y^2+2*1/8*y^2*y*(1-u)+2*1/16*y^2*(1-u)^2+1/4*y*(1-u)*y*(1-u)+2*1/8*y*(1-u)*(1-u)^2+1/16*(1-u)^2*(1-u)^2
  q[6,5] <- 1/4*y*(1-u)*y*(1-u)+2*1/4*y*(1-u)*(1-u)^2+1/4*(1-u)^2*(1-u)^2
  q[6,6] <- (1-u)^2*(1-u)^2
  
  
  nu <- rep(0, 6)
  
  nu[1] <- num[1,1]           #the number of sibships with genotype (AA,AA)
  nu[2] <- num[1,2] +num[2,1] #the number of sibships with genotype (AA, Aa) or (Aa, AA)
  nu[3] <- num[1,3] +num[3,1]
  nu[4] <- num[2,2]
  nu[5] <- num[2,3] +num[3,2]
  nu[6] <- num[3,3]
  
  nu
  
  
  m<- 100 # the number of iterations
  
  pt1<- matrix(0, m, 3)
  w<- matrix(0, nrow = 6, ncol = 6)
  pt1[1, ] <- c(0.1, 0.7, 0.2) # initial value
  t<- 1
  while (t<m)
  {
    
    for (i in 1:6)
    {
      
      w[i,1]<- pt1[t,1]*pt1[t,1]*q[i,1]
      w[i,2]<- 2*pt1[t,1]*pt1[t,2]*q[i,2]
      w[i,3]<- 2*pt1[t,1]*pt1[t,3]*q[i,3]
      w[i,4]<- pt1[t,2]*pt1[t,2]*q[i,4]
      w[i,5]<- 2*pt1[t,2]*pt1[t,3]*q[i,5]
      w[i,6]<- pt1[t,3]*pt1[t,3]*q[i,6]
      
    }
    
    w1<- rep(0, 6)
    
    for (i in 1:6)
    {
      w1[i]<- sum(w[i,])
    }
    
    w2<-rep(0, 6)
    for (i in 1:6)
    {
      w2[i]=nu[i]*(2*w[i,1]+w[i,2]+w[i,3])/w1[i]
    }
    
    p1<- sum(w2)
    
    w3<-rep(0, 6)
    for (i in 1:6)
    {
      w3[i]=nu[i]*(2*w[i,4]+w[i,2]+w[i,5])/w1[i]
    }
    p2<- sum(w3)
    
    w4<-rep(0, 6)
    for (i in 1:6)
    {
      w4[i]=nu[i]*(2*w[i,6]+w[i,3]+w[i,5])/w1[i]
    }
    
    p3<- sum(w4)
    
    s<-sum(p1+p2+p3)
    
    t= t+1
    
    pt1[t,]<- c(p1/s, p2/s, p3/s)
    
  }
  
  pt1[m,] # the estimate of genotype frequences
  
  ## EM algorithm for alelles frequences
  
  
  ma<- 100 # the number of iterations
  
  wa<- matrix(0, nrow = 6, ncol = 6)
  pt1a<- matrix(0, ma, 2)
  pt1a[1, ] <- c(0.1, 0.9) # initial value
  t<- 1
  while (t < ma)
  {
    
    
    for (i in 1:6)
    {
      
      wa[i,1]<- (pt1a[t,1]^2)*(pt1a[t,1]^2)*q[i,1]
      wa[i,2]<- 2*(pt1a[t,1]^2)*(2*pt1a[t,1]*pt1a[t,2])*q[i,2]
      wa[i,3]<- 2*(pt1a[t,1]^2)*(pt1a[t,2]^2)*q[i,3]
      wa[i,4]<- (2*pt1a[t,1]*pt1a[t,2])*(2*pt1a[t,1]*pt1a[t,2])*q[i,4]
      wa[i,5]<- 2*(2*pt1a[t,1]*pt1a[t,2])*(pt1a[t,2]^2)*q[i,5]
      wa[i,6]<- (pt1a[t,2]^2)*(pt1a[t,2]^2)*q[i,6]
      
    }
    
    w1a <- rep(0, 6)
    
    for (i in 1:6)
    {
      w1a[i] <- sum(wa[i,])
    }
    
    w2a<-rep(0, 6)
    for (i in 1:6) # i = gi
    {
      w2a[i]=nu[i]*(4*wa[i,1]+3*wa[i,2]+2*wa[i,3]+2*wa[i,4]+wa[i,5])/w1a[i]
    }
    
    p1a<- sum(w2a)
    
    w3a<-rep(0, 6)
    for (i in 1:6)
    {
      w3a[i]=nu[i]*(wa[i,2]+2*wa[i,3]+2*wa[i,4]+3*wa[i,5]+4*wa[i,6])/w1a[i]
    }
    p2a<- sum(w3a)
    
    sa<-sum(p1a+p2a)
    
    t= t+1
    
    pt1a[t,]<- c(p1a/sa, p2a/sa)
    
  }
  
  pt1a[ma,] # the estimate of alleles frequences
  
  ##### Likelihood ratio test
  
  ### caculate unconstrained log - likelihood function
  
  l <- matrix(0, nrow = 6, ncol = 6)
  
  for (i in 1:6)
  {
    
    l[i,1]<- pt1[m,1]*pt1[m,1]*q[i,1]
    l[i,2]<- 2*pt1[m,1]*pt1[m,2]*q[i,2]
    l[i,3]<- 2*pt1[m,1]*pt1[m,3]*q[i,3]
    l[i,4]<- pt1[m,2]*pt1[m,2]*q[i,4]
    l[i,5]<- 2*pt1[m,2]*pt1[m,3]*q[i,5]
    l[i,6]<- pt1[m,3]*pt1[m,3]*q[i,6]
    
  }
  
  l1<- rep(0, 6)
  
  for (i in 1:6)
  {
    l1[i]<- log(sum(l[i,]))
  }
  
  
  ln1 <- sum(nu*l1) # unconstrained log - likelihood function
  
  ### calculate log - likelihood function under HWE
  
  la <- matrix(0, nrow = 6, ncol = 6)
  
  for (i in 1:6)
  {
    
    la[i,1]<- (pt1a[ma,1]^2)*(pt1a[ma,1]^2)*q[i,1]
    la[i,2]<- 2*(pt1a[ma,1]^2)*(2*pt1a[ma,1]*pt1a[ma,2])*q[i,2]
    la[i,3]<- 2*(pt1a[ma,1]^2)*(pt1a[ma,2]^2)*q[i,3]
    la[i,4]<- (2*pt1a[ma,1]*pt1a[ma,2])*(2*pt1a[ma,1]*pt1a[ma,2])*q[i,4]
    la[i,5]<- 2*(2*pt1a[ma,1]*pt1a[ma,2])*(pt1a[ma,2]^2)*q[i,5]
    la[i,6]<- (pt1a[ma,2]^2)*(pt1a[ma,2]^2)*q[i,6]
    
  }
  
  l1a <- rep(0, 6)
  
  for (i in 1:6)
  {
    l1a[i] <- log(sum(la[i,]))
  }
  
  
  ln0 <- sum(nu*l1a) # constrained log - likelihood function 
  
  kai <- -2*(ln0 - ln1)
  
  pv  <- 1 - pchisq(kai, df = 1)
  
  return(list(kai, pv))}