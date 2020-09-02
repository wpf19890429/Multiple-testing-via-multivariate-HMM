bwfw_mvHMM<-function(z, A, mu_0, Sigma_0, pc, mu_1, Sigma_1)
{


## DETAILS
 # bwfw_mvHMM calculates values for backward, forward variables, probabilities of hidden states, 
 # --mvLIS values and etc.

## VALUES
 # alpha: rescaled forward variables
 # beta: rescaled backward variables
 # mvLIS: multivariate LIS values


################################################################## 

## Initialize
    NUM<-ncol(z)
    p<-nrow(z)
    pii<-c(0.5, 0.5)

## Densities
    f0z<-rep(0, NUM)
    for(i in 1:NUM)
    {
        f0z[i]<-dmvnorm(z[, i], mean = mu_0, sigma = Sigma_0, log = FALSE) 
    }
    f1z<-rep(0, NUM)
    nc<-length(pc)
    for(i in 1:NUM)
    {
        for(c in 1:nc)
        {
              f1z[i]<-f1z[i]+pc[c]*dmvnorm(z[, i], mean = mu_1[c, ], sigma = Sigma_1[, , c], log = FALSE) 
        }
    }


## the backward-forward procedure

# a. the forward variables
# --rescaled 
    alpha<-matrix(0, NUM, 2, byrow=TRUE)

# scaling variable c_0
    c0<-rep(0, NUM)
    alpha[1, 1]<-pii[1]*f0z[1]
    alpha[1, 2]<-pii[2]*f1z[1]

# rescaling alpha
    c0[1]<-1/sum(alpha[1, ])
    alpha[1, ]<-c0[1]*alpha[1, ]

    for (k in 1:(NUM-1))
    { 
      alpha[k+1, 1]<-(alpha[k, 1]*A[1, 1]+alpha[k, 2]*A[2, 1])*f0z[k+1]
      alpha[k+1, 2]<-(alpha[k, 1]*A[1, 2]+alpha[k, 2]*A[2, 2])*f1z[k+1]
# rescaling alpha
      c0[k+1]<-1/sum(alpha[k+1, ])
      alpha[k+1, ]<-c0[k+1]*alpha[k+1, ]
    }

# b. the forward variables
# --rescaled
    beta<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)

    beta[NUM, 1]<-c0[NUM]
    beta[NUM, 2]<-c0[NUM]

    for (k in (NUM-1):1)
    { 
      beta[k, 1]<-A[1, 1]*f0z[k+1]*beta[k+1, 1]+A[1, 2]*f1z[k+1]*beta[k+1, 2]
      beta[k, 2]<-A[2, 1]*f0z[k+1]*beta[k+1, 1]+A[2, 2]*f1z[k+1]*beta[k+1, 2]
# rescaling beta
# using the same scaling factors as alpha 
      beta[k, ]<-c0[k]*beta[k, ]
    }

# c. cmLIS values
# --original
# --the same formulae hold for the rescaled alpha and beta

    mLIS<-rep(0, NUM)
    for (k in 1:NUM)
    { 
      q1<-alpha[k, 1]*beta[k, 1]
      q2<-alpha[k, 2]*beta[k, 2]
      mLIS[k]<-q1/(q1+q2)
    }

# d. return the results of the bwfw proc.

    bwfw.res<-list(bw=alpha, fw=beta, mLIS=mLIS, c0=c0)
    return(bwfw.res)
  
}



