rdata_mvHMM<-function(NUM, pii, A, mu_0, Sigma_0, pc, mu_1, Sigma_1)
{

## USAGE
 # rdata_mvHMM(n, pi, A, ...)


## Initialize
    theta<-rep(0, NUM)
    p<-length(mu_0)
    z<-matrix(0, nrow=p, ncol=NUM)

## Generating the states
 # initial state
    theta[1]<-rbinom(1, 1, pii[2])
 # other states
    for (i in 2:NUM)
    {
      if (theta[i-1]==0)
         theta[i]<-rbinom(1, 1, A[1, 2])
      else
         theta[i]<-rbinom(1, 1, A[2, 2])
    }

## Generating the observations

    nc<-length(pc)
    for (i in 1:NUM)
    {
      if (theta[i]==0)
      {
         z[, i]<-mvrnorm(1, mu_0, Sigma_0)
      }
      else
      { 
         c<-sample(1:nc, 1, prob=pc)
         mu_1c<-mu_1[c, ]
         Sigma_1c<-Sigma_1[, , c]
         z[, i]<-mvrnorm(1, mu_1c, Sigma_1c)       
      }
    }
    data<-list(s=theta, o=z)
    return(data)

}
