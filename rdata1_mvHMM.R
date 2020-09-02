rdata1_mvHMM<-function(NUM, pii, A, mu_0, Sigma_0, mu_1, Sigma_1)
{

## USAGE
 # rdata1_mvHMM(n, pi, A, ...)



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
    for (i in 1:NUM)
    {
      if (theta[i]==0)
      {
         z[, i]<-mvrnorm(1, mu_0, Sigma_0)
      }
      else
      { 
         z[, i]<-mvrnorm(1, mu_1, Sigma_1)       
      }
    }
    data<-list(s=theta, o=z)
    return(data)

}
