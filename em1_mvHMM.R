em1_mvHMM<-function(z, maxiter=200)
{

############################################################


    NUM<-ncol(z)
# precision tolerance level
    ptol<-1e-4
    niter<-0


### initializing model parameters
    pii.new<-c(0, 1)
    A.new<-matrix(0.5, 2, 2)
    mu_0<-c(0, 0)
    Sigma_0<-matrix(c(1, 0, 0, 1), 2, 2, byrow=TRUE)
##  Sigma_0.new<-matrix(c(1, 0.2, 0.2, 1), 2, 2, byrow=TRUE)
    mu_1.new<-c(1, 1)
    Sigma_1.new<-matrix(c(1, 0.2, 0.2, 1), 2, 2, byrow=TRUE)

    diff<-10
    Loglikelihood.new<--10000

### The E-M Algorithm
    while(diff>ptol && niter<maxiter)
    {
        niter<-niter+1
        pii.old<-pii.new
        A.old<-A.new
##      Sigma_0.old<-Sigma_0.new
        mu_1.old<-mu_1.new
        Sigma_1.old<-Sigma_1.new
        Loglikelihood.old<-Loglikelihood.new

## updating the weights and probabilities of hidden states
        bwfw.res<-bwfw1_mvHMM(z, A.old, mu_0, Sigma_0, mu_1.old, Sigma_1.old)


# the backward-forward variable
        alpha<-bwfw.res$bw
        beta<-bwfw.res$fw


## Densities
        f0z<-rep(0, NUM)
        for(i in 1:NUM)
        {
             f0z[i]<-dmvnorm(z[, i], mean = mu_0, sigma = Sigma_0, log = FALSE) 
        }
        f1z<-rep(0, NUM)
        for(i in 1:NUM)
        {
             f1z[i]<-dmvnorm(z[, i], mean = mu_1.old, sigma = Sigma_1.old, log = FALSE) 
        }

### example: xi[i, p, q]=Pr(theta_i=p-1, theta_{i+1}=q-1 | z_1^m)
        xi<-array(0, c(NUM-1, 2, 2))
        b1<-rep(0, NUM-1)
        for(i in 1:(NUM-1))
        {
             xi[i, 1, 1]<-alpha[i, 1]*beta[i+1, 1]*f0z[i+1]*A.old[1, 1]
             xi[i, 1, 2]<-alpha[i, 1]*beta[i+1, 2]*f1z[i+1]*A.old[1, 2]
             xi[i, 2, 1]<-alpha[i, 2]*beta[i+1, 1]*f0z[i+1]*A.old[2, 1]
             xi[i, 2, 2]<-alpha[i, 2]*beta[i+1, 2]*f1z[i+1]*A.old[2, 2]
             b1[i]<-1/(xi[i, 1, 1]+xi[i, 1, 2]+xi[i, 2, 1]+xi[i, 2, 2])
             xi[i, 1, 1]<-b1[i]*xi[i, 1, 1]
             xi[i, 1, 2]<-b1[i]*xi[i, 1, 2]             
             xi[i, 2, 1]<-b1[i]*xi[i, 2, 1]
             xi[i, 2, 2]<-b1[i]*xi[i, 2, 2]
        }
        
# b. transition matrix A 
        for(i in 1:2)
            for(j in 1:2)
            {
                 q1<-sum(xi[, i, j])
                 q2<-sum(xi[, i, 1])+sum(xi[, i, 2])
                 A.new[i, j]<-q1/q2
            }

# c. non-null distribution 
        ptheta<-rep(0, NUM)
        ptheta<-alpha[, 2]*beta[, 2]/(alpha[, 1]*beta[, 1]+alpha[, 2]*beta[, 2])
        q1<-sum(ptheta)
        q2<-sum(ptheta*z[1, ])
        q3<-sum(ptheta*z[2, ])
        mu_1.new<-c(q2/q1, q3/q1)
##      Sigma_0.new<-matrix(0, 2, 2)
##      for(i in 1:NUM)
##      {
##           Sigma_0.new<-Sigma_0.new+(1-ptheta[i])*matrix(z[, i]-mu_0, ncol=1)%*%matrix(z[, i]-mu_0, nrow=1)/(NUM-q1)
##      }
        Sigma_1.new<-matrix(0, 2, 2)
        for(i in 1:NUM)
        {
             Sigma_1.new<-Sigma_1.new+ptheta[i]*matrix(z[, i]-mu_1.new, ncol=1)%*%matrix(z[, i]-mu_1.new, nrow=1)/q1
        }                           
        c0<-bwfw.res$c0
        Loglikelihood.new<--sum(log(c0))
        df1<-abs(Loglikelihood.old-Loglikelihood.new)
        diff<-df1
    }


# g. return the results of the E-M algorithm
    em.var<-list(pii=pii.new, A=A.new, mu_1=mu_1.new, Sigma_1=Sigma_1.new, ni=niter)
    return(em.var)

}
