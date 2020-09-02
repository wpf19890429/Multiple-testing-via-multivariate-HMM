em_mvHMM<-function(z, L=2, maxiter=200)
{

############################################################


    NUM<-ncol(z)
# precision tolerance level
    ptol<-1e-4
    niter<-0


### initializing model parameters
    A.new<-matrix(0.5, 2, 2)
    mu_0<-c(0, 0)
    Sigma_0<-matrix(c(1, 0, 0, 1), 2, 2, byrow=TRUE)
#   Sigma_0.new<-matrix(c(1, 0.2, 0.2, 1), 2, 2, byrow=TRUE)
    pc.new<-rep(1/L, L)
#   mu_1.new<-matrix(1, nrow=L, ncol=2)
    mu_1.new<-matrix(c(2, 2, 2, 2), nrow = 2, ncol = 2, byrow = TRUE)
    Sigma_1.new<-array(rep(c(1, 0.2, 0.2, 1), L), dim=c(2, 2, L))
    diff<-10
    Loglikelihood.new<--10000

### The E-M Algorithm
    while(diff>ptol && niter<maxiter)
    {
        niter<-niter+1
        A.old<-A.new
#       Sigma_0.old<-Sigma_0.new
        pc.old<-pc.new
        mu_1.old<-mu_1.new
        Sigma_1.old<-Sigma_1.new
        Loglikelihood.old<-Loglikelihood.new

## updating the weights and probabilities of hidden states
        bwfw.res<-bwfw_mvHMM(z, A.old, mu_0, Sigma_0, pc.old, mu_1.old, Sigma_1.old)


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
             for(c in 1:L)
             {
                  f1z[i]<-f1z[i]+pc.old[c]*dmvnorm(z[, i], mean = mu_1.old[c, ], sigma = Sigma_1.old[, , c], log = FALSE) 
             }
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

# c. intermediate variable
        ptheta<-rep(0, NUM)
        ptheta<-alpha[, 2]*beta[, 2]/(alpha[, 1]*beta[, 1]+alpha[, 2]*beta[, 2])


# c. weight variables
        omega<-matrix(0, nrow=NUM, ncol=L, byrow=TRUE)
        for (c in 1:L)
        { 
          for (i in 1:NUM)
          {

                 f1z.c<-dmvnorm(z[, i], mean = mu_1.old[c, ], sigma = Sigma_1.old[, , c], log = FALSE)
                 omega[i, c]<-ptheta[i]*pc.old[c]*f1z.c/f1z[i]
          }
        }

# updating Sigma_0
#       Sigma_0.new<-matrix(0, 2, 2)
#       for(i in 1:NUM)
#       {
#            Sigma_0.new<-Sigma_0.new+(1-ptheta[i])*matrix(z[, i]-mu_0, ncol=1)%*%matrix(z[, i]-mu_0, nrow=1)/(NUM-sum(ptheta))
#       }

# d. updating the non-null distribution 

 
     for (c in 1:L)
     {
# (i). probability weights
        q1<-sum(omega[, c])
        q2<-sum(ptheta)
        pc.new[c]<-q1/q2


# (ii). means
        q3<-sum(omega[, c]*z[1, ])
        q4<-sum(omega[, c]*z[2, ])
        mu_1.new[c, ]<-c(q3/q1, q4/q1)

# (iii). sds
        Sigma_1.new[, , c]<-matrix(0, nrow=2, ncol=2)
        for(i in 1:NUM)
        {
             Sigma_1.new[, , c]<-Sigma_1.new[, , c]+omega[i, c]*matrix(z[, i]-mu_1.new[c, ], ncol=1)%*%matrix(z[, i]-mu_1.new[c, ], nrow=1)/q1
        }  
     }                         
     c0<-bwfw.res$c0
     Loglikelihood.new<--sum(log(c0))
     df1<-abs(Loglikelihood.old-Loglikelihood.new)
     diff<-df1
    }

# g. return the results of the E-M algorithm
    em.var<-list(A=A.new, pc=pc.new, mu_1=mu_1.new, Sigma_1=Sigma_1.new, ni=niter)
    return(em.var)

}
