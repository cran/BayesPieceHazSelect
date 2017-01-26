#' PiecewiseBayesSelect
#' @param Y1  Vector Containing  event times (or censoring time due to death/censoring)
#' @param I1 Vector Containing  event indicators (1 if l event for a patient, 0 otherwise)
#' @param X Matrix of Patient Covariates, the last inc are left out of the selection procedure
#' @param hyperparameters  List containing 11 hyperparameters and four starting values. In order they are: psi-the swap rate of the SVSS algorithm.
#'  c-parameter involved in Sigma matrix for selection. z1a, z1b - beta hyper parameters on probability of inclusion for each of the three hazard functions.
#'  a1,b1- hyperparameters on sigma_lambda.
#'   clam1- spatial dependency of baseline hazard (between 0 and 1) for the  hazard function.
#'    Alpha1 - The parameter for the number of split points in the hazard (must be whole number).
#'    J1max - Maximum number of split points allowed (must be whole number).
#'    J1- Starting number of split points.  cl1 -Tuning parameter for log baseline hazard height sampler.
#' @param beta1start  Starting Values for Beta1
#' @param B  Number of iterations
#' @param inc  Number of variables left out of selection
#' @param Path  Where to save posterior samples
#' @param burn  percent of posterior sample to burn in (burn*B must be a whole number)
#'@import graphics
#'@import stats
#'@import mvtnorm
#'@import utils
#' @examples
#' ##Randomly Generate Semicompeting Risks Data
#' ####Generates random patient time, indicator and covariates.
#' n=100
#' Y1=runif(n,0,100)
#' I1=rbinom(n,1,.5)
#' library(mvtnorm)
#' X=rmvnorm(n,rep(0,13),diag(13))
#' ####Read in Hyperparameters
#' ##Swap Rate
#' psi=.5
#' c=20
#' ###Eta Beta function probabilities
#' z1a=.4
#' z1b=1.6
#' ####Hierarchical lam params
#' ###Sigma^2 lambda_ hyperparameters
#' a1=.7
#' b1=.7
#' ##Spacing dependence c in [0,1]
#' clam1=1
#' #####NumSplit
#' alpha1=3
#' J1max=10
#' ####Split Point Starting Value ###
#' J1=3
#' ##Tuning parameter for lambda
#' cl1=.25
#' ###Beta Starting Values
#' beta1start=c(0,0,-1,0,0,0,1,1,1,1,1,-1,-1)
#' hyper=c(psi,c,z1a,z1b,a1,b1,clam1,alpha1,J1max,J1,cl1)
#' ###Number of iterations and output location
#' B=200
#'Path=tempdir()
#'inc=2
#'burn=.4
#' PiecewiseBayesSelect(Y1,I1,X,hyper,beta1start,B,inc,Path,burn)
#' @export
PiecewiseBayesSelect=function(Y1,I1,X,hyperparameters,beta1start,B,inc,Path,burn){



  if(inc%%1>0){
    cat("inc must be a natural number")
  }else{




    ####Hyperparameters##
    ##Swap Rate
    psi=hyperparameters[1]
    ##
    c=hyperparameters[2]
    ###Eta Beta function probabilities
    z1a=hyperparameters[3]
    z1b=hyperparameters[4]

    ####Hierarchical lam params
    ###Siglam
    a1=hyperparameters[5]
    b1=hyperparameters[6]

    ##Spacing dependence c in [0,1]
    clam1=hyperparameters[7]

    ##Lamsampler params
    #####NumSplit
    alpha1=hyperparameters[8]

    J1max=hyperparameters[9]

    ####Split Points###
    J1=hyperparameters[10]



    cl1=hyperparameters[11]




    p1=ncol(X)-inc

    n=length(Y1)











    #####In program
    ###Make Acceptance Matrices
    ###Beta/Eta###
    beta1=matrix(rep(1,B*(p1+inc)),nrow=B)
    eta1=matrix(rep(1,B*p1),nrow=B)
    ####Frailty Matrix###
    ###
    Mulam1=rep(0,B)
    Siglam1=rep(1,B)






    ###Make Eta1Start
    beta1[1,]=beta1start

    ##
    eta1start=rep(1,p1)


    for(i in 1:p1){
      if(beta1start[i]==0){
        eta1start[i]=0
      }
    }


    eta1[1,]=eta1start



    m1 = max(Y1[I1==1])+.001




    ####Acceptance Matrices


    Acceptlam1=matrix(rep(NA,B*(J1max+1)),nrow=B)

    accepts1=rep(0,B)

    Indmix1=rep(0,B)

    sum1=rep(0,B)

    split1=rep(0,B)


    Indcond1=matrix(rep(NA,p1*B),nrow=B)



    #########################S Matrices!!!
    #Reset up lam and S1 matrices
    s1=matrix(rep(NA,B*(J1max+2)),nrow=B)
    s1[1,1:(J1+2)]=sort(seq(0,m1,length.out = J1+2))

    lam1=matrix(rep(NA,B*(J1max+1)),nrow=B)

    lam1[1,1:(J1+1)]=rep(0,J1+1)

    ###Acceptance
    split1=rep(0,B)

    IndB1=rep(0,B)

    ###Death
    IndD1=rep(0,B)

    Indeta1=rep(0,B)

    Ind1s=rep(0,B)








    n=length(Y1)
    G1=J1+1






    #####
    LK1L=function(Y1,I1,X,Beta1,s1,lam1){

      LOGBH=0
      et1=X%*%Beta1


      for(k in 1:G1){


        Del=pmax(0,pmin(Y1,s1[k+1])-s1[k])



        LOGBH=LOGBH-sum(Del*exp(lam1[k])*exp(et1))

        zu=Y1<=s1[k+1]
        zl=Y1>s1[k]
        LOGBH=LOGBH+sum(zu*zl*I1)*lam1[k]
      }



      return(LOGBH)

    }
    ###Haz 2


    ###



    #####
    LK1=function(Y1,I1,X,Beta1,s1,lam1){

      LOGBH=0
      et1=X%*%Beta1


      for(k in 1:G1){


        Del=pmax(0,pmin(Y1,s1[k+1])-s1[k])



        LOGBH=LOGBH-sum(Del*exp(lam1[k])*exp(et1))

      }

      LOGBH=LOGBH+sum(I1*et1)



      return(LOGBH)

    }
    ###Haz 2







    if(inc>1){
      cat("More than One Variable Included", "


          ")

      ###Set Up Additional Acceptance Matrix

      IncCond1=matrix(rep(0,B*inc),nrow=B)



      iter=c(0,0)


      ##Sampler


      for(b in 2:B){





        if(b%%10000==0){cat(b, "iterations",date(), "  ")}else{
          if(b%%5000==0){cat(b, " iterations ")}}

        U=runif(1,0,1)


        iter[1]="etabeta1"

        ###eta1,beta1
        eta1[b,]=eta1[b-1,]
        beta1[b,]=beta1[b-1,]

        if(sum(eta1[b-1,])==0|sum(eta1[b-1,])==p1){
          if(sum(eta1[b-1,])==0){
            ###Add Automatically
            iter[2]="Add"
            Ind=sample(1:p1,1)
            eta1[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta1[b-1,c(includednew, (p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta1[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta1[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])


            alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))

            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }
          }
          if(sum(eta1[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            iter[2]="delete"
            eta1[b,Ind]=0
            beta1[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta1[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }}
        }else{

          U=runif(1,0,1)

          if(U<psi){
            ###Swapper
            includedold=rep(0,p1)
            iter[2]="swap"
            for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            ones=includedold
            zeros=rep(0,p1)
            for(k in 1:p1){if(eta1[b-1,k]==0){zeros[k]=k}}
            zeros=zeros[zeros != 0]
            ###Sample swap indices###
            if(length(ones)==1){
              Indone=ones}else{
                Indone=sample(ones,1)}
            if(length(zeros)==1){Indzero=zeros}else{
              Indzero=sample(zeros,1)}
            ####Change Beta/eta
            eta1[b,Indone]=0
            eta1[b,Indzero]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
            spot1=max(spotold)
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ###Generate new vector##
            beta1[b,Indone]=0

            ##meannew,varnew##
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta1[b-1,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta1[b,Indzero]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta1[b,Indzero],meannew,varnew))
            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta1[b-1,Indone],meanold,varold))
            ######Accept reject###
            Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

            alphab1=Liken-Likeo+dn-do
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{
              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}

            }

          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta1[b-1,Ind]==1){
              ##delete##
              iter[2]="delete"
              eta1[b,Ind]=0
              beta1[b,Ind]=0
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)


              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta1[b-1,Ind],meanold,varold))
              ######Accept reject###
              Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
              Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

              alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{
                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}

              }

            }else{
              ###Add###
              eta1[b,Ind]=1

              iter[2]="add"
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta1[b-1,c(includednew,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta1[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta1[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
              Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{
                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}
              }


            }

          }}



        iter[1]="Beta1"
        iter[2]="Included"

        if(sum(eta1[b,])==0){
          ##Sample Included
          Sigmanew= c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])
          zeta1n=beta1[b,(p1+1):(p1+inc)]
          for(k in 1:inc){
            zeta1=zeta1n
            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetano=zeta1[-k]
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            zeta1[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1n[k],meannew,varnew))
            beta=beta1[b,]
            beta[(p1+1):(p1+inc)]=zeta1
            Likeo=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta,s1[b-1,],lam1[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))

            if(is.finite(alphab1m)==FALSE){
              IncCond1[b,k]=0
            }else{
              if(U>alphab1m){
                IncCond1[b,k]=0
              }else{IncCond1[b,k]=1
              beta1[b,]=beta
              zeta1n=zeta1
              }}
            ##End Inc Sampler
          } }else{
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            zeta1n=beta1[b,c(includednew,(p1+1):(p1+inc))]

            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            p=length(includednew)+inc
            ####Update All included variables

            for(k in (length(includednew)+1):(length(includednew)+inc)){
              zeta1=zeta1n
              V1 = Sigmanew[k,k]
              V2 = as.matrix(Sigmanew[-k,-k])
              V12 = as.matrix(Sigmanew[k,-k])
              thetano = as.matrix(zeta1[-k])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              zeta1[k]=rnorm(1,meannew,varnew)
              dn=log(dnorm(zeta1[k],meannew,varnew))
              ###density old
              do=log(dnorm(beta1[b,(p1+k-length(includednew))],meannew,varnew))


              ######Accept reject###
              Likeo=LK1(Y1,I1,X,c(beta1[b,1:p1],zeta1n[(length(zeta1n)-inc+1):length(zeta1n)]), s1[b-1,],lam1[b-1,])
              Liken=LK1(Y1,I1,X,c(beta1[b,1:p1],zeta1[(length(zeta1n)-inc+1):length(zeta1n)]), s1[b-1,],lam1[b-1,])

              alphab1s=Liken-Likeo+dn -do
              U=log(runif(1,0,1))

              if(is.finite(alphab1s)==FALSE){
                IncCond1[b,(k-p1)]=0

              }else{

                if(U>alphab1s){



                  IncCond1[b,(k-p1)]=0

                }else{IncCond1[b,(k-p1)]=1
                zeta1n=zeta1
                beta1[b,]=c(beta1[b,1:p1],zeta1[(length(zeta1)-inc+1):length(zeta1)])

                }
              }

            }

            ###End included sampler###
          }



        #####Conditional Sampler for Included!###


        if(sum(eta1[b,])>0){

          iter[2]="Conditional Inclusion"
          ##Jointly Update nonzero betas
          zeta1=beta1[b,]
          zeta1=zeta1[zeta1!=0]
          zeta1n=zeta1
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])



          ###############
          ####

          for(k in 1:length(includednew)){


            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetab=beta1[b,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1n[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1n[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1[k],meannew,varnew))






            beta=beta1[b,]
            beta[c(includednew,(p1+1):(p1+inc))]=zeta1n


            Likeo=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta,s1[b-1,],lam1[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond1[b,k]=0
            }else{
              if(U>alphab1m){
                Indcond1[b,includednew[k]]=0
                zeta1n[k]=zeta1[k]
              }else{Indcond1[b,includednew[k]]=1
              beta1[b,]=beta
              zeta1[k]=zeta1n[k]
              }}

          }



          ##Jointly Update nonzero betas
          iter[2]="mixing"
          zeta1n=beta1[b,]
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])


          zeta1n[c(includednew,(p1+1):(p1+inc))]=rmvnorm(1,rep(0,length(includednew)+inc),Sigmanew)

          beta=beta1[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[c(includednew,(p1+1):(p1+inc))],rep(0,length(includednew)+inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,length(includednew)+inc),Sigmanew))

          ######Accept reject###
          Likeo=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])
          Liken=LK1(Y1,I1,X,zeta1n, s1[b-1,],lam1[b-1,])

          alphamix1=Liken-Likeo+dn -do
          U=log(runif(1,0,1))
          if(is.finite(alphamix1)==FALSE){
            Indmix1[b]=0
          }else{
            if(U>alphamix1){

              Indmix1[b]=0
            }else{Indmix1[b]=1
            beta1[b,]=zeta1n
            }}

        }else{
          ##Jointly Update nonzero betas
          iter[2]="mixing No eta"
          zeta1n=beta1[b,]
          Sigmanew=c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])


          zeta1n[(p1+1):(p1+inc)]=rmvnorm(1,rep(0,inc),Sigmanew)

          beta=beta1[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[(p1+1):(p1+inc)],rep(0,inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,inc),Sigmanew))

          ######Accept reject###
          Likeo=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])
          Liken=LK1(Y1,I1,X,zeta1n,s1[b-1,],lam1[b-1,])

          alphamix1=Liken-Likeo+dn -do
          U=log(runif(1,0,1))

          if(is.finite(alphamix1)==FALSE){
            Indmix1[b]=0}else{

              if(U>alphamix1){



                Indmix1[b]=0
              }else{Indmix1[b]=1
              beta1[b,]=zeta1n
              }}

        }



        S1=s1[b-1,]
        S1=S1[!is.na(S1)]


        L1=lam1[b-1,]
        L1=as.matrix(L1[!is.na(L1)])


        ############################################
        #####Start LogBH Samplers###################
        ############################################
        ####Lam1####

        iter[1]="LogBH1"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
        Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


        length1=rep(0,J1+1)




        for(j in 1:length(length1)){
          length1[j]=s1[b-1,j+1]-s1[b-1,j]
        }


        if(J1<2){
          if(J1==1){
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            SigLam1=solve(diag(J1+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m1))
            SigLam1=Q1
          }
        }else{


          for(j in 2:J1){
            W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


          SigLam1=solve(diag(J1+1)-W1)%*%Q1

        }



        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam

        if(J1>0){

          Mulam1[b]=rnorm(1,(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%L1)/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1))),sqrt(Siglam1[b-1]/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1)))))


          Siglam1[b]=1/rgamma(1,a1+(J1+1)/2,b1+.5*(t(as.matrix(rep(Mulam1[b],J1+1))-L1)%*%solve(SigLam1)%*%(as.matrix(rep(Mulam1[b],J1+1))-L1)))


          ##Siglam

          iter[2]="Sigma"
        }else{



          Mulam1[b]=rnorm(1,lam1[b-1,1],sqrt(Siglam1[b-1]))


          Siglam1[b]=1/rgamma(1,a1+1/2,b1+.5*(Mulam1[b]-lam1[b-1,1])^2)



        }


        #lambda1
        iter[2]="lam1"
        lam1[b,]=lam1[b-1,]
        #######

        for(m in 1:(J1+1)){



          lam=lam1[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl1,cl1)


          if(J1==0){
            do=log(dnorm(lambda[m],Mulam1[b],sqrt(Siglam1[b])))
            dn=log(dnorm(lam[m],Mulam1[b],sqrt(Siglam1[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
            do=dmvnorm(lam,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
          }

          Likeo=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b-1,],lam1[b,])
          Liken=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b-1,],lam)




          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam1[b,m]=lam1[b-1,m]
            Acceptlam1[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam1[b,m]=1
              lam1[b,m]=lam[m]
            }else{Acceptlam1[b,m]=0}
          }


        }




        #####################################################
        ###################################################

        iter[1]="Haz1"
        iter[2]="Birth"

        ###Random Perturbation###
        U1=runif(1,0,1)
        #####

        s=s1[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J1max){
          Birth=runif(1,0,m1)

          s1[b,1:(J1+3)]=sort(c(s,Birth))

          for(k in 2:(J1+2)){
            if(Birth>s1[b-1,k-1] & Birth<s1[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J1+2)

          if(Ind==1 | Ind==J1+1){
            if(Ind==1){
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
            }else{
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam1[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b-1,],lam1[b,])

          if(J1>0){
            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))
          }else{
            do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))
          }

          prior=((2*J1+3)*(2*J1+2)*(Birth-s1[b-1,Ind])*(s1[b-1,Ind+1]-Birth))/((m1^2)*(s1[b-1,Ind+1]-s1[b-1,Ind]))

          G1=G1+1
          J1=J1+1

          Ln=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b,],lam)


          ##Make SigLam1



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1





            }else{

              SigLam1n=2/m1
            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

          }



          dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),Siglam1[b]*SigLam1n))




          alpha=Ln-Lo+dn-do-log(U1*(1-U1)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB1[b]=0
            s1[b,]=s1[b-1,]
            J1=J1-1
            G1=G1-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB1[b]=1
              lam1[b,1:(J1+1)]=lam
            }else{
              s1[b,]=s1[b-1,]
              IndB1[b]=0
              J1=J1-1
              G1=G1-1
            }

          }


        }else{
          s1[b,]=s1[b-1,]
          IndB1[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U1=runif(1,0,1)

        if(J1==0){
          IndD1[b]=0
          s1[b,]=s1[b-1,]
        }else{

          if(J1==1){
            Ind=2
          }else{

            Ind=sample(2:(J1+1),1)
          }


          s=s1[b,]
          s=s[-Ind]

          lam=lam1[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s1[b,Ind]-s1[b,Ind-1])*lam1[b,Ind-1]+(s1[b,Ind+1]-s1[b,Ind])*lam1[b,Ind])/(s1[b,Ind+1]-s1[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])



          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1=solve(diag(J1+1)-W1)%*%Q1

              do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))




            }else{


              do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))

            }
          }else{

            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1=solve(diag(J1+1)-W1)%*%Q1

            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))


          }
          #############################################
          #############################################

          Lo=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b,],lam1[b,])


          prior=((m1^2)*(s1[b,Ind+1]-s1[b,Ind-1]))/((2*J1+1)*(2*J1)*(s1[b,Ind]-s1[b,Ind-1])*(s1[b,Ind+1]-s1[b,Ind]))


          G1=G1-1
          J1=J1-1


          Ln=LK1L(Y1,I1,X,as.matrix(beta1[b,]), s,lam)

          ###Make siglam matrix



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)



          length1=rep(0,J1+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1


              dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))



            }else{

              SigLam1n=2/m1
              dn=log(dpois(J1,alpha1))+log(dnorm(lam,Mulam1[b],Siglam1[b]))

            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

            dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U1*(1-U1))

          if(is.nan(alpha)==TRUE){
            IndD1[b]=0
            J1=J1+1
            G1=G1+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s1[b,]=c(s,NA)
              IndD1[b]=1
              lam1[b,1:(J1+1)]=lam
              lam1[b,(J1+2):J1max]=rep(NA,J1max-J1-1)
            }else{
              IndD1[b]=0
              J1=J1+1
              G1=G1+1
            }
          }

          ####End else
        }
        ##






        split1[b]=J1

        ##
        sum1[b]=sum(eta1[b,])



      }





      ################End Samplers
      cat(c,z1a,z1b,"




          ", "



          ", "


          ", "Posterior Inclusion Probabilities after half Burnin", "


          ", "Hazard 1", "

          ", colMeans(eta1[(B*burn+1):B,])*100, "


          ", "IndEta",mean(Indeta1[(B*burn+1):B])*100,"


          ","IndMix",mean(Indmix1[(B*burn+1):B])*100,"


          ", "Included Acceptance", "

          ", "Haz1", "
          ", colMeans(IncCond1[(B*burn+1):B,])*100, "

          ", colMeans(Indcond1[(B*burn+1):B,],na.rm=TRUE)*100,"


          ","Survival","

          ","IndDeath",mean(IndD1[(B*burn+1):B])*100,"

          ","IndBirth",mean(IndB1[(B*burn+1):B])*100,"

          ","Lambda","

          ", "Lam1",
          colMeans(Acceptlam1[(B*burn+1):B,],na.rm=TRUE)*100)





      Path1= paste0(Path,"/IncCond1.txt")

      write.table(IncCond1[(burn*B+1):B,], Path1, sep="\t")







    }




    if(inc==1){
      cat("One Variable Included")

      Ind1s=rep(0,B)


      for(b in 2:B){



        if(b%%10000==0){cat(b, "iterations",date(), "  ")}else{
          if(b%%5000==0){cat(b, " iterations ")}}

        U=runif(1,0,1)


        iter[1]="etabeta1"

        ###eta1,beta1
        eta1[b,]=eta1[b-1,]
        beta1[b,]=beta1[b-1,]

        if(sum(eta1[b-1,])==0|sum(eta1[b-1,])==p1){
          if(sum(eta1[b-1,])==0){
            ###Add Automatically
            iter[2]="Add"
            Ind=sample(1:p1,1)
            eta1[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ####
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta1[b-1,c(includednew, (p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta1[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta1[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))

            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }
          }
          if(sum(eta1[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            iter[2]="delete"
            eta1[b,Ind]=0
            beta1[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta1[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }}
        }else{

          U=runif(1,0,1)

          if(U<psi){
            ###Swapper
            includedold=rep(0,p1)
            iter[2]="swap"
            for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            ones=includedold
            zeros=rep(0,p1)
            for(k in 1:p1){if(eta1[b-1,k]==0){zeros[k]=k}}
            zeros=zeros[zeros != 0]
            ###Sample swap indices###
            if(length(ones)==1){
              Indone=ones}else{
                Indone=sample(ones,1)}
            if(length(zeros)==1){Indzero=zeros}else{
              Indzero=sample(zeros,1)}
            ####Change Beta/eta
            eta1[b,Indone]=0
            eta1[b,Indzero]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
            spot1=max(spotold)
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])
            Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
            ###Generate new vector##
            beta1[b,Indone]=0

            ##meannew,varnew##
            V1 = Sigmanew[spot2,spot2]
            V2 = as.matrix(Sigmanew[-spot2,-spot2])
            V12 = as.matrix(Sigmanew[spot2,-spot2])
            thetab=beta1[b-1,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot2])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            beta1[b,Indzero]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta1[b,Indzero],meannew,varnew))
            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta1[b-1,Indone],meanold,varold))
            ######Accept reject###
            Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

            alphab1=Liken-Likeo+dn-do
            U=log(runif(1,0,1))
            if(is.finite(alphab1)==FALSE){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{
              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}

            }

          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta1[b-1,Ind]==1){
              ##delete##
              iter[2]="delete"
              eta1[b,Ind]=0
              beta1[b,Ind]=0
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)


              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,c(includedold,(p1+1):(p1+inc))])%*%X[,c(includedold,(p1+1):(p1+inc))])

              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta1[b-1,c(includedold,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta1[b-1,Ind],meanold,varold))
              ######Accept reject###
              Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
              Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

              alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{
                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}

              }

            }else{
              ###Add###
              eta1[b,Ind]=1

              iter[2]="add"
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta1[b-1,c(includednew,(p1+1):(p1+inc))]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta1[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta1[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
              Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
              U=log(runif(1,0,1))
              if(is.finite(alphab1)==FALSE){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{
                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}
              }


            }

          }}





        ###End SVSS


        ###INCLUDED SAMPLERS

        iter[1]="Beta1"
        iter[2]="Included"

        if(sum(eta1[b,])==0){
          ##Sample Included
          Sigmanew= c*solve(t(X[,(p1+1):(p1+inc)])%*%X[,(p1+1):(p1+inc)])
          zeta1n=beta1[b,(p1+1):(p1+inc)]
          meannew=0
          varnew = sqrt(Sigmanew)
          zeta1=rnorm(1,meannew,varnew)
          dn=log(dnorm(zeta1,meannew,varnew))
          ###density old
          do=log(dnorm(zeta1n,meannew,varnew))
          beta=beta1[b,]
          beta[(p1+1):(p1+inc)]=zeta1
          Likeo=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])
          Liken=LK1(Y1,I1,X,beta,s1[b-1,],lam1[b-1,])

          alphab1m=Liken-Likeo+dn -do
          U=log(runif(1,0,1))

          if(is.finite(alphab1m)==FALSE){
            Ind1s[b]=0
          }else{
            if(U>alphab1m){
              Ind1s[b]=0
            }else{Ind1s[b]=1
            beta1[b,]=beta
            zeta1n=zeta1
            }}
          ##End Inc Sampler
        }else{
          includednew=rep(0,p1)
          for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
          includednew=includednew[includednew != 0]
          zeta1n=beta1[b,c(includednew,(p1+1):(p1+inc))]

          ###Make sigma matrices##
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])
          ####
          p=length(includednew)+inc
          ####Update All included variables

          for(k in (length(includednew)+1):(length(includednew)+inc)){
            zeta1=zeta1n
            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetano = as.matrix(zeta1[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1[k],meannew,varnew))
            ###density old
            do=log(dnorm(beta1[b,(p1+k-length(includednew))],meannew,varnew))


            ######Accept reject###
            Likeo=LK1(Y1,I1,X,c(beta1[b,1:p1],zeta1n[(length(zeta1n)-inc+1):length(zeta1n)]),s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,c(beta1[b,1:p1],zeta1[(length(zeta1n)-inc+1):length(zeta1n)]),s1[b-1,],lam1[b-1,])

            alphab1s=Liken-Likeo+dn -do
            U=log(runif(1,0,1))

            if(is.finite(alphab1s)==FALSE){
              Ins1s[b]=0

            }else{

              if(U>alphab1s){



                Ind1s[b]=0

              }else{Ind1s[b]=1
              zeta1n=zeta1
              beta1[b,]=c(beta1[b,1:p1],zeta1[(length(zeta1)-inc+1):length(zeta1)])

              }
            }

          }

          ###End included sampler###
        }



        #####Conditional Sampler for Included!###


        if(sum(eta1[b,])>0){

          iter[2]="Conditional Inclusion"
          ##Jointly Update nonzero betas
          zeta1=beta1[b,]
          zeta1=zeta1[zeta1!=0]
          zeta1n=zeta1
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])



          ###############
          ####

          for(k in 1:length(includednew)){


            V1 = Sigmanew[k,k]
            V2 = as.matrix(Sigmanew[-k,-k])
            V12 = as.matrix(Sigmanew[k,-k])
            thetab=beta1[b,c(includednew,(p1+1):(p1+inc))]
            thetano = as.matrix(thetab[-k])
            meannew = t(V12)%*%solve(V2)%*%thetano
            varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            ##################
            zeta1n[k]=rnorm(1,meannew,varnew)
            dn=log(dnorm(zeta1n[k],meannew,varnew))
            ###density old
            do=log(dnorm(zeta1[k],meannew,varnew))






            beta=beta1[b,]
            beta[c(includednew,(p1+1):(p1+inc))]=zeta1n


            Likeo=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta,s1[b-1,],lam1[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond1[b,k]=0
            }else{
              if(U>alphab1m){
                Indcond1[b,includednew[k]]=0
                zeta1n[k]=zeta1[k]
              }else{Indcond1[b,includednew[k]]=1
              beta1[b,]=beta
              zeta1[k]=zeta1n[k]
              }}

          }



          ##Jointly Update nonzero betas
          iter[2]="mixing"
          zeta1n=beta1[b,]
          Sigmanew=c*solve(t(X[,c(includednew,(p1+1):(p1+inc))])%*%X[,c(includednew,(p1+1):(p1+inc))])


          zeta1n[c(includednew,(p1+1):(p1+inc))]=rmvnorm(1,rep(0,length(includednew)+inc),Sigmanew)

          beta=beta1[b,]
          beta=beta[beta!=0]

          dn=log(dmvnorm(zeta1n[c(includednew,(p1+1):(p1+inc))],rep(0,length(includednew)+inc),Sigmanew))
          ###density old
          do=log(dmvnorm(beta,rep(0,length(includednew)+inc),Sigmanew))

          ######Accept reject###
          Likeo=LK1(Y1,I1,X,beta1[b,], s1[b-1,],lam1[b-1,])
          Liken=LK1(Y1,I1,X,zeta1n,s1[b-1,],lam1[b-1,])

          alphamix1=Liken-Likeo+dn -do
          U=log(runif(1,0,1))
          if(is.finite(alphamix1)==FALSE){
            Indmix1[b]=0
          }else{
            if(U>alphamix1){

              Indmix1[b]=0
            }else{Indmix1[b]=1
            beta1[b,]=zeta1n
            }}

        }

        S1=s1[b-1,]
        S1=S1[!is.na(S1)]

        L1=lam1[b-1,]
        L1=as.matrix(L1[!is.na(L1)])


        iter[1]="LogBH1"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
        Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


        length1=rep(0,J1+1)




        for(j in 1:length(length1)){
          length1[j]=s1[b-1,j+1]-s1[b-1,j]
        }


        if(J1<2){
          if(J1==1){
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            SigLam1=solve(diag(J1+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m1))
            SigLam1=Q1
          }
        }else{


          for(j in 2:J1){
            W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


          SigLam1=solve(diag(J1+1)-W1)%*%Q1

        }



        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam

        if(J1>0){

          Mulam1[b]=rnorm(1,(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%L1)/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1))),sqrt(Siglam1[b-1]/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1)))))


          Siglam1[b]=1/rgamma(1,a1+(J1+1)/2,b1+.5*(t(as.matrix(rep(Mulam1[b],J1+1))-L1)%*%solve(SigLam1)%*%(as.matrix(rep(Mulam1[b],J1+1))-L1)))


          ##Siglam

          iter[2]="Sigma"
        }else{



          Mulam1[b]=rnorm(1,lam1[b-1,1],sqrt(Siglam1[b-1]))


          Siglam1[b]=1/rgamma(1,a1+1/2,b1+.5*(Mulam1[b]-lam1[b-1,1])^2)



        }

        #if(is.finite(Mulam1[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam1[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}



        #lambda1
        iter[2]="lam1"
        lam1[b,]=lam1[b-1,]
        #######

        for(m in 1:(J1+1)){



          lam=lam1[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl1,cl1)


          if(J1==0){
            do=log(dnorm(lambda[m],Mulam1[b],sqrt(Siglam1[b])))
            dn=log(dnorm(lam[m],Mulam1[b],sqrt(Siglam1[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
            do=dmvnorm(lam,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
          }

          Likeo=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b-1,],lam1[b,])
          Liken=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b-1,],lam)




          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam1[b,m]=lam1[b-1,m]
            Acceptlam1[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam1[b,m]=1
              lam1[b,m]=lam[m]
            }else{Acceptlam1[b,m]=0}
          }


        }





        #####################################################
        ###################################################

        iter[1]="Haz1"
        iter[2]="Birth"

        ###Random Perturbation###
        U1=runif(1,0,1)
        #####

        s=s1[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J1max){
          Birth=runif(1,0,m1)

          s1[b,1:(J1+3)]=sort(c(s,Birth))

          for(k in 2:(J1+2)){
            if(Birth>s1[b-1,k-1] & Birth<s1[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J1+2)

          if(Ind==1 | Ind==J1+1){
            if(Ind==1){
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
            }else{
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam1[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b-1,],lam1[b,])

          if(J1>0){
            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))
          }else{
            do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))
          }

          prior=((2*J1+3)*(2*J1+2)*(Birth-s1[b-1,Ind])*(s1[b-1,Ind+1]-Birth))/((m1^2)*(s1[b-1,Ind+1]-s1[b-1,Ind]))

          G1=G1+1
          J1=J1+1

          Ln=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b,],lam)



          ##Make SigLam1



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1





            }else{

              SigLam1n=2/m1
            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

          }



          dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),Siglam1[b]*SigLam1n))




          alpha=Ln-Lo+dn-do-log(U1*(1-U1)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB1[b]=0
            s1[b,]=s1[b-1,]
            J1=J1-1
            G1=G1-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB1[b]=1
              lam1[b,1:(J1+1)]=lam
            }else{
              s1[b,]=s1[b-1,]
              IndB1[b]=0
              J1=J1-1
              G1=G1-1
            }

          }


        }else{
          s1[b,]=s1[b-1,]
          IndB1[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U1=runif(1,0,1)

        if(J1==0){
          IndD1[b]=0
          s1[b,]=s1[b-1,]
        }else{

          if(J1==1){
            Ind=2
          }else{

            Ind=sample(2:(J1+1),1)
          }


          s=s1[b,]
          s=s[-Ind]

          lam=lam1[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s1[b,Ind]-s1[b,Ind-1])*lam1[b,Ind-1]+(s1[b,Ind+1]-s1[b,Ind])*lam1[b,Ind])/(s1[b,Ind+1]-s1[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])



          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1=solve(diag(J1+1)-W1)%*%Q1

              do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))




            }else{


              do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))

            }
          }else{

            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1=solve(diag(J1+1)-W1)%*%Q1

            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))


          }
          #############################################
          #############################################

          Lo=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b,],lam1[b,])



          prior=((m1^2)*(s1[b,Ind+1]-s1[b,Ind-1]))/((2*J1+1)*(2*J1)*(s1[b,Ind]-s1[b,Ind-1])*(s1[b,Ind+1]-s1[b,Ind]))


          G1=G1-1
          J1=J1-1


          Ln=LK1L(Y1,I1,X,as.matrix(beta1[b,]), s,lam)


          ###Make siglam matrix



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)



          length1=rep(0,J1+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1


              dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))



            }else{

              SigLam1n=2/m1
              dn=log(dpois(J1,alpha1))+log(dnorm(lam,Mulam1[b],Siglam1[b]))

            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

            dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U1*(1-U1))

          if(is.nan(alpha)==TRUE){
            IndD1[b]=0
            J1=J1+1
            G1=G1+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s1[b,]=c(s,NA)
              IndD1[b]=1
              lam1[b,1:(J1+1)]=lam
              lam1[b,(J1+2):J1max]=rep(NA,J1max-J1-1)
            }else{
              IndD1[b]=0
              J1=J1+1
              G1=G1+1
            }
          }

          ####End else
        }
        ##





        split1[b]=J1

        sum1[b]=sum(eta1[b,])


      }





      ################End Samplers
      cat(c,z1a,z1b,"




          ", "



          ", "


          ", "Posterior Inclusion Probabilities after half Burnin", "


          ", "Hazard 1", "

          ", colMeans(eta1[(B*burn+1):B,])*100, "


          ", "IndEta",mean(Indeta1[(B*burn+1):B])*100,"


          ","IndMix",mean(Indmix1[(B*burn+1):B])*100,"


          ", "Included Acceptance", "

          ", "Haz1", "
          ", mean(Ind1s[(B*burn+1):B])*100, "

          ", colMeans(Indcond1[(B*burn+1):B,],na.rm=TRUE)*100,"




          ","Survival","

          ","IndDeath",mean(IndD1[(B*burn+1):B])*100,"

          ","IndBirth",mean(IndB1[(B*burn+1):B])*100,"

          ","Lambda","

          ", "Lam1",
          colMeans(Acceptlam1[(B*burn+1):B,],na.rm=TRUE)*100)




      Path1= paste0(Path,"/Ind1s.txt")

      write.table(Ind1s[(burn*B+1):B], Path1, sep="\t")




      par(mfrow=c(2,1))

      plot(1:B,sum1,type="l",xlab="",ylab="Haz: # Included", main="Traceplot: # Included")



      plot(1:B,split1,type="l",xlab="",ylab="Haz: Split #", main="Traceplot: # Split points")



    }



    ###If 0 inc



    if(inc==0){

      cat("No Variables Included")

      for(b in 2:B){


        if(b%%10000==0){cat(b, "iterations",date(), "  ")}else{
          if(b%%5000==0){cat(b, " iterations ")}}




        ###eta1,beta1
        eta1[b,]=eta1[b-1,]
        beta1[b,]=beta1[b-1,]

        if(sum(eta1[b-1,])==0|sum(eta1[b-1,])==p1){
          if(sum(eta1[b-1,])==0){
            ###Add Automatically
            Ind=sample(1:p1,1)
            eta1[b,Ind]=1
            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]
            spotnew=rep(0,length(includednew))
            for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
            spot2=max(spotnew)
            ###Make sigma matrices##
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
            ####

            meannew = 0
            varnew = sqrt(Sigmanew)
            ##################
            beta1[b,Ind]=rnorm(1,meannew,varnew)
            dn=log(dnorm(beta1[b,Ind],meannew,varnew))
            ######Accept reject###

            Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

            alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))

            if(U>alphab1){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{Indeta1[b]=1}
          }
          if(sum(eta1[b-1,])==p1){
            ###Delete Automatically
            Ind=sample(1:p1,1)
            eta1[b,Ind]=0
            beta1[b,Ind]=0
            includedold=rep(0,p1)
            for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
            includedold=includedold[includedold != 0]
            spotold=rep(0,length(includedold))
            for(k in 1:length(includedold)){if(includedold[k]==Ind){spotold[k]=k}}
            spot1=max(spotold)


            ###Make sigma matrices##
            Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

            ###Old density###
            V1 = Sigmaold[spot1,spot1]
            V2 = as.matrix(Sigmaold[-spot1,-spot1])
            V12 = as.matrix(Sigmaold[spot1,-spot1])
            thetab=beta1[b-1,includedold]
            thetano = as.matrix(thetab[-spot1])
            meanold = t(V12)%*%solve(V2)%*%thetano
            varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
            do=log(dnorm(beta1[b-1,Ind],meanold,varold))
            ######Accept reject###
            Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

            alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
            U=log(runif(1,0,1))

            if(U>alphab1){
              eta1[b,]=eta1[b-1,]
              beta1[b,]=beta1[b-1,]
              Indeta1[b]=0
            }else{Indeta1[b]=1}
          }
        }else{

          U=runif(1,0,1)

          if(U<psi){

            if(sum(eta1[b-1,])==1){

              includedold=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              ones=includedold
              zeros=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==0){zeros[k]=k}}
              zeros=zeros[zeros != 0]
              ###Sample swap indices###
              Indone=ones
              Indzero=sample(zeros,1)
              ####Change Beta/eta
              eta1[b,Indone]=0
              eta1[b,Indzero]=1
              beta1[b,Indone]=0
              ##

              Sigmaold=c*solve(t(X[,Indone])%*%X[,Indone])
              Sigmanew=c*solve(t(X[,Indzero])%*%X[,Indzero])

              meannew = 0
              varnew = sqrt(Sigmanew)
              ##################
              beta1[b,Indzero]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta1[b,Indzero],meannew,varnew))
              ###Old density###

              meanold = 0
              varold = sqrt(Sigmaold)
              do=log(dnorm(beta1[b-1,Indone],meanold,varold))
              ######Accept reject###
              Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
              Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

              alphab1=Liken-Likeo+dn-do
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }else{



              ###Swapper
              includedold=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
              includedold=includedold[includedold != 0]
              ones=includedold
              zeros=rep(0,p1)
              for(k in 1:p1){if(eta1[b-1,k]==0){zeros[k]=k}}
              zeros=zeros[zeros != 0]
              ###Sample swap indices###
              Indone=sample(ones,1)
              Indzero=sample(zeros,1)
              ####Change Beta/eta
              eta1[b,Indone]=0
              eta1[b,Indzero]=1
              includednew=rep(0,p1)
              for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotold=rep(0,length(includedold))
              for(k in 1:length(includedold)){if(Indone==includedold[k]){spotold[k]=k}}
              spot1=max(spotold)
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Indzero==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])
              Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
              ###Generate new vector##
              beta1[b,Indone]=0

              ##meannew,varnew##
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta1[b-1,includednew]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta1[b,Indzero]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta1[b,Indzero],meannew,varnew))
              ###Old density###
              V1 = Sigmaold[spot1,spot1]
              V2 = as.matrix(Sigmaold[-spot1,-spot1])
              V12 = as.matrix(Sigmaold[spot1,-spot1])
              thetab=beta1[b-1,includedold]
              thetano = as.matrix(thetab[-spot1])
              meanold = t(V12)%*%solve(V2)%*%thetano
              varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              do=log(dnorm(beta1[b-1,Indone],meanold,varold))
              ######Accept reject###
              Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
              Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

              alphab1=Liken-Likeo+dn-do
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}
            }

          }else{
            ###Add/Delete
            Ind=sample(1:p1,1)
            if(eta1[b-1,Ind]==1){
              ##delete##

              if(sum(eta1[b-1,])==1){


                eta1[b,Ind]=0
                beta1[b,Ind]=0
                includedold=rep(0,p1)
                for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
                includedold=includedold[includedold != 0]
                spotold=rep(0,length(includedold))
                for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
                spot1=max(spotold)


                ###Make sigma matrices##
                Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

                ###Old density###
                V1 = Sigmaold[spot1,spot1]
                V2 = as.matrix(Sigmaold[-spot1,-spot1])
                V12 = as.matrix(Sigmaold[spot1,-spot1])
                thetab=beta1[b-1,includedold]
                thetano = as.matrix(thetab[-spot1])
                meanold = 0
                varold = sqrt(Sigmaold)
                do=log(dnorm(beta1[b-1,Ind],meanold,varold))
                ######Accept reject###
                Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
                Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

                alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
                U=log(runif(1,0,1))

                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}



              }else{


                eta1[b,Ind]=0
                beta1[b,Ind]=0
                includedold=rep(0,p1)
                for(k in 1:p1){if(eta1[b-1,k]==1){includedold[k]=k}}
                includedold=includedold[includedold != 0]
                spotold=rep(0,length(includedold))
                for(k in 1:length(includedold)){if(Ind==includedold[k]){spotold[k]=k}}
                spot1=max(spotold)


                ###Make sigma matrices##
                Sigmaold=c*solve(t(X[,includedold])%*%X[,includedold])

                ###Old density###
                V1 = Sigmaold[spot1,spot1]
                V2 = as.matrix(Sigmaold[-spot1,-spot1])
                V12 = as.matrix(Sigmaold[spot1,-spot1])
                thetab=beta1[b-1,includedold]
                thetano = as.matrix(thetab[-spot1])
                meanold = t(V12)%*%solve(V2)%*%thetano
                varold = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
                do=log(dnorm(beta1[b-1,Ind],meanold,varold))
                ######Accept reject###
                Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
                Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

                alphab1=Liken-Likeo-do + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
                U=log(runif(1,0,1))

                if(U>alphab1){
                  eta1[b,]=eta1[b-1,]
                  beta1[b,]=beta1[b-1,]
                  Indeta1[b]=0
                }else{Indeta1[b]=1}

              }


            }else{
              ###Add###



              eta1[b,Ind]=1


              includednew=rep(0,p1)
              for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
              includednew=includednew[includednew != 0]
              spotnew=rep(0,length(includednew))
              for(k in 1:length(includednew)){if(Ind==includednew[k]){spotnew[k]=k}}
              spot2=max(spotnew)
              ###Make sigma matrices##
              Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])
              ####
              V1 = Sigmanew[spot2,spot2]
              V2 = as.matrix(Sigmanew[-spot2,-spot2])
              V12 = as.matrix(Sigmanew[spot2,-spot2])
              thetab=beta1[b-1,includednew]
              thetano = as.matrix(thetab[-spot2])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              beta1[b,Ind]=rnorm(1,meannew,varnew)
              dn=log(dnorm(beta1[b,Ind],meannew,varnew))
              ######Accept reject###
              Likeo=LK1(Y1,I1,X,beta1[b-1,],s1[b-1,],lam1[b-1,])
              Liken=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])

              alphab1=Liken-Likeo+dn + log(beta(sum(eta1[b,])+z1a,p1-sum(eta1[b,])+z1b)) - log(beta(sum(eta1[b-1,])+z1a,p1-sum(eta1[b-1,])+z1b))
              U=log(runif(1,0,1))

              if(U>alphab1){
                eta1[b,]=eta1[b-1,]
                beta1[b,]=beta1[b-1,]
                Indeta1[b]=0
              }else{Indeta1[b]=1}



            }

          }

        }

        ##End Eta Beta

        includednew=rep(0,p1)
        for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
        includednew=includednew[includednew != 0]

        if(sum(eta1[b,])>0){


          if(sum(eta1[b,])==1){

            iter[2]="Conditional Inclusion"


            includednew=rep(0,p1)
            for(k in 1:p1){if(eta1[b,k]==1){includednew[k]=k}}
            includednew=includednew[includednew != 0]

            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])



            meannew = 0
            varnew = sqrt(Sigmanew)

            beta=beta1[b,]

            ##################
            beta[includednew]=rnorm(1,meannew,varnew)

            dn=log(dnorm(beta[includednew],meannew,varnew))
            ###density old
            do=log(dnorm(beta1[b,includednew],meannew,varnew))







            Likeo=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,beta,s1[b-1,],lam1[b-1,])

            alphab1m=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphab1m)==FALSE){
              Indcond1[b,includednew]=0
            }else{
              if(U>alphab1m){
                Indcond1[b,includednew]=0
              }else{Indcond1[b,includednew]=1
              beta1[b,]=beta
              }}

          }else{

            iter[2]="Conditional Inclusion"
            ##Jointly Update nonzero betas
            zeta1=beta1[b,]
            zeta1=zeta1[zeta1!=0]
            zeta1n=zeta1
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])



            ###############
            ####

            for(k in 1:length(includednew)){


              V1 = Sigmanew[k,k]
              V2 = as.matrix(Sigmanew[-k,-k])
              V12 = as.matrix(Sigmanew[k,-k])
              thetab=beta1[b,includednew]
              thetano = as.matrix(thetab[-k])
              meannew = t(V12)%*%solve(V2)%*%thetano
              varnew = sqrt(V1 - t(V12)%*%solve(V2)%*%V12)
              ##################
              zeta1n[k]=rnorm(1,meannew,varnew)
              dn=log(dnorm(zeta1n[k],meannew,varnew))
              ###density old
              do=log(dnorm(zeta1[k],meannew,varnew))






              beta=beta1[b,]
              beta[includednew]=zeta1n


              Likeo=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])
              Liken=LK1(Y1,I1,X,beta,s1[b-1,],lam1[b-1,])

              alphab1m=Liken-Likeo+dn -do
              U=log(runif(1,0,1))
              if(is.finite(alphab1m)==FALSE){
                Indcond1[b,includednew[k]]=0
              }else{
                if(U>alphab1m){
                  Indcond1[b,includednew[k]]=0
                  zeta1n[k]=zeta1[k]
                }else{Indcond1[b,includednew[k]]=1
                beta1[b,]=beta
                zeta1[k]=zeta1n[k]
                }}

            }



            ##Jointly Update nonzero betas
            iter[2]="mixing"
            zeta1n=beta1[b,]
            Sigmanew=c*solve(t(X[,includednew])%*%X[,includednew])


            zeta1n[includednew]=rmvnorm(1,rep(0,length(includednew)),Sigmanew)

            beta=beta1[b,]
            beta=beta[beta!=0]

            dn=log(dmvnorm(zeta1n[includednew],rep(0,length(includednew)),Sigmanew))
            ###density old
            do=log(dmvnorm(beta,rep(0,length(includednew)),Sigmanew))

            ######Accept reject###
            Likeo=LK1(Y1,I1,X,beta1[b,],s1[b-1,],lam1[b-1,])
            Liken=LK1(Y1,I1,X,zeta1n,s1[b-1,],lam1[b-1,])

            alphamix1=Liken-Likeo+dn -do
            U=log(runif(1,0,1))
            if(is.finite(alphamix1)==FALSE){
              Indmix1[b]=0
            }else{
              if(U>alphamix1){

                Indmix1[b]=0
              }else{Indmix1[b]=1
              beta1[b,]=zeta1n
              }}

          }



        }





        S1=s1[b-1,]
        S1=S1[!is.na(S1)]


        L1=lam1[b-1,]
        L1=as.matrix(L1[!is.na(L1)])


        ############################################
        #####Start LogBH Samplers###################
        ############################################
        ####Lam1####

        iter[1]="LogBH1"
        iter[2]="matrixsetup"

        W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
        Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


        length1=rep(0,J1+1)




        for(j in 1:length(length1)){
          length1[j]=s1[b-1,j+1]-s1[b-1,j]
        }


        if(J1<2){
          if(J1==1){
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            SigLam1=solve(diag(J1+1)-W1)%*%Q1





          }else{

            Q1=as.matrix(2/(m1))
            SigLam1=Q1
          }
        }else{


          for(j in 2:J1){
            W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
            Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
          }


          Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


          SigLam1=solve(diag(J1+1)-W1)%*%Q1

        }



        iter[2]="Mu"
        ##Lambda1 Hierarchical Sampler
        ##Mulam

        if(J1>0){

          Mulam1[b]=rnorm(1,(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%L1)/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1))),sqrt(Siglam1[b-1]/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1)))))


          Siglam1[b]=1/rgamma(1,a1+(J1+1)/2,b1+.5*(t(as.matrix(rep(Mulam1[b],J1+1))-L1)%*%solve(SigLam1)%*%(as.matrix(rep(Mulam1[b],J1+1))-L1)))


          ##Siglam

          iter[2]="Sigma"
        }else{



          Mulam1[b]=rnorm(1,lam1[b-1,1],sqrt(Siglam1[b-1]))


          Siglam1[b]=1/rgamma(1,a1+1/2,b1+.5*(Mulam1[b]-lam1[b-1,1])^2)



        }

        #if(is.finite(Mulam1[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}
        #if(is.finite(Siglam1[b])==FALSE){stop("Adjust Hierarchical Hyper-Parameters")}



        #lambda1
        iter[2]="lam1"
        lam1[b,]=lam1[b-1,]
        #######

        for(m in 1:(J1+1)){



          lam=lam1[b,]
          lam=lam[is.na(lam)==FALSE]
          lambda=lam

          lam[m]=lambda[m]+runif(1,-cl1,cl1)


          if(J1==0){
            do=log(dnorm(lambda[m],Mulam1[b],sqrt(Siglam1[b])))
            dn=log(dnorm(lam[m],Mulam1[b],sqrt(Siglam1[b])))
          }else{


            #do=-(t(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lambda)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])

            #dn=-(t(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1)))%*%solve(SigLam1)%*%(as.matrix(lam)-as.matrix(rep(Mulam1[b],J1+1))))/(2*Siglam1[b])


            do=dmvnorm(lambda,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
            do=dmvnorm(lam,rep(Mulam1[b],J1+1),Siglam1[b]*SigLam1)
          }

          Likeo=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b-1,],lam1[b,])
          Liken=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b-1,],lam)


          U=log(runif(1,0,1))
          alphalam=Liken-Likeo+dn-do

          if(is.nan(alphalam)==TRUE){
            lam1[b,m]=lam1[b-1,m]
            Acceptlam1[b,m]=0
          }else{

            if(U<alphalam){
              Acceptlam1[b,m]=1
              lam1[b,m]=lam[m]
            }else{Acceptlam1[b,m]=0}
          }


        }






        #############################################
        ###################################################

        iter[1]="Haz1"
        iter[2]="Birth"

        ###Random Perturbation###
        U1=runif(1,0,1)
        #####

        s=s1[b-1,]
        s=s[!is.na(s)]

        if(length(s)<J1max){
          Birth=runif(1,0,m1)

          s1[b,1:(J1+3)]=sort(c(s,Birth))

          for(k in 2:(J1+2)){
            if(Birth>s1[b-1,k-1] & Birth<s1[b-1,k]){
              Ind=k-1
            }
          }

          lam=rep(0,J1+2)

          if(Ind==1 | Ind==J1+1){
            if(Ind==1){
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
            }else{
              lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
              lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            }
          }else{
            lam[Ind]=lam1[b,Ind] - ((s1[b-1,Ind+1]-Birth)/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[Ind+1]=lam1[b,Ind] + ((Birth-s1[b-1,Ind])/(s1[b-1,Ind+1]-s1[b-1,Ind]))*log((1-U1)/U1)
            lam[1:(Ind-1)]=lam1[b,1:(Ind-1)]
            lam[(Ind+2):length(lam)]=lam1[b,(Ind+1):(J1+1)]
          }

          lam=lam[!is.na(lam)]

          lambda=lam1[b,]
          lambda=lambda[!is.na(lambda)]

          Lo=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b-1,],lam1[b,])



          if(J1>0){
            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))
          }else{
            do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))
          }

          prior=((2*J1+3)*(2*J1+2)*(Birth-s1[b-1,Ind])*(s1[b-1,Ind+1]-Birth))/((m1^2)*(s1[b-1,Ind+1]-s1[b-1,Ind]))

          G1=G1+1
          J1=J1+1

          Ln=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b,],lam)



          ##Make SigLam1



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1





            }else{

              SigLam1n=2/m1
            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

          }



          dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),Siglam1[b]*SigLam1n))




          alpha=Ln-Lo+dn-do-log(U1*(1-U1)) + log(prior)

          if(is.nan(alpha)==TRUE){
            IndB1[b]=0
            s1[b,]=s1[b-1,]
            J1=J1-1
            G1=G1-1
          }else{

            U=log(runif(1,0,1))

            if(U<alpha){
              IndB1[b]=1
              lam1[b,1:(J1+1)]=lam
            }else{
              s1[b,]=s1[b-1,]
              IndB1[b]=0
              J1=J1-1
              G1=G1-1
            }

          }


        }else{
          s1[b,]=s1[b-1,]
          IndB1[b]=0
        }


        #########################################################
        ###################Death Sampler#########################
        ##########################################################
        iter[2]="Death"

        U1=runif(1,0,1)

        if(J1==0){
          IndD1[b]=0
          s1[b,]=s1[b-1,]
        }else{

          if(J1==1){
            Ind=2
          }else{

            Ind=sample(2:(J1+1),1)
          }


          s=s1[b,]
          s=s[-Ind]

          lam=lam1[b,]
          lambda=lam[!is.na(lam)]

          lam=lam[!is.na(lam)]
          lam=lam[-Ind]

          lam[Ind-1]=((s1[b,Ind]-s1[b,Ind-1])*lam1[b,Ind-1]+(s1[b,Ind+1]-s1[b,Ind])*lam1[b,Ind])/(s1[b,Ind+1]-s1[b,Ind-1])



          #############################################
          ####Sets up SigLam1 matrix for old density###
          #############################################


          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


          length1=diff(s1[b,])



          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1=solve(diag(J1+1)-W1)%*%Q1

              do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))




            }else{


              do=log(dpois(J1,alpha1))+log(dnorm(lambda,Mulam1[b],Siglam1[b]))

            }
          }else{

            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1=solve(diag(J1+1)-W1)%*%Q1

            do=log(dpois(J1,alpha1))+log(dmvnorm(lambda,rep(Mulam1[b],length(lambda)),SigLam1*Siglam1[b]))


          }
          #############################################
          #############################################

          Lo=LK1L(Y1,I1,X,as.matrix(beta1[b,]),s1[b,],lam1[b,])



          prior=((m1^2)*(s1[b,Ind+1]-s1[b,Ind-1]))/((2*J1+1)*(2*J1)*(s1[b,Ind]-s1[b,Ind-1])*(s1[b,Ind+1]-s1[b,Ind]))


          G1=G1-1
          J1=J1-1


          Ln=LK1L(Y1,I1,X,as.matrix(beta1[b,]), s,lam)


          ###Make siglam matrix



          W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
          Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)



          length1=rep(0,J1+1)




          for(j in 1:length(length1)){
            length1[j]=s[j+1]-s[j]
          }


          if(J1<2){
            if(J1==1){
              W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
              W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
              Q1[1,1]=2/(2*length1[1]+length1[2])
              Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
              SigLam1n=solve(diag(J1+1)-W1)%*%Q1


              dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))



            }else{

              SigLam1n=2/m1
              dn=log(dpois(J1,alpha1))+log(dnorm(lam,Mulam1[b],Siglam1[b]))

            }
          }else{


            for(j in 2:J1){
              W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
              Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
            }


            Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
            Q1[1,1]=2/(2*length1[1]+length1[2])
            W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
            W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


            SigLam1n=solve(diag(J1+1)-W1)%*%Q1

            dn=log(dpois(J1,alpha1))+log(dmvnorm(lam,rep(Mulam1[b],length(lam)),SigLam1n*Siglam1[b]))


          }
          ####




          alpha=Ln-Lo+dn-do+log(prior)+log(U1*(1-U1))

          if(is.nan(alpha)==TRUE){
            IndD1[b]=0
            J1=J1+1
            G1=G1+1
          }else{

            U=log(runif(1,0,1))

            iter[2]="AcceptRejDeath"

            if(U<alpha){
              s1[b,]=c(s,NA)
              IndD1[b]=1
              lam1[b,1:(J1+1)]=lam
              lam1[b,(J1+2):J1max]=rep(NA,J1max-J1-1)
            }else{
              IndD1[b]=0
              J1=J1+1
              G1=G1+1
            }
          }

          ####End else
        }
        ##






        #######################
        #####End of Death sampler
        ######################



        split1[b]=J1

        ##
        sum1[b]=sum(eta1[b,])



        ##End Sampler
      }
      ###End of Sampler




      ################End Samplers
      cat(c,z1a,z1b,"




          ", "



          ", "


          ", "Posterior Inclusion Probabilities after half Burnin", "


          ", "Hazard", "

          ", colMeans(eta1[(B*burn+1):B,])*100, "


          ", "IndEta",mean(Indeta1[(B*burn+1):B])*100,"


          ","IndMix",mean(Indmix1[(B*burn+1):B])*100,"


          ", "Included Acceptance", "

          ", "Hazard", "
          ", "

          ", colMeans(Indcond1[(B*burn+1):B,],na.rm=TRUE)*100,"



          ","Survival","

          ","IndDeath",mean(IndD1[(B*burn+1):B])*100,"

          ","IndBirth",mean(IndB1[(B*burn+1):B])*100,"

          ","Lambda","

          ", "Lam1",
          colMeans(Acceptlam1[(B*burn+1):B,],na.rm=TRUE)*100 )




      ##End


    }





    ###Return Values



    Path1= paste0(Path,"/beta1.txt")

    write.table(beta1[(burn*B+1):B,], Path1, sep="\t")



    Path1= paste0(Path,"/eta1.txt")

    write.table(eta1[(burn*B+1):B,], Path1, sep="\t")





    Path1= paste0(Path,"/lam1.txt")

    write.table(lam1[(burn*B+1):B,], Path1, sep="\t")



    Path1= paste0(Path,"/s1.txt")

    write.table(s1[(burn*B+1):B,], Path1, sep="\t")



    Path1= paste0(Path,"/sum1.txt")

    write.table(sum1[(burn*B+1):B], Path1, sep="\t")



    Path1= paste0(Path,"/split1.txt")

    write.table(split1[(burn*B+1):B], Path1, sep="\t")



    Path1= paste0(Path,"/siglam1.txt")

    write.table(Siglam1[(burn*B+1):B], Path1, sep="\t")


    Path1= paste0(Path,"/mulam1.txt")

    write.table(Mulam1[(burn*B+1):B], Path1, sep="\t")



    Path1= paste0(Path,"/Indeta1.txt")

    write.table(Indeta1[(burn*B+1):B], Path1, sep="\t")



    Path1= paste0(Path,"/IndD1.txt")

    write.table(IndD1[(burn*B+1):B], Path1, sep="\t")

    Path1= paste0(Path,"/IndB1.txt")

    write.table(IndB1[(burn*B+1):B], Path1, sep="\t")


    Path1= paste0(Path,"/Acceptlam1.txt")

    write.table(Acceptlam1[(burn*B+1):B,], Path1, sep="\t")








    Path1= paste0(Path,"/Indmix1.txt")


    write.table(Indmix1[(burn*B+1):B], Path1, sep="\t")









    par(mfrow=c(2,1))

    plot(1:B,sum1,type="l",xlab="",ylab=" # Included", main="Traceplot: # Included")



    plot(1:B,split1,type="l",xlab="",ylab="Hazard: Split #", main="Traceplot: # Split points")


  }
}







