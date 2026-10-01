#' Refine cell-type purity estimates using a deconvolution EM algorithm
#'
#' Refines prior sample-level purity estimates using bulk transcriptomic data.
#' The function models each gene as a mixture of two latent expression
#' components and alternates between estimating gene-specific distribution
#' parameters and updating the sample-specific purity values.
#'
#' Gene-level parameter estimation is performed using \code{foreach} and
#' \code{\%dopar\%}. When a parallel backend has been registered, genes are
#' processed in parallel; otherwise, \code{foreach} may run sequentially.
#'
#' @param exprM A numeric N by G matrix of bulk expression values, where rows
#'   represent samples and columns represent genes. Sample and gene identifiers
#'   may be supplied as row and column names, respectively.
#' @param purity A numeric vector of length N containing prior purity or
#'   cell-type proportion estimates for the samples. Values should lie strictly
#'   between zero and one and must follow the row ordering of \code{exprM}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Epurity}{A numeric vector of length N containing the refined
#'     sample-level purity estimates.}
#'   \item{paraM}{A G by 4 numeric matrix of fitted gene-specific parameters.
#'     The columns \code{uy} and \code{uz} contain the estimated component
#'     means, while \code{vy} and \code{vz} contain the estimated component
#'     variances.}
#'   \item{convergence}{A character string indicating whether the outer
#'     expectation-maximization algorithm converged.}
#'   \item{time}{The processing time reported by \code{proc.time}.}
#' }
#'
#' @details
#' The columns of \code{exprM} are standardized before model fitting. For each
#' gene, the algorithm estimates the means and variances of two latent
#' expression components using an inner expectation-maximization procedure.
#' Given these gene-level parameters, sample-specific purity estimates and the
#' concentration parameter of their beta distribution are updated by numerical
#' optimization. These steps are repeated until the purity estimates satisfy
#' the convergence criterion or the maximum number of iterations is reached.
#'
#' The function requires the internal helper functions \code{Estep.1d()} and
#' \code{Mstep.1d()}.
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom stats dbeta optim var
#' @export
deNet_purity<-function(exprM,purity)
{


  t<-proc.time()

  options(warn=-1)
  convergecutoff=0.1 ### convergence criterion for EM
  exprM<-t(scale((exprM)))
  geneNames<-rownames(exprM);
  curgeneNames<-rownames(exprM);


  for(j in 1:1000)  # --- over EM iterations
    {

      if(j==1){
        cur.P<-purity;
        cur.b<-mean(cur.P)*(1-mean(cur.P))/var(cur.P)-1
      }else{
        # print(paste("Epurity",j,":",cor(purity,cur.P)))
      }

      ## -- Estimate u and v (mean and variance of latent variables, i.e., uy uz vy vz)

      paraM<-foreach(pi= 1:nrow(exprM))%dopar% { # -- for each gene

          x.v=exprM[pi,]
          for(i in 1:10000) # -- (A Loop) over EM iterations
          {
            if(i==1){
              rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)
            }
            r<-rnew;
            temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,cur.P)
            rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
            if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01))
            {break()}
          }
          c(r$uy,r$uz,r$vy,r$vz)
      }

      paraM<-do.call("rbind",paraM)
      colnames(paraM)<-c("uy","uz","vy","vz")
      # print("EM1d")
      before.P<-cur.P;
      ## -- Estimate purity (P) and delta parameter in beta distribution
      for(i in 1:10000) {

        new.P<-sapply(1:ncol(exprM),function(ni){
            f<-function(x) {
              V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
              U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
              0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)- log(dbeta(purity[ni],x*cur.b,(1-x)*cur.b)) # purity
            }
            optim(cur.P[ni],f,lower=0,upper=1,method="Brent")$par})
        f<-function(x) -sum(log(dbeta(purity,new.P*x,(1-new.P)*x)))
        new.b<-optim(cur.b,f,lower=0)$par


        if(all(abs(c(new.P-cur.P))<0.1))
        {#print(paste("converge in step",i));
          break()}
        cur.b<-new.b;
        cur.P<-new.P
        # print(cur.b)
        # print(mean(abs(cur.P-purity)))
      }
      # print(paste("cur.b",cur.b))


      if(all((before.P-cur.P)<0.1))break();
      print(mean(abs(cur.P-before.P)))
    }

    if(j==1000){ stop("EM algorithm did not converge"); convergence="P not converge"}else{
      print("EM algorithm converged"); convergence="P converged"}

    Epurity<-cur.P

  return(list(Epurity,paraM,convergence,time))

}





## -- function to be used
Estep.1d<-function(uy,uz,vy,vz, c,X,A)
{
  P<-c*A
  Ux<-P*uy+(1-P)*uz;
  Vx<-P^2*vy+(1-P)^2*vz
  Sxy<-P*vy
  Ey<-uy+(X-Ux)*Sxy/Vx
  Vy<-vy-Sxy^2/Vx

  Ey2<-Vy+Ey^2

  Sxz<-(1-P)*vz
  Ez<-uz+(X-Ux)*Sxz/Vx
  Vz<-vz-Sxz^2/Vx
  Ez2<-Vz+Ez^2

  list(Ey=Ey,Ey2=Ey2, Ez=Ez, Ez2=Ez2)
}

## -- function to be used
Mstep.1d<-function(Ey,Ey2,Ez,Ez2)
{
  uy<-mean(Ey)
  uz<-mean(Ez)
  vy<-mean(Ey2)-uy^2
  vz<-mean(Ez2)-uz^2
  list(uy=uy,uz=uz,vy=vy,vz=vz)
}

## -- function to be used
calculateLL.1d<-function(uy,uz,vy,vz,c,X,A)
{
  P<-c*A
  Ux<-P*uy+(1-P)*uz;
  Vx<-P^2*vy+(1-P)^2*vz
  sum(log( sapply(1:length(X),function(i)dnorm(X[i],Ux[i],sqrt(Vx[i])))))
}
