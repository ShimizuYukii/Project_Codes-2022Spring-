library('forecast')

bf <-
  function(y, case=c("poly", "categ", "fourier", "pcont", "pdisc"), degree=1, nslices=1, scale=FALSE)
  {
    Indicator<-function(x, H) return(ifelse((x %in% H), 1, 0))
    
    case <- match.arg(case); nobs=length(y); 
    
    if (case=="categ")
    {
      bins.y<-unique(sort(y)); 
      r<- length(unique(sort(y)))-1; fy<-array(rep(0), c(r, nobs));
      for (i in 1:r){ fy[i,]<-sapply(y, function(x) (x==bins.y[i]))} 
    }
    else if (case=="fourier")
    {
      fy<-array(rep(0), c(2*degree, nobs)); 
      for(i in 1:degree)
      {
        fy[2*i-1, 1:nobs]<- cos(2*pi*y*i);
        fy[2*i, 1:nobs]<-  sin(2*pi*y*i); 
      }
    }
    else if (case=="poly") 
    {
      if (degree==0) stop("This case is not defined");
      fy <- array(rep(0), c(degree, nobs));
      for (k in 1:degree) fy[k, ] <- y^k; 
    }
    else if (case=="pdisc")
    {
      if ((nslices==0) | (nslices==1)){message("The minimum number of slices is 2"); nslices=2;}
      r <- (degree + 1) * nslices - 1; 
      fy <- array(rep(0), c(r, nobs));
      slicing <- ldr.slices(y,nslices);
      bins.y <- slicing$bins;
      
      if (degree==0)	# Piecewise constant discontinuous
      {
        for(i in 1:r) fy[i,] <- Indicator(y, bins.y[[i]]);
      }	
      else if (degree==1) # Piecewise linear discontinuous
      { 
        for(i in 1:(nslices-1))
        {
          fy[2*i-1, ] <- Indicator(y, bins.y[[i]]);
          fy[2*i, ] <- Indicator(y, bins.y[[i]])*(y-bins.y[[i]][1]);
        }
        fy[2*nslices-1, ] <- Indicator(y, bins.y[[nslices]])*(y-bins.y[[nslices]][1]);
      } 
      else if (degree==2) # Piecewise quadratic discontinuous
      { 
        for(i in 1:(nslices-1))
        {
          fy[3*(i-1)+1, ] <- Indicator(y, bins.y[[i]]);
          fy[3*(i-1)+2, ] <- Indicator(y, bins.y[[i]])*(y-bins.y[[i]][1]);
          fy[3*(i-1)+3, ] <- Indicator(y, bins.y[[i]])*(y-bins.y[[i]][1])**2;
        }
        fy[3*nslices-2, ] <- Indicator(y, bins.y[[nslices]])*(y-bins.y[[nslices]][1]);
        fy[3*nslices-1, ] <- Indicator(y, bins.y[[nslices]])*(y-bins.y[[nslices]][1])**2;
      }
      else if (degree==3)# Piecewise cubic discontinuous
      { 
        for(i in 1:(nslices-1))
        {
          fy[4*(i-1)+1, ] <- Indicator(y, bins.y[[i]]);
          fy[4*(i-1)+2, ] <- Indicator(y, bins.y[[i]])*(y-bins.y[[i]][1]);
          fy[4*(i-1)+3, ] <- Indicator(y, bins.y[[i]])*(y-bins.y[[i]][1])**2;
          fy[4*(i-1)+4, ] <- Indicator(y, bins.y[[i]])*(y-bins.y[[i]][1])**3;
        }
        fy[4*nslices-3, ] <- Indicator(y, bins.y[[nslices]])*(y-bins.y[[nslices]][1]);
        fy[4*nslices-2, ] <- Indicator(y, bins.y[[nslices]])*(y-bins.y[[nslices]][1])**2;
        fy[4*nslices-1, ] <- Indicator(y, bins.y[[nslices]])*(y-bins.y[[nslices]][1])**3;
      } 
    } 
    else if (case=="pcont")
    {
      if ((nslices==0) | (nslices==1)){message("The minimum number of slices is 2"); nslices=2;}
      if (degree==0) stop("Piecewise Constant Continuous is not defined.");
      
      r <- nslices*degree+1; 
      fy <- array(rep(0), c(r, nobs));
      slicing <- ldr.slices(y, nslices);
      bins.y <- slicing$bins;
      
      if (degree==1)# Piecewise linear continuous
      {  
        fy[1,] <- Indicator(y, bins.y[[1]]);
        if (r>1) for(i in 1:nslices) fy[i+1,] <- Indicator(y, bins.y[[i]])*(y - bins.y[[i]][1]);
      } 
      else if (degree==2)# Piecewise quadratic continuous
      { 
        fy[1,] <- Indicator(y, bins.y[[1]]);
        for(i in 1:nslices)
        {
          fy[2*i,] <- Indicator(y, bins.y[[i]])*(y - bins.y[[i]][1]);
          fy[2*i+1,] <- Indicator(y, bins.y[[i]])*(y - bins.y[[i]][1])**2;
        }
      } 
      else if (degree==3)# Piecewise cubic continuous
      { 
        fy[1,] <- Indicator(y, bins.y[[1]]);
        
        for(i in 1:nslices)
        {
          fy[3*i-1,] <- Indicator(y, bins.y[[i]])*(y - bins.y[[i]][1]);
          fy[3*i,] <- Indicator(y, bins.y[[i]])*(y - bins.y[[i]][1])**2;
          fy[3*i+1,] <- Indicator(y, bins.y[[i]])*(y - bins.y[[i]][1])**3;
        }
      } 	
    }
    return( scale(t(Re(fy)), center=TRUE, scale=scale))
  }

orthonorm <- function(u) 
{
  if (is.null(u)) return(NULL)
  if (!(is.matrix(u))) u <- as.matrix(u);
  dd <- dim(u); n <- dd[1]; p <-dd[2];
  
  if (prod(abs(La.svd(u)$d) > 1e-08) == 0) stop("collinears vectors in orthonorm")
  if (n < p)
  {
    warning("There are too much vectors to orthonormalize in orthonorm.")
    u <- as.matrix(u[, 1:p])
    n <- p
  }
  v <- u;
  if (p > 1)
  {
    for (i in 2:p)
    {
      coef.proj <- c(crossprod(u[, i], v[, 1:(i - 1)]))/diag(crossprod(v[, 1:(i - 1)]));
      v[, i] <- u[, i] - matrix(v[, 1:(i - 1)], nrow = n) %*% matrix(coef.proj, nrow = i - 1)
    }
  }
  coef.proj <- 1/sqrt(diag(crossprod(v)))
  return(t(t(v) * coef.proj))
}

pfc <-
function(X, y, fy=NULL, numdir=NULL, structure=c("iso", "aniso", "unstr", "unstr2"), eps_aniso=1e-3, numdir.test=FALSE,...)
{
	"%^%"<-function(M, pow) 
	{ 
		if (prod(dim(M)==list(1,1))) return( as.matrix(M^pow) )
		eigenM = eigen(M) 
		return(eigenM$vectors%*%diag(c(eigenM$values)^pow)%*%t(eigenM$vectors))  
	}
  
	Trace<-function(X)
	{	
		if (!is.matrix(X)) stop("Argument to Trace is not a matrix in pfc");
		return(sum(diag(X))) 
	}
	onepfc = function(X, y, fy)
	{
		# X is univariate predictor
		nobs <- length(X); r <- dim(fy)[2]

		P_F <- fy%*%solve(t(fy)%*%fy)%*%t(fy)
		Xc <- scale(X, TRUE, FALSE)
		Sigmahat_fit <- (1/nobs)*t(Xc)%*%P_F%*%(Xc)
		ev.fit <- eigen(Sigmahat_fit)

		temp.dat<-data.frame(cbind(X, fy)); xnam<-paste("xx", 1:r, sep="");
		names(temp.dat)<-c("yy", xnam);
		fm.lm<- as.formula( paste("yy ~ ", paste(xnam, collapse= "+")));
		summary.fm <- summary(lm(fm.lm, data=temp.dat));

		Betahat <- matrix(summary.fm$coefficients[2:(r+1),1], ncol=r);
		Gammahat <- matrix(1, ncol=1, nrow=1); 
		Deltahat <- matrix(summary.fm$sigma^2, ncol=1, nrow=1);
		Muhat <- matrix(summary.fm$coefficients[1,1], ncol=1);

		loglik <- - 0.5*n*(1+log(2*pi*summary.fm$sigma^2));
		numpar <- p  + dim(fy)[2] + 1;
		aic <- -2*loglik + 2*numpar;
		bic <- -2*loglik + log(n)*numpar;

		ans <- list(R=X, Muhat=Muhat, Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, 
				loglik=loglik, aic=aic, bic=bic, numpar=numpar, numdir=1, model="pfc", 
				call=match.call(), structure="iso", y=y, fy=fy, Xc=Xc, numdir.test=numdir.test);

		class(ans)<- "pfc";

		return(ans)
	}
	
	if (is.null(fy)) {fy <- scale(y, TRUE, TRUE); numdir <- 1}
	r <- dim(fy)[2]; X <- as.matrix(X)
	op <- dim(X); n <- op[1]; p <- op[2] 
	eff.numdir <- min(numdir, r, p)	

	vnames <- dimnames(X)[[2]]
	if (is.null(vnames)) vnames <- paste("X", 1:p, sep="")
	if (p==1) return(onepfc(X=X, y=y, fy=fy))
	
	Muhat <- apply(X, 2, mean)
	Xc <-  scale(X, TRUE, FALSE)

	P_F <- fy%*%solve(t(fy)%*%fy)%*%t(fy)
  
	Sigmahat <- cov(X)
	Sigmahat_fit <- cov(P_F%*%X)
	Sigmahat_res <- Sigmahat - Sigmahat_fit

	if (structure=="iso")
	{
		iso <- function(i)
		{
			ev <- eigen(Sigmahat);	
			ev.fit <- eigen(Sigmahat_fit); 
			all_evalues <-ev.fit$values
			evalues <- all_evalues[1:i]
			sigma2hat <- Re(sum(ev$values)/p); 

			Gammahat <- Re(matrix(ev.fit$vectors[,1:i], ncol=i));
			dimnames(Gammahat) <- list(vnames, paste("Dir", 1:i, sep=""))
			Betahat <-Re(t(Gammahat)%*%t(Xc)%*%fy%*%solve(t(fy)%*%fy));
			sigma2hat <- Re((sum(ev$values)-sum(evalues))/p); 
			Deltahat <- sigma2hat*diag(1, p); 
			dimnames(Deltahat) <- list(vnames, vnames)

			loglik <- - 0.5*n*p*(1+log(2*pi*sigma2hat));
			numpar <- p + (p-i)*i + i*dim(fy)[2] + 1;
			aic <- -2*loglik + 2*numpar;
			bic <- -2*loglik + log(n)*numpar;

			return(list(Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, 
					evalues=evalues, loglik=loglik, aic=aic, bic=bic, numpar=numpar));
		}

		if (identical(numdir.test, FALSE))
		{
			out <- iso(eff.numdir);

			ans <- list(R=X%*%orthonorm(out$Gammahat), Muhat=Muhat, Betahat=out$Betahat, Gammahat=out$Gammahat, 
					Deltahat=out$Deltahat, loglik=out$loglik, aic=out$aic, bic=out$bic, numpar=out$numpar, 
					numdir=eff.numdir, evalues=out$evalues, structure="iso", y=y, fy=fy,  
					Xc=Xc, call=match.call(expand.dots=TRUE), numdir.test=numdir.test);

			class(ans) <- "pfc";	
			return(ans);		
		}

		if (identical(numdir.test, TRUE))
		{
			aic <- bic <- numpar <- loglik <- vector(length=eff.numdir+1);
			Betahat <- Deltahat <- Gammahat <-vector("list"); 

			# No fitting values (eff.numdir=0)
			ev <- eigen(Sigmahat); 
			sigma2hat <- sum(ev$values)/p; 
			loglik[1] <- - 0.5*n*p*(1+log(2*pi*sigma2hat));
			numpar[1] <- p + 1;
			aic[1] <- -2*loglik[1] + 2*numpar[1];
			bic[1] <- -2*loglik[1] + log(n)*numpar[1];
		
			for (i in 1:eff.numdir)
			{
				fit <- iso(i);
				Betahat[[i]] <-fit$Betahat; 
				Gammahat[[i]] <-fit$Gammahat; 
				Deltahat[[i]] <- fit$Deltahat;
				loglik[i+1] <- fit$loglik; 
				numpar[i+1] <- fit$numpar;
				aic[i+1] <- fit$aic; 
				bic[i+1] <- fit$bic;	
			}
			ans <- list(R=X%*%orthonorm(Gammahat[[eff.numdir]]), Muhat=Muhat, Betahat=Betahat, Gammahat=Gammahat, 
					Deltahat=Deltahat, loglik=loglik, aic=aic, bic=bic, numpar=numpar, 
					numdir=eff.numdir, model="pfc", evalues=fit$evalues, structure="iso", 
					y=y, fy=fy,  Xc=Xc, call=match.call(), numdir.test=numdir.test);

			class(ans)<- "pfc";
			return(ans)
		} 
	}

	if (structure=="aniso")
	{
		aniso = function(X, y, fy, d, eps_aniso=1e-3, numdir.test)
		{
			vnames <- dimnames(X)[[2]]
			if (is.null(vnames)) vnames <- paste("X", 1:ncol(X), sep="")
			op <- dim(X); n <- op[1]; p <- op[2]

			# Initial Step
			fit <- pfc(X=X, y=y, fy=fy, numdir=d, structure="iso", numdir.test=numdir.test)

			if (identical(numdir.test, FALSE))
			{
				Betahatx <- fit$Betahat; Gammahatx <- fit$Gammahat
				Xc <- scale(X, TRUE, FALSE) - fy%*%t(Gammahatx%*%Betahatx)
				deltahat <- diag(cov(Xc))

				repeat
				{
					Xnew = X%*%((1/sqrt(deltahat))*diag(p))
					fit <- pfc(X=Xnew, y=y, fy=fy, numdir=d, structure="iso", numdir.test=FALSE)
					Betahatx <- fit$Betahat; Gammahatx <- (diag(p)*sqrt(deltahat))%*%fit$Gammahat 
					Xc <- scale(X, TRUE, FALSE) - fy%*%t(Gammahatx%*%Betahatx)
					deltahat0 <- diag(t(Xc)%*%(Xc)/n)
					if (sum(abs(deltahat-deltahat0)) < eps_aniso) break
					deltahat <- deltahat0 
				}
				dimnames(Gammahatx) <- list(vnames, paste("Dir", 1:d, sep=""))
				Deltahat <- deltahat*diag(p)
				dimnames(Deltahat) <- list(vnames, vnames)

				loglik <- - 0.5*n*p*(1+log(2*pi)) - 0.5*n*log(prod(deltahat))
				numpar <- p + d*(p-d) + ncol(fy)*d + p 
				aic <- -2*loglik + 2*numpar
				bic <- -2*loglik + log(n)*numpar

				ans <- list(Betahat=Betahatx, Gammahat=orthonorm(Gammahatx), Deltahat=Deltahat, evalues=fit$evalues, 
						loglik=loglik, aic=aic, bic=bic, numpar=numpar, numdir.test=numdir.test)
		
				return(ans)
			}

			Deltahat <- Betahat <- Gammahat <- vector("list")
			aic <- bic <- numpar <- loglik <- vector(length=eff.numdir+1)

			# No fitting values (eff.numdir=0)
			ev <- eigen(Sigmahat); 
			loglik[1] <- - 0.5*n*p*(1+log(2*pi)) - 0.5*n*log(prod(ev$values))
			numpar[1] <- p + p
			aic[1] <- -2*loglik[1] + 2*numpar[1]
			bic[1] <- -2*loglik[1] + log(n)*numpar[1]

			for (i in 1:eff.numdir)
			{
				Betahatx <- fit$Betahat[[i]]; Gammahatx <- fit$Gammahat[[i]] 
				Xc <- scale(X, TRUE, FALSE) - fy%*%t(Gammahatx%*%Betahatx)
				deltahat <- diag(t(Xc)%*%(Xc)/n)

				repeat
				{
					Xnew = X%*%((1/sqrt(deltahat))*diag(p))
					fit2 <- pfc(X=Xnew, y=y, fy=fy, numdir=i, structure="iso", numdir.test=FALSE)
					Betahatx <- fit2$Betahat; Gammahatx <- (diag(p)*sqrt(deltahat))%*%fit2$Gammahat 
					Xc <- scale(X, TRUE, FALSE) - fy%*%t(Gammahatx%*%Betahatx)
					deltahat0 <- diag(t(Xc)%*%(Xc)/n)
					if (sum(abs(deltahat-deltahat0)) < eps_aniso) break
					deltahat <- deltahat0 
				}

				Deltahat[[i]] <- deltahat*diag(p)
				dimnames(Deltahat[[i]]) <- list(vnames, vnames)
				
				loglik[i+1] = - 0.5*n*p*(1+log(2*pi)) - 0.5*n*log(prod(deltahat))
				numpar[i+1] <- p + (p-i)*i + i*dim(fy)[2] + p
				aic[i+1] <- -2*loglik[i+1] + 2*numpar[i+1]
				bic[i+1] <- -2*loglik[i+1] + log(n)*numpar[i+1]
				Betahat[[i]] <- Betahatx
				Gammahat[[i]] <- orthonorm(Gammahatx)
				dimnames(Gammahat[[i]]) <- list(vnames, paste("Dir", 1:i, sep=""))
			}
			ans <- list(Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, evalues=fit2$evalues, 
						loglik=loglik, aic=aic, bic=bic, numpar=numpar, numdir.test=numdir.test)

			return(ans)	
		}

		fit <- aniso(X=X, y=y, fy=fy, d=eff.numdir, eps_aniso=eps_aniso, numdir.test=numdir.test)

		ans <- list(Muhat=Muhat, Betahat=fit$Betahat, Gammahat=fit$Gammahat, Deltahat=fit$Deltahat, model="pfc",
				loglik=fit$loglik, aic=fit$aic, bic=fit$bic, numpar=fit$numpar, numdir=eff.numdir, 
				evalues=fit$evalues, structure="aniso", Xc=Xc, y=y, fy=fy,  call=match.call(), numdir.test=fit$numdir.test)
		
		if (numdir.test==FALSE) ans$R <- X%*%orthonorm(((fit$Deltahat)%^%(-1))%*%fit$Gammahat) else 
			ans$R <- X%*%orthonorm(((fit$Deltahat[[eff.numdir]])%^%(-1))%*%fit$Gammahat[[eff.numdir]])

		class(ans)<- "pfc"

		return(ans)
	}

	if (structure=="unstr")
	{
		unstr<-function(i)
		{
			sqrt_Sigmahat_res <- Sigmahat_res%^%0.5 
			Inv_Sqrt_Sigmahat_res <- solve(sqrt_Sigmahat_res)
			lf_matrix <- Inv_Sqrt_Sigmahat_res%*%Sigmahat_fit%*%Inv_Sqrt_Sigmahat_res
			all_evalues <- eigen(lf_matrix, symmetric=T)$values
			evalues <- all_evalues[1:i]

			Vhat <- eigen(lf_matrix, symmetric=T)$vectors
			Vhati <- matrix(Vhat[,1:i], ncol=i)
			Gammahat <- (Sigmahat_res%^%0.5)%*%Vhati%*%solve((t(Vhati)%*%Sigmahat_res%*%Vhati)%^%0.5)  
			dimnames(Gammahat)<- list(vnames, paste("Dir", 1:i, sep=""))

			Khat<-diag(0, p) 
			if (i < min(ncol(fy),p)) {diag(Khat)[(i+1):min(ncol(fy), p )]<- all_evalues[(i+1):min(ncol(fy), p)]}
			Deltahat <- sqrt_Sigmahat_res%*%Vhat%*%(diag(p)+Khat)%*%t(Vhat)%*%sqrt_Sigmahat_res
			dimnames(Deltahat) <- list(vnames, vnames)
			Betahat <- ((t(Vhati)%*%Sigmahat_res%*%Vhati)%^%0.5)%*%t(Vhati)%*%solve(Sigmahat_res%^%0.5)%*%t(Xc)%*%fy%*% solve(t(fy)%*%fy)

			temp0 <- -(n*p/2)*(1 + log(2*pi))
			temp1 <- -(n/2)*log(det(Sigmahat_res)) 
			temp2 <- 0; 

			if (i < min(ncol(fy),p)) temp2 <- -(n/2)*sum(log(1 + all_evalues[(i+1):p]))
			loglik <- temp0 + temp1 + temp2
			numpar <- p + (p-i)*i + i*ncol(fy) + p*(p+1)/2
			aic <- -2*loglik + 2*numpar
			bic <- -2*loglik + log(n)*numpar

			return(list(Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, evalues=evalues, 
				loglik=loglik, aic=aic, bic=bic, numpar=numpar))
		}

		if (identical(numdir.test, FALSE))
		{
			out <- unstr(eff.numdir)

			ans <- list(R=X%*%orthonorm(solve(out$Deltahat)%*%out$Gammahat), Muhat=Muhat, Betahat=out$Betahat, Gammahat=out$Gammahat, 
					Deltahat=out$Deltahat, evalues=out$evalues, loglik=out$loglik, aic=out$aic, bic=out$bic, numpar=out$numpar, 
					numdir=eff.numdir, model="pfc", structure="unstr", y=y, fy=fy, Xc=Xc,  call=match.call(), numdir.test=numdir.test)

			class(ans) <- "pfc"	
			return(ans);	
		}

		aic <- bic <- numpar <- loglik <- vector(length=eff.numdir+1)
		evalues <- vector(length=eff.numdir)
		Betahat <- Deltahat <- Gammahat <-vector("list") 
 
		loglik[1] <- - 0.5*n*p*(1+log(2*pi)) - 0.5*n*log(det(Sigmahat)) 
		numpar[1] <- p + p*(p+1)/2
    aic[1] <- -2*loglik[1] + 2*numpar[1]
		bic[1] <- -2*loglik[1] + log(n)*numpar[1]
		Deltahat[[1]] <- Sigmahat
    dimnames(Deltahat[[1]]) <- list(vnames, vnames) 

		for (i in 1:eff.numdir)
		{
			fit <- unstr(i);
			Betahat[[i]] <-fit$Betahat 
			Gammahat[[i]] <-fit$Gammahat 
			Deltahat[[i]] <- fit$Deltahat
			loglik[i+1] <- fit$loglik 
			numpar[i+1] <- fit$numpar
			aic[i+1] <- fit$aic
			bic[i+1] <- fit$bic	
		}
		ans <- list(R=X%*%orthonorm(solve(Deltahat[[eff.numdir]])%*%Gammahat[[eff.numdir]]), Muhat=Muhat, 
				Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, evalues=fit$evalues, loglik=loglik, 
				aic=aic, bic=bic, numpar=numpar, numdir=eff.numdir, model="pfc", structure="unstr", y=y, 
				fy=fy, Xc=Xc, call=match.call(), numdir.test=numdir.test)

		class(ans)<- "pfc"
		return(ans)
	}
	else if (structure=="unstr2")
	{
		unstr2 <- function(i)
		{
			objfun <- function(W)
			{
				Qt <- W$Qt; dc <- W$dim[1] 
				p <- ncol(Qt); S <- W$Sigmas
				U <- matrix(Qt[,1:dc], ncol=dc)	
				V <- matrix(Qt[,(dc+1):p], ncol=(p-dc))
				value <- -(n/2)*(p*log(2*pi)+p+log(det(t(V)%*%S$Sigmahat%*%V))+ log(det(t(U)%*%S$Sigmahat_res%*%U)))

				terme1 <- solve(t(U)%*%S$Sigmahat_res%*%U)%*%(t(U)%*%S$Sigmahat_res%*%V);
				terme2 <- (t(U)%*%S$Sigmahat%*%V)%*%solve(t(V)%*%S$Sigmahat%*%V)

				gradient <- 2*(terme1 - terme2)
				return(list(value=value, gradient=gradient))
			}
			sigmas <- list(Sigmahat=Sigmahat, Sigmahat_fit=Sigmahat_fit, Sigmahat_res=Sigmahat_res, p=p, n=n)

			W <- list(Qt = svd(Sigmahat_fit)$u, dim=c(numdir, p), Sigmas=list(Sigmahat=Sigmahat, Sigmahat_fit=Sigmahat_fit, 
						Sigmahat_res=Sigmahat_res, p=p, n=n))

			objfun <- assign("objfun", objfun, envir=.BaseNamespaceEnv) 
			grassoptim <- GrassmannOptim(objfun, W,...)

			Gammahat <- matrix(grassoptim$Qt[,1:i], ncol=i, dimnames=list(vnames, paste("Dir", 1:i, sep=""))) 

			Gammahat0 <- matrix(grassoptim$Qt[, (i+1):p], ncol=p-i, dimnames=list(vnames, paste("Dir", (i+1):p, sep="")))

			Betahat <- t(Gammahat)%*%t(Xc)%*%fy%*%solve(t(fy)%*%fy)
			Omegahat <- t(Gammahat)%*%Sigmahat_res%*%Gammahat
			Omegahat0 <-t(Gammahat0)%*%Sigmahat%*%Gammahat0
			Deltahat <- Gammahat%*%Omegahat%*%t(Gammahat) + Gammahat0%*%Omegahat0%*%t(Gammahat0)
			dimnames(Deltahat) <- list(vnames, vnames)

			temp0 <- -(n*p/2)*(1+log(2*pi))
			temp1 <- -(n/2)*log(det(t(Gammahat)%*%Sigmahat_res%*%Gammahat)) 
			temp2 <- -(n/2)*log(det(t(Gammahat0)%*%Sigmahat%*%Gammahat0)) 

			loglik <- temp0 + temp1 + temp2 
			numpar <- p + (p-i)*i + i*dim(fy)[2] + i*(i+1)/2 + (p-i)*(p-i+1)/2 
			aic <- -2*loglik + 2*numpar
			bic <- -2*loglik + log(n)*numpar
			ev.fit <- eigen(Sigmahat_fit) 
			evalues <- ev.fit$values[1:i]

			return(list(Betahat=Betahat, Gammahat=Gammahat, Gammahat0=Gammahat0, Omegahat=Omegahat, Omegahat0=Omegahat0, 
					Deltahat=Deltahat, evalues=evalues, loglik=loglik, aic=aic, bic=bic, numpar=numpar))
		}

		if (identical(numdir.test, FALSE))
		{
			out <- unstr2(numdir)

			ans <- list(R = X%*%out$Gammahat, Muhat=Muhat, Betahat=out$Betahat, Gammahat=out$Gammahat, 
					Gammahat0=out$Gammahat0, Omegahat=out$Omegahat, Omegahat0=out$Omegahat0, 
					Deltahat=out$Deltahat, evalues=out$evalues, loglik=out$loglik, aic=out$aic,
					bic=out$bic, numpar=out$numpar, numdir=numdir, model="pfc", structure="unstr2", 
					y=y, fy=fy, Xc=Xc,  call=match.call(), numdir.test=numdir.test)

			class(ans) <- "pfc"	
			return(ans)
		}

		aic <- bic <- numpar <- loglik <- vector(length=numdir+1)
		Betahat <- Deltahat <- Gammahat <- Gammahat0 <- Omegahat <- Omegahat0 <- vector("list") 

		loglik[1] <- -(n*p/2)*(log(2*pi) + (1+log(Trace(Sigmahat)/p)))
		numpar[1] <- p + p*(p+1)/2
		aic[1] <- -2*loglik[1] + 2*numpar[1]
		bic[1] <- -2*loglik[1] + log(n)*numpar[1]

		for(m in 1:numdir)
		{
			fit <- unstr2(m);
			Betahat[[m]] <-fit$Betahat 
			Gammahat[[m]] <-fit$Gammahat 
			Omegahat[[m]] <- fit$Omegahat 
			Omegahat0[[m]] <- fit$Omegahat0
			Deltahat[[m]] <- fit$Deltahat 
			loglik[m+1] <- fit$loglik 
			numpar[m+1] <- fit$numpar
			aic[m+1] <- fit$aic 
			bic[m+1] <- fit$bic	
		}
		ans <- list(R = X%*%Gammahat[[numdir]], evalues=fit$evalues, loglik =loglik, aic=aic, bic=bic, numdir=numdir, 
				numpar=numpar, Muhat=Muhat, Betahat=Betahat, Gammahat=Gammahat,	
				Gammahat0=Gammahat0, Omegahat=Omegahat, Omegahat0=Omegahat0, 
				Deltahat=Deltahat, model="pfc", structure="unstr2", 
				y=y, fy=fy, Xc=Xc,  call=match.call(), numdir.test=numdir.test)

		class(ans)<- "pfc";
		return(ans)
	}
}

summary.pfc <-
  function(object,...)
  {
    "%^%"<-function(M, pow) 
    { 
      if (prod(dim(M)==list(1,1))) return( as.matrix(M^pow) )
      eigenM = eigen(M) 
      return(eigenM$vectors%*%diag(c(eigenM$values)^pow)%*%t(eigenM$vectors))  
    }
    
    lrtest <- function(object)
    {
      Tests <- Dfs <- Pvals <-cnames <- vector(length=object$numdir)
      for (w in 1:object$numdir)
      {
        Tests[w] <- 2*(object$loglik[object$numdir+1]-object$loglik[w])
        Dfs[w] <- object$numpar[object$numdir+1]-object$numpar[w]
        Pvals[w] <- 1-pchisq(Tests[w], df=Dfs[w])
        cnames[w] <- paste(paste(w-1, "D vs >= ", sep=""), paste(w, "D", sep=""), sep="")
      }
      ans <- data.frame(cbind(Tests, Dfs, Pvals)) 
      rownames(ans) <- cnames
      colnames(ans) <- c("Stat", "df", "p.value")
      return(ans)
    }
    
    ics <- function(object)
    {
      mat <- data.frame(rbind(object$aic, object$bic))
      v1 <- c("aic", "bic") 
      v2 <-paste("d=", 0:(object$numdir), sep="")
      dimnames(mat) = list(v1, v2)
      return(mat)
    }
    
    r2 <- function(object)
    {
      R2 <- cnames <- vector(length=object$numdir)
      
      for (w in 1:object$numdir)
      {
        if (identical(object$numdir.test, FALSE)) tempd <-data.frame(cbind(object$y-mean(object$y), 
                                                                           object$Xc%*%((object$Deltahat)%^%(-1))%*%matrix(object$Gammahat[, 1:w], ncol=w))) 
        if (identical(object$numdir.test, TRUE)) tempd <-data.frame(cbind(object$y-mean(object$y), 
                                                                          object$Xc%*%((object$Deltahat[[w]])%^%(-1))%*%object$Gammahat[[w]])) 
        colnames(tempd)[1] <- "y"
        R2[w] <- summary(lm(y~.-1, data=tempd))$r.squared
      }
      Rsq <- data.frame(rbind(object$evalues, R2))
      v1 <- c("Eigenvalues", "R^2(OLS|pfc)") 
      v2 <-paste("Dir", 1:object$numdir, sep="")
      dimnames(Rsq) = list(v1, v2)
      return(Rsq)
    }
    ans <- object; class(ans) <-"summarypfc"
    
    if (identical(object$numdir.test, TRUE))
    {	
      ans$LRT <- lrtest(object)
      ans$IC <- ics(object)
      if (is.numeric(object$y)) ans$Rsq <- r2(object)
      return(ans)
    }
    return(ans)
  }


plot.pfc <-
  function(x,...)
  {
    if (is.numeric(x$y))
    {
      dd <- NCOL(x$R);
      
      if (dd==1)
      { 
        plot(x$y, x$R, pch=20, ylab="First Component", xlab="y", ...); 
        lines(lowess(x$y, x$R), lwd=2, col="red");
        
        title(main="Response against the reduction")
      }
      
      if (dd==2) 
      {
        par(mfrow=c(1,2), oma=c(0,0,2,0)); 
        
        plot(x$y, x$R[,1], pch=20, ylab="First Component", xlab="y",...); 
        lines(lowess(x$y, x$R[,1]), lwd=2, col="red");
        
        plot(x$y, x$R[,2], pch=20, ylab="Second Component", xlab="y",...);
        lines(lowess(x$y, x$R[,2]), lwd=2, col="red");
        
        title(main="Components of the reduction vs Response", outer=TRUE)		
      }
      
      if (dd==3) 
      {
        par(mfrow=c(2,2), oma=c(0,0,2,0)); 
        
        plot(x$y, x$R[,1], pch=20, ylab="First Component", xlab="y",...); 
        lines(lowess(x$y, x$R[,1]), lwd=2, col="red");
        
        plot(x$y, x$R[,2], pch=20, ylab="Second Component", xlab="y",...); 
        lines(lowess(x$y, x$R[,2]), lwd=2, col="red");
        
        plot(x$y, x$R[,3], pch=20, ylab="Third Component", xlab="y",...);
        lines(lowess(x$y, x$R[,3]), lwd=2, col="red");
        
        title(main="Components of the reduction vs Response", outer=TRUE)
      }
      
      if (dd>3) 
      {
        par(mfrow=c(2,2), oma=c(0,0,2,0)); 
        
        plot(x$y, x$R[,1], pch=20, ylab="First Component", xlab="y",...); 
        lines(lowess(x$y, x$R[,1]), lwd=2, col="red");
        
        plot(x$y, x$R[,2], pch=20, ylab="Second Component", xlab="y",...); 
        lines(lowess(x$y, x$R[,2]), lwd=2, col="red");
        
        plot(x$y, x$R[,3], pch=20, ylab="Third Component", xlab="y",...); 
        lines(lowess(x$y, x$R[,3]), lwd=2, col="red");
        
        plot(x$y, x$R[,4], pch=20, ylab="Fourth Component", xlab="y",...);
        lines(lowess(x$y, x$R[,4]), lwd=2, col="red");
        
        title(main="Components of the reduction vs Response", outer=TRUE)
      }
    }
    
    if (is.factor(x$y))
    {
      mycolors <- c("black", "blue", "red", "green", "yellow", "gray", "cyan", "magenta")
      pchs <- c(20, 21, 22, 17, 1, 2, 3, 4)
      nl <- length(unique(x$y))
      mycol <- mypch <- as.integer(factor(x$y, levels=unique(x$y)))
      for (i in 1:nl)  { mycol[mycol==i]<- mycolors[i]; mypch[mypch==i] <- pchs[i]}
      
      if (NCOL(x$R)==1)
      {
        par(mfrow=c(1,2), oma=c(0,0,2,0)) 
        plot(x$R[,1], xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1.5, xaxt="n", yaxt="n", axes=FALSE)
        legend(0.2, 0.8, unique(x$y), border = "blank", cex = 1, title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])
        
        plot(x$R[,1], xlab="index", ylab="PFC - Dir1", pch=mypch, col=mycol, cex=1.5)
        title(main="Scatterplot of the sufficient reduction", outer=TRUE)
      }
      
      
      if (NCOL(x$R)==2)
      {
        par(mfrow=c(1,2), oma=c(0,0,2,0)) 
        plot(x$R[,1], xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1.5, xaxt="n", yaxt="n", axes=FALSE)
        legend(0.2, 0.8, unique(x$y), border = "blank", cex = 1, title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])
        
        plot(x$R[,1], x$R[,2], xlab="PFC - Dir1", ylab="PFC - Dir2", pch=mypch, col=mycol, cex=1)
        title(main="Scatterplot of the components of the sufficient reduction", outer=TRUE)
      }
      
      if (NCOL(x$R) == 3)
      {
        par(mfrow=c(2,2), oma=c(0,0,2,0)); 
        plot(x$R[,1], xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1.5, xaxt="n", yaxt="n", axes=FALSE)
        legend(0.2, 0.8, unique(x$y), border = "blank", cex = 1, title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])
        
        plot(x$R[,1], x$R[,2], xlab="PFC - Dir1", ylab="PFC - Dir2", pch=mypch, col=mycol, cex=1)
        plot(x$R[,1], x$R[,3], xlab="PFC - Dir1", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
        plot(x$R[,2], x$R[,3], xlab="PFC - Dir2", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
        title(main="Scatterplots of the components of the sufficient reduction", outer=TRUE)
      }
      
      if (NCOL(x$R) > 3)
      {
        par(mfrow=c(3,3), oma=c(0,0,2,0));  
        plot(x$R[,1], xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1.5, xaxt="n", yaxt="n", axes=FALSE)
        legend(0.2, 0.8, unique(x$y), border = "blank", cex = 1, title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])
        
        plot(x$R[,1], x$R[,2], xlab="PFC - Dir1", ylab="PFC - Dir2", pch=mypch, col=mycol, cex=1)
        plot(x$R[,1], x$R[,3], xlab="PFC - Dir1", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
        plot(x$R[,1], x$R[,4], xlab="PFC - Dir2", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
        plot(x$R[,2], x$R[,3], xlab="PFC - Dir2", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
        plot(x$R[,2], x$R[,4], xlab="PFC - Dir2", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
        title(main="Scatterplots of the components of the sufficient reduction", outer=TRUE)
      }
      
    }
  }


library('plot3D')



par(mfrow =c(1,1))
setwd('E:\\ALLDocs\\2022(1)Projs\\NomparaProj')
data=read.table("BostonHousePrice.txt",head=TRUE)


train.x=data[1:355,c(1,6,5,11,10,7,13)]
train.y=data[1:355,14]
test.x=data[356:506,c(1,6,5,11,10,7,13)]
test.y=data[356:506,14]

Y=data[,14]
X=data[,c(1,6,5,11,10,7,13)]

fit1 <- pfc(X=X, y=Y, fy=bf(y=Y, case="poly",degree=3),numdir=3, structure="iso")


pfc_fitted <- function(model,rs=NULL,x=NULL,mu=NULL,std=NULL,weighted=FALSE){
  Reduc=model$R
  y=model$y
  if(is.null(rs)){
    rs=as.matrix(x)%*%orthonorm(solve(model$Deltahat)%*%model$Gammahat)
  }
  ans=vector('numeric',dim(rs)[1])
  for (j in 1:dim(rs)[1]) {
    w=vector('numeric',dim(Reduc)[1])
    for (i in 1:dim(Reduc)[1]) {
      w[i]<-exp(-0.5*(Reduc[i,]-rs[j,])%*%(Reduc[i,]-rs[j,]))
      if(weighted==TRUE){
        w[i]<-(dnorm(y[i],mu,std))^0*exp(-0.5*(Reduc[i,]-rs[j,])%*%(Reduc[i,]-rs[j,]))
      }
    }
    w <- w/sum(w)
    ans[j]<-w%*%y
  }
  return(ans)
}




X.copy<-X
Y<-BoxCox(Y,0)
X[,1]<-BoxCox(X[,1],0)
X[,2]<-BoxCox(X[,2],1)
X[,3]<-BoxCox(X[,3],0)
X[,4]<-BoxCox(X[,4],3)
X[,5]<-BoxCox(X[,5],-1)
X[,6]<-BoxCox(X[,6],0)
X[,7]<-BoxCox(X[,7],0)

hist(Y)
par(new=TRUE)
plot(seq(1.5,4,0.1),dnorm(seq(1.5,4,0.1),mean(Y),sd(Y)),type = 'l')

X <- scale(X)
X.copy<-scale(X.copy)


t1=proc.time()
res=0
for (i in 1:length(Y)) {
  fit<-pfc(X=X[-i,], y=exp(Y)[-i], fy=bf(y=exp(Y)[-i], case="poly",degree = 4),numdir=2, structure="iso")
  r=t(as.matrix(X[i,]))%*%orthonorm(solve(fit$Deltahat)%*%fit$Gammahat)
  yhat=pfc_fitted(fit,r)
  res<-res+(yhat-exp(Y[i]))^2
}
res/506
t2=proc.time()
t2-t1



fit<-pfc(X=train.x,y=train.y, fy=bf(y=train.y, case="poly",degree=4),numdir=2, structure="aniso")
(pfc_fitted(fit,x=test.x)-test.y)%*%(pfc_fitted(fit,x=test.x)-test.y)/151

plot(pfc_fitted(fit1,SIRhat),Y)
lines(c(0,4),c(0,4))

fit<-pfc(X=X,y=exp(Y), fy=bf(y=exp(Y), case="poly",degree=4),numdir=2, structure="aniso")
PFChat<-fitpfc$R
X1=seq(min(PFChat[,1]),max(PFChat[,1]),0.4)
X2=seq(min(PFChat[,2]),max(PFChat[,2]),0.4)
Yhat=matrix(0,length(X1),length(X2))
pfc_fitted(fit,t(as.matrix(c(X1[1],X2[2]))))
for (i in 1:length(X1)) {
  for (j in 1:length(X2)) {
    Yhat[i,j]=pfc_fitted(fit,t(as.matrix(c(X1[i],X2[j]))))
  }
}
plot3d(PFChat[,1],PFChat[,2],exp(Y),col = 'red')
surface3d(X1,X2,Yhat,type='s',alpha=0.4,lit=FALSE,front='lines',back='lines')
rgl.postscript('PFC_Raw.pdf', fmt = 'pdf')



for (i in 1:length(X1)) {
  for (j in 1:length(X2)) {
    toplot[1,(i-1)*length(X2)+j]=X1[i]
    toplot[2,(i-1)*length(X2)+j]=X2[j]
    toplot[3,(i-1)*length(X2)+j]=exp(Yhat[i,j])
  }
}



write.table(cbind(Rhat,exp(Y)),file ='TrueResult' ,sep =',',row.names=FALSE,col.names = FALSE)
write.table(as.matrix(t(toplot)),file ='PFCResult' ,sep =',',row.names=FALSE,col.names = FALSE)






res
resnew



library(regpro)
  t1=proc.time()
  resnw=0
  for (i in 1:length(Y)) {
    yhat=kernesti.regr(as.matrix(scale(X)[i,]),as.matrix(scale(X)[-i,]),Y[-i], h=1.06*506^(-1/6), kernel="gauss", g=NULL, gernel="gauss", vect=FALSE)
    resnw<-resnw+(exp(yhat)-exp(Y[i]))^2
  }
  resnw/506
  t2=proc.time()
  t2-t1
  


    t1=proc.time()
  resnw=0
  for (i in 1:length(Y)) {
    s0<-dr(MEDV~CRIM+RM+NOX+PTRATIO+TAX +AGE+LSTAT,data = rawfulldata[-i,], slice.function=dr.slices.arc,nslices=8, numdir = 2, method = "sir")
    r0=scale(as.matrix(X.copy)%*%s0$evectors[,1:2])
    fit<-loess(exp(Y)~Dir1+Dir2,data = as.data.frame(cbind(r0,Y)))
    yhat=predict(fit,  t(as.data.frame(r0[i,])))
    reslp<-reslp+(yhat-exp(Y[i]))^2
  }
  resnw/506
  t2=proc.time()
  t2-t1
  
  PFChat
  
  t1=proc.time()
  reslp=0
  for (i in 1:length(Y)) {
    fit<-pfc(X=X[-i,], y=exp(Y)[-i], fy=bf(y=exp(Y)[-i], case="poly",degree=4),numdir=2, structure="aniso")
    r0=scale(as.matrix(X)%*%orthonorm(solve(fit$Deltahat)%*%fit$Gammahat))
    fit<-loess(exp(Y)~Dir1+Dir2,data = as.data.frame(cbind(r0,Y)))
    yhat=predict(fit,  t(as.data.frame(r0[i,])))
    reslp<-reslp+(yhat-exp(Y[i]))^2
  }
  reslp/506
  t2=proc.time()
  t2-t1

t1=proc.time()
t1=proc.time()
  resnw=0
  for (i in 1:length(Y)) {
    fit<-pfc(X=X[-i,], y=exp(Y)[-i], fy=bf(y=exp(Y)[-i], case="poly",degree=4),numdir=2, structure="iso")
    r0=scale(as.matrix(X)%*%orthonorm(solve(fit$Deltahat)%*%fit$Gammahat))
    yhat=kernesti.regr(r0[i,],r0[-i,],Y[-i], h=1.06*506^(-1/6), kernel="gauss", g=NULL, gernel="gauss", vect=FALSE)
    resnw<-resnw+(exp(yhat)-exp(Y[i]))^2
  }
  resnw/506
  t2=proc.time()
  t2-t1
  
fulldata<-as.data.frame(cbind(X,exp(Y)))
colnames(fulldata)[8]<-'MEDV'
library(dr)
colnames(rawfulldata)[8]<-'MEDV'
rawfulldata


fitpfc<-pfc(X=X,y=exp(Y), fy=bf(y=exp(Y), case="poly",degree=4),numdir=2, structure="aniso")
PFChat<-fitpfc$R
fit<-sm.regression(PFChat,exp(Y),h=c(1.06*sd(PFChat[,1])*506^(-1/6),1.06*sd(PFChat[,2])*506^(-1/6)))
plot3d(PFChat[,1],PFChat[,2],exp(Y),col = 'red')
surface3d(fit$eval.points[,1],fit$eval.points[,2],fit$estimate,type='s',alpha=0.4,lit=FALSE,front='lines',back='lines')
rgl.postscript('PFC_NW.pdf', fmt = 'pdf')


s0<-dr(MEDV~CRIM+RM+NOX+PTRATIO+TAX +AGE+LSTAT,data = rawfulldata, slice.function=dr.slices.arc,nslices=5, numdir = 2, method = "sir")
SIRhat<-as.matrix(X.copy)%*%s0$evectors[,1:2]
fit<-sm.regression(SIRhat,exp(Y),h=c(1.06*sd(SIRhat[,1])*506^(-1/6),1.06*sd(SIRhat[,2])*506^(-1/6)))
library(rgl)
plot3d(SIRhat[,1],SIRhat[,2],exp(Y),col = 'red')
surface3d(fit$eval.points[,1],fit$eval.points[,2],fit$estimate,type='s',alpha=0.4,lit=FALSE,front='lines',back='lines')
rgl.postscript('SIR_NW.pdf', fmt = 'pdf')

fitpfc<-pfc(X=X,y=exp(Y), fy=bf(y=exp(Y), case="poly",degree=3),numdir=2, structure="aniso")
PFChat<-fitpfc$R
fit<-loess(exp(Y)~Dir1+Dir2,data = as.data.frame(cbind(PFChat,Y)))
X1=seq(0.8*min(PFChat[,1]),0.8*max(PFChat[,1]),0.4)
X2=seq(0.8*min(PFChat[,2]),0.8*max(PFChat[,2]),0.4)
Yhat=matrix(0,length(X1),length(X2))
for (i in 1:length(X1)) {
  for (j in 1:length(X2)) {
    temp=t(as.data.frame(c(X1[i],X2[j])))
    colnames(temp)<-c('Dir1','Dir2')
    Yhat[i,j]=predict(fit,temp)
  }
}
plot3d(PFChat[,1],PFChat[,2],exp(Y),col = 'red')
surface3d(X1,X2,Yhat,type='s',alpha=0.4,lit=FALSE,front='lines',back='lines')
rgl.postscript('PFC_LPR.pdf', fmt = 'pdf')


s0<-dr(MEDV~CRIM+RM+NOX+PTRATIO+TAX +AGE+LSTAT,data = rawfulldata, slice.function=dr.slices.arc,nslices=5, numdir = 2, method = "sir")
SIRhat<-as.matrix(X.copy)%*%s0$evectors[,1:2]
fit<-loess(exp(Y)~Dir1+Dir2,data = as.data.frame(cbind(SIRhat,Y)))
X1=seq(0.8*min(SIRhat[,1]),0.8*max(SIRhat[,1]),0.4)
X2=seq(0.8*min(SIRhat[,2]),0.8*max(SIRhat[,2]),0.4)
Yhat=matrix(0,length(X1),length(X2))
for (i in 1:length(X1)) {
  for (j in 1:length(X2)) {
    temp=t(as.data.frame(c(X1[i],X2[j])))
    colnames(temp)<-c('Dir1','Dir2')
    Yhat[i,j]=predict(fit,temp)
  }
}
library(rgl)
plot3d(SIRhat[,1],SIRhat[,2],exp(Y),col = 'red')
surface3d(X1,X2,Yhat,type='s',alpha=0.4,lit=FALSE,front='lines',back='lines')
rgl.postscript('SIR_LPR.pdf', fmt = 'pdf')

write.table(cbind(scale(X),scale(SIRhat),scale(fitpfc$R),exp(Y)),file ='Tolocal.txt' ,sep =',',row.names=FALSE,col.names = FALSE)


summary(fit)
library('glmnet')

t1=proc.time()
reslm=0
for (i in 1:length(Y)) {
  fit<-lm(log(MEDV)~.-INDUS-AGE,data=data[-i,])
  yhat=as.numeric(predict.lm(fit,newdata = as.data.frame(data[i,])))
  reslm<-reslm+(exp(yhat)-exp(Y[i]))^2
}
reslm/506
t2=proc.time()
t2-t1


install.packages('gplm')
library(gplm)

t1=proc.time()
for (i in 1:506) {
  fit<-kbackfit(SIRhat,exp(Y),h=1)
}
t2=proc.time()
t2-t1
fit$rss/506

install.packages('gam')
library('gam')

temp
fit<-gam(exp(Y)~Dir1+Dir2,data = as.data.frame(cbind(r0,Y)))
predict(fit,as.data.frame(SIRhat))

t1=proc.time()
reslp=0
for (i in 1:length(Y)) {
  s0<-dr(MEDV~CRIM+RM+NOX+PTRATIO+TAX +AGE+LSTAT,data = rawfulldata[-i,], slice.function=dr.slices.arc,nslices=8, numdir = 2, method = "sir")
  r0=scale(as.matrix(X.copy)%*%s0$evectors[,1:2])
  fit<-gam(exp(Y)~lo(Dir1)+lo(Dir2),data = as.data.frame(cbind(r0,Y)))
  yhat=predict(fit, as.data.frame(t(r0[i,])))
  reslp<-reslp+(yhat-exp(Y[i]))^2
}
reslp/506
t2=proc.time()
t2-t1

s0<-dr(MEDV~CRIM+RM+NOX+PTRATIO+TAX +AGE+LSTAT,data = rawfulldata, slice.function=dr.slices.arc,nslices=5, numdir = 2, method = "sir")
SIRhat<-as.matrix(X.copy)%*%s0$evectors[,1:2]
fit<-gam(exp(Y)~lo(Dir1)+lo(Dir2),data = as.data.frame(cbind(SIRhat,Y)))
X1=seq(0.8*min(SIRhat[,1]),0.8*max(SIRhat[,1]),0.4)
X2=seq(0.8*min(SIRhat[,2]),0.8*max(SIRhat[,2]),0.4)
Yhat=matrix(0,length(X1),length(X2))
for (i in 1:length(X1)) {
  for (j in 1:length(X2)) {
    temp=as.data.frame(t(c(X1[i],X2[j])))
    colnames(temp)<-c('Dir1','Dir2')
    Yhat[i,j]=predict(fit,temp)
  }
}
library(rgl)
plot3d(SIRhat[,1],SIRhat[,2],exp(Y),col = 'red')
surface3d(X1,X2,Yhat,type='s',alpha=0.4,lit=FALSE,front='lines',back='lines')
rgl.postscript('SIR_LPR_AM.pdf', fmt = 'pdf')

fitpfc<-pfc(X=X,y=exp(Y), fy=bf(y=exp(Y), case="poly",degree=4),numdir=2, structure="aniso")
PFChat<-fitpfc$R
fit<-gam(exp(Y)~lo(Dir1)+lo(Dir2),data = as.data.frame(cbind(PFChat,Y)))
X1=seq(min(PFChat[,1]),max(PFChat[,1]),0.4)
X2=seq(min(PFChat[,2]),max(PFChat[,2]),0.4)
Yhat=matrix(0,length(X1),length(X2))
for (i in 1:length(X1)) {
  for (j in 1:length(X2)) {
    temp=as.data.frame(t(c(X1[i],X2[j])))
    colnames(temp)<-c('Dir1','Dir2')
    Yhat[i,j]=predict(fit,temp)
  }
}
plot3d(PFChat[,1],PFChat[,2],exp(Y),col = 'red')
surface3d(X1,X2,Yhat,type='s',alpha=0.4,lit=FALSE,front='lines',back='lines')
rgl.postscript('PFC_LPR_AM.pdf', fmt = 'pdf')

fulldata
t1=proc.time()
resnw=0
for (i in 1:length(Y)) {
  fit<-gam(log(MEDV)~ lo(CRIM)+lo(RM)+lo(NOX)+lo(PTRATIO)+lo(TAX)+lo(AGE)+lo(LSTAT),data=fulldata[-i,])
  yhat=as.numeric(predict(fit,newdata = as.data.frame(fulldata[i,])))
  resnw<-resnw+(exp(yhat)-exp(Y[i]))^2
}
resnw/506
t2=proc.time()
t2-t1