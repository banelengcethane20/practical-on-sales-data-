#(i)
#Finding ANOVA SUMS
library(datarium)
t(tail(marketing, n=8))
X_matrix=cbind(c(1,1,1,1,1,1,1,1),
                 c(20.64,200.16,179.64,45.84,113.04,212.40,340.32,278.52),
                   c(4.92,50.40,42.72,4.44,5.88,11.16,50.40,10.32),
                     c(37.92,4.32,7.20,16.56,9.72,7.68,79.44,10.44))
Y_vector=cbind(c(7.08,23.52,20.76,9.12,11.64,15.36,30.60,16.08))
X_matrix
Y_vector
n = 8
p = 4
J_matrix = matrix(1,8,8)
J_matrix

#calculation Anova Sums
y_transpose_y = t(Y_vector)%*%Y_vector
y_transpose_y
yT_J_matrix_y = (1/n)%*%t(Y_vector)%*%J_matrix%*%Y_vector
yT_J_matrix_y
SST = y_transpose_y-yT_J_matrix_y
SST
y_trans_x = t(Y_vector)%*%X_matrix
y_trans_x
x_trans_y = t(y_trans_x)
x_trans_y
x_trans_x = t(X_matrix)%*%X_matrix
x_trans_x
x_trans_x_Cofactors=function(M){
  stopifnot(length(unique(dim(M)))==1)
  co_factor=M
  for(i in 1:dim(M)[1]){
    for(j in 1:dim(M)[2]){
      co_factor[i,j]=(det(M[-i,-j])*(-1)^(i+j))
    }
  }
  return(co_factor)
}

# Co-Factors
x_trans_x_Cofactors(x_trans_x)

#Adjoint Matrix 
t(x_trans_x_Cofactors(x_trans_x))

det_x_trans_x = det(x_trans_x)
det_x_trans_x

#inverse of X_TRANS_X
x_trans_x_inverse=(1/det_x_trans_x)*t(x_trans_x_Cofactors(x_trans_x))
x_trans_x_inverse

yTx_x_trans_x_inverse_xTy= y_trans_x%*%x_trans_x_inverse%*%x_trans_y
yTx_x_trans_x_inverse_xTy

SSE = y_transpose_y - yTx_x_trans_x_inverse_xTy
SSE

SSR= yTx_x_trans_x_inverse_xTy - yT_J_matrix_y
SSR

# (ii)
#calculation of MS, F, Rsquared and Rsquared adj
X_matrix=cbind(c(1,1,1,1,1,1,1,1),
                 c(20.64,200.16,179.64,45.84,113.04,212.40,340.32,278.52),
                   c(4.92,50.40,42.72,4.44,5.88,11.16,50.40,10.32),
                     c(37.92,4.32,7.20,16.56,9.72,7.68,79.44,10.44))
Y_vector=cbind(c(7.08,23.52,20.76,9.12,11.64,15.36,30.60,16.08))
n = 8
p = 4
J_matrix = matrix(1,8,8)
y_transpose_y = t(Y_vector)%*%Y_vector
yT_J_matrix_y = (1/n)%*%t(Y_vector)%*%J_matrix%*%Y_vector
SST = y_transpose_y-yT_J_matrix_y
y_trans_x = t(Y_vector)%*%X_matrix
x_trans_y = t(y_trans_x)
x_trans_x = t(X_matrix)%*%X_matrix
x_trans_x_Cofactors=function(M){
  stopifnot(length(unique(dim(M)))==1)
  co_factor=M
  for(i in 1:dim(M)[1]){
    for(j in 1:dim(M)[2]){
      co_factor[i,j]=(det(M[-i,-j])*(-1)^(i+j))
    }
  }
  return(co_factor)
}
# Co-Factors
x_trans_x_Cofactors(x_trans_x)
#Adjoint Matrix 
t(x_trans_x_Cofactors(x_trans_x))
det_x_trans_x = det(x_trans_x)
x_trans_x_inverse=(1/det_x_trans_x)*t(x_trans_x_Cofactors(x_trans_x))
yTx_x_trans_x_inverse_xTy= y_trans_x%*%x_trans_x_inverse%*%x_trans_y
SSE = y_transpose_y - yTx_x_trans_x_inverse_xTy
SSR= yTx_x_trans_x_inverse_xTy - yT_J_matrix_y

MSR = SSR/(p-1)
MSR

MSE = SSE/(n-p)
MSE

MST = SST/(n-1)
MST

F_cal = MSR/MSE
F_cal

R_SQUARED = SSR/SST
R_SQUARED

R_SQUARED_ADJ = 1 - (MSE/MST)
R_SQUARED_ADJ

#ANOVA TABLE
Anova_Table = function (object, reg_collapse=TRUE,...) 
{
  if (length(list(object, ...)) > 1L) 
    return(anova.lmlist(object, ...))
  if (!inherits(object, "lm")) 
    warning("calling anova.lm(<fake-lm-object>) ...")
  w <- object$weights
  ssr <- sum(if (is.null(w)) object$residuals^2 else w * object$residuals^2)
  mss <- sum(if (is.null(w)) object$fitted.values^2 else w * 
               object$fitted.values^2)
  if (ssr < 1e-10 * mss) 
    warning("ANOVA F-tests on an essentially perfect fit are unreliable")
  dfr <- df.residual(object)
  p <- object$rank
  if (p > 0L) {
    p1 <- 1L:p
    comp <- object$effects[p1]
    asgn <- object$assign[stats:::qr.lm(object)$pivot][p1]
    nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
    tlabels <- nmeffects[1 + unique(asgn)]
    ss <- c(vapply(split(comp^2, asgn), sum, 1), ssr)
    df <- c(lengths(split(asgn, asgn)), dfr)
    if(reg_collapse){
      if(attr(object$terms, "intercept")){
        collapse_p<-2:(length(ss)-1)
        ss<-c(ss[1],sum(ss[collapse_p]),ss[length(ss)])
        df<-c(df[1],sum(df[collapse_p]),df[length(df)])
        tlabels<-c(tlabels[1],"Source")
      } else{
        collapse_p<-1:(length(ss)-1)
        ss<-c(sum(ss[collapse_p]),ss[length(ss)])
        df<-c(df[1],sum(df[collapse_p]),df[length(df)])
        tlabels<-c("Regression")
      }
    }
  }else {
    ss <- ssr
    df <- dfr
    tlabels <- character()
    if(reg_collapse){
      collapse_p<-1:(length(ss)-1)
      ss<-c(sum(ss[collapse_p]),ss[length(ss)])
      df<-c(df[1],sum(df[collapse_p]),df[length(df)])
    }
  }
  
  ms <- ss/df
  f <- ms/(ssr/dfr)
  P <- pf(f, df, dfr, lower.tail = FALSE)
  table <- data.frame(df, ss, ms, f, P)
  table <- rbind(table, 
                 colSums(table))
  table$ms[nrow(table)]<-table$ss[nrow(table)]/table$df[nrow(table)]
  table[length(P):(length(P)+1), 4:5] <- NA
  dimnames(table) <- list(c(tlabels, "Error","Total"), 
                          c("Df","SS", "MS", "F", 
                            "P"))
  if (attr(object$terms, "intercept")){
    table <- table[-1, ]
    table$MS[nrow(table)]<-table$MS[nrow(table)]*(table$Df[nrow(table)])/(table$Df[nrow(table)]-1)
    table$Df[nrow(table)]<-table$Df[nrow(table)]-1
  }
  structure(table, heading = c("Analysis of Variance Table\n"), 
            class = c("anova", "data.frame"))
}

youtube= c(20.64,200.16,179.64,45.84,113.04,212.40,340.32,278.52)
facebook=c(4.92,50.40,42.72,4.44,5.88,11.16,50.40,10.32)
newspaper=c(37.92,4.32,7.20,16.56,9.72,7.68,79.44,10.44)
sales=c(7.08,23.52,20.76,9.12,11.64,15.36,30.60,16.08)

fit_model=lm(sales ~ youtube + facebook + newspaper)

Anova_Table(fit_model)

#(iii)
#calculate betas
X_matrix=cbind(c(1,1,1,1,1,1,1,1),
               c(20.64,200.16,179.64,45.84,113.04,212.40,340.32,278.52),
               c(4.92,50.40,42.72,4.44,5.88,11.16,50.40,10.32),
               c(37.92,4.32,7.20,16.56,9.72,7.68,79.44,10.44))
Y_vector=cbind(c(7.08,23.52,20.76,9.12,11.64,15.36,30.60,16.08))
y_trans_x = t(Y_vector)%*%X_matrix
x_trans_y = t(y_trans_x)
x_trans_y
x_trans_x = t(X_matrix)%*%X_matrix
x_trans_x
x_trans_x_Cofactors=function(M){
  stopifnot(length(unique(dim(M)))==1)
  co_factor=M
  for(i in 1:dim(M)[1]){
    for(j in 1:dim(M)[2]){
      co_factor[i,j]=(det(M[-i,-j])*(-1)^(i+j))
    }
  }
  return(co_factor)
}

# Co-Factors
x_trans_x_Cofactors(x_trans_x)

#Adjoint Matrix 
t(x_trans_x_Cofactors(x_trans_x))

det_x_trans_x = det(x_trans_x)
det_x_trans_x

#inverse of X_TRANS_X
x_trans_x_inverse=(1/det_x_trans_x)*t(x_trans_x_Cofactors(x_trans_x))
x_trans_x_inverse
Beta_hat= x_trans_x_inverse%*%x_trans_y # OR summary(fit_model)
Beta_hat
# estimate of Sigma Square is S squared which is equal to MSE
S_squared = MSE
S_squared

#since B_o is not the co-efficient we will check the significance of B_1,B_2 and B_3
Beta_hat_1 = 0.03223192
Beta_hat_2 = 0.22972964
Beta_hat_3 = 0.02868940

t_crit = 2.776 # at alpha=0.05/2 = 0.025 and n-p degrees of freedom = 8-4=4
#Hypothesis 1
# H0 : Beta_1 = 0 vs H1 : Beta_1 != 0
#using test statistic t
t_1 = (Beta_hat_1-0)/sqrt(MSE*1.914941e-05)
t_1
#The null hypothesis H0 : Beta_1 = 0 is reject since
## |t_1| = 7.084912 > t_crit = 2.776
#and conclude that variable Youtube contributes significantly to the model

# Hypothesis 2
# H0 : Beta_2 = 0 vs H1 : Beta_2 != 0
# using test statistic t
t_2 = (Beta_hat_2-0)/sqrt(MSE*5.060408e-04)
t_2
# The null hypothesis H0 : Beta_2 = 0 is reject since
## |t_1| = 9.823138 > t_crit = 2.776
# and conclude that variable Facebook contributes significantly to the model

#Hypothesis 3
# H0 : Beta_3 = 0 vs H1 : Beta_3 != 0
#using test statistic t
t_3 = (Beta_hat_3-0)/sqrt(MSE*2.465407e-04)
t_3
#We Fail to reject the null hypothes H0 : Beta_3 = 0 since
## |t_1| = 1.757532 !> t_crit = 2.776
# we then conclude that variable Newspaper is Zero

#(iv)
# CI for regression coefficients
youtube= c(20.64,200.16,179.64,45.84,113.04,212.40,340.32,278.52)
facebook=c(4.92,50.40,42.72,4.44,5.88,11.16,50.40,10.32)
newspaper=c(37.92,4.32,7.20,16.56,9.72,7.68,79.44,10.44)
sales=c(7.08,23.52,20.76,9.12,11.64,15.36,30.60,16.08)

fit_model=lm(sales ~ youtube + facebook + newspaper)

confint(fit_model,'youtube', level=0.95)
confint(fit_model,'facebook', level=0.95)
confint(fit_model,'newspaper', level=0.95)

#(v)
Anova_Table = function (object, reg_collapse=TRUE,...) 
{
  if (length(list(object, ...)) > 1L) 
    return(anova.lmlist(object, ...))
  if (!inherits(object, "lm")) 
    warning("calling anova.lm(<fake-lm-object>) ...")
  w <- object$weights
  ssr <- sum(if (is.null(w)) object$residuals^2 else w * object$residuals^2)
  mss <- sum(if (is.null(w)) object$fitted.values^2 else w * 
               object$fitted.values^2)
  if (ssr < 1e-10 * mss) 
    warning("ANOVA F-tests on an essentially perfect fit are unreliable")
  dfr <- df.residual(object)
  p <- object$rank
  if (p > 0L) {
    p1 <- 1L:p
    comp <- object$effects[p1]
    asgn <- object$assign[stats:::qr.lm(object)$pivot][p1]
    nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
    tlabels <- nmeffects[1 + unique(asgn)]
    ss <- c(vapply(split(comp^2, asgn), sum, 1), ssr)
    df <- c(lengths(split(asgn, asgn)), dfr)
    if(reg_collapse){
      if(attr(object$terms, "intercept")){
        collapse_p<-2:(length(ss)-1)
        ss<-c(ss[1],sum(ss[collapse_p]),ss[length(ss)])
        df<-c(df[1],sum(df[collapse_p]),df[length(df)])
        tlabels<-c(tlabels[1],"Source")
      } else{
        collapse_p<-1:(length(ss)-1)
        ss<-c(sum(ss[collapse_p]),ss[length(ss)])
        df<-c(df[1],sum(df[collapse_p]),df[length(df)])
        tlabels<-c("Regression")
      }
    }
  }else {
    ss <- ssr
    df <- dfr
    tlabels <- character()
    if(reg_collapse){
      collapse_p<-1:(length(ss)-1)
      ss<-c(sum(ss[collapse_p]),ss[length(ss)])
      df<-c(df[1],sum(df[collapse_p]),df[length(df)])
    }
  }
  
  ms <- ss/df
  f <- ms/(ssr/dfr)
  P <- pf(f, df, dfr, lower.tail = FALSE)
  table <- data.frame(df, ss, ms, f, P)
  table <- rbind(table, 
                 colSums(table))
  table$ms[nrow(table)]<-table$ss[nrow(table)]/table$df[nrow(table)]
  table[length(P):(length(P)+1), 4:5] <- NA
  dimnames(table) <- list(c(tlabels, "Error","Total"), 
                          c("Df","SS", "MS", "F", 
                            "P"))
  if (attr(object$terms, "intercept")){
    table <- table[-1, ]
    table$MS[nrow(table)]<-table$MS[nrow(table)]*(table$Df[nrow(table)])/(table$Df[nrow(table)]-1)
    table$Df[nrow(table)]<-table$Df[nrow(table)]-1
  }
  structure(table, heading = c("Analysis of Variance Table\n"), 
            class = c("anova", "data.frame"))
}
youtube= c(20.64,200.16,179.64,45.84,113.04,212.40,340.32,278.52)
sales=c(7.08,23.52,20.76,9.12,11.64,15.36,30.60,16.08)
fit_model_1=lm(sales ~ youtube)
Anova_Table(fit_model_1)

SSE_r = 119.57
SSE_ur = 4.323227
degrees_Of_freedom = (8-2) - (8-4)
degrees_Of_freedom

F_cal = ((SSE_r - SSE_ur)/degrees_Of_freedom)/(SSE_ur/4)
F_cal

F_crit = 6.94 #at alpha = 0.05 2 degrees of freedom N and 4 degrees of freedom D

# the Null hypothesis H0 : Beta_2 = Beta_3 = 0
#Since F_cal = 53.31516 > F_cri = 6.94, we Reject the Null hypothesis
#and conclude that atleast one of the parameters in the unrestricted model is not Zero
#and atleast one of the variables contributes significantly to the model.

###########################################################################################