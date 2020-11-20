m = 10
n = 10
maxrank = 2
p = .2
#W = matrix(sample(maxrank,m*n,replace=TRUE),nrow=m)
v1 = sample(maxrank,m,replace=TRUE)
v2 = sample(maxrank,m,replace=TRUE)

W = v1 %*% t(v1) + v2 %*% t(v2)
omega = matrix(rbinom(m*n,1,p),nrow=m)

print("W")
print(W)
print("omega")
print(omega)

P = function(Y){
  return(Y*omega)
}
Pp = function(Y){
  return(Y- P(Y))
}

S = function(lambda,W){
  decomp = svd(W)
  U = decomp$u
  V = decomp$v
  D = decomp$d
  r = length(D)
  Dlambda = diag(D) - diag(rep(lambda,r))
  return(U%*%Dlambda%*%t(V))
}

lambda = .1
eps = .0001
stop = FALSE

numreps = 10000
Zold = matrix(rep(0,m*n),nrow=m)
error = NULL
while(stop == FALSE){
  Znew = S(lambda,P(W)+Pp(Zold))
  print(norm(Znew - Zold,"F")^2)
  if(norm(Znew - Zold,"F")^2 < eps ){ #* norm(Zold,"F")^2
    stop = TRUE
  }
  else{Zold = Znew}
  error <- append(error,norm(Znew-W,"F"))
}

plot(error)
print("znew")
print(Znew)
print("difference")
print(Znew-W)
