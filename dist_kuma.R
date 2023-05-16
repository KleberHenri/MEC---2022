####################################################
# Programa: dist_kuma.R
#
# Estimação por máxima verossimilhança os 
# os parâmetros da distribuição kumaraswamy
#
# Autor: Kleber Henrique
#
# Versão: 1.0.0
#
# Data da última modificação: 06/07/2022
#
####################################################

# Pacote utilizado
library(moments)


# Função contendo a log-verossimilhança, vetor score e 
# a otimização via BFGS
Kuma.fit <- function(y){
  
  n <- length(y)  #Tamanho do vetor y
  
  # Função log-verossimilhança
  loglik_kuma <- function(theta){
    alpha <-  theta[1]
    beta <-  theta[2]
    
    loglik <- n*log(alpha*beta) + (alpha-1)*sum(log(y)) + (beta-1)*sum(log(1-(y^alpha)))
    loglik
   
  }
  
  # Vetor score
  scoreFn <- function(theta){
    alpha <-  theta[1]
    beta <-  theta[2]
    
    Ualpha <- (n/alpha)+sum(log(y))+(1-beta)*sum(((y^alpha)*log(y))/(1-(y^alpha)))
    Ubeta <- (n/beta) + sum(log(1-(y^alpha)))
    
    derivada <- c(Ualpha, Ubeta)
    derivada
    
  }
  
  # Chutes inicias
  ini <- c(3,3)
  
  # Maximização via BFGS
  opt = optim(ini, fn = loglik_kuma, gr = scoreFn, method = "BFGS", hessian = TRUE
             control = list(fnscale = -1, reltol = 1e-12))
  
  if (opt$conv != 0)
    warning("FUNCTION DID NOT CONVERGE!")
  
  # Adicionando elementos da otimização a uma lista para ser retornada ao 
  # final da função "kuma.fit"
  k <- c()
  coef <- opt$par
  k$coef <- coef
  k$conv <- opt$conv
  k$loglik <- opt$value
  k$counts <- as.numeric(opt$counts[1])
  
  k$alpha <- coef[1]
  k$beta <- coef[2]
  
  # Matriz de informação de Fisher 
  a_barra <- 1 + (k$beta/(k$beta-2))*((digamma(k$beta)-digamma(2))^2 -
                                      (trigamma(k$beta)-trigamma(2)))
  b_barra <- -(1/(k$beta-1))*(digamma(k$beta+1) - digamma(2))

  Ibb <- (n*a_barra)/(k$alpha^2)
  Ibs <- (n*b_barra)/k$alpha
  Isb <- t(Ibs)
  Iss <- n/(k$beta^2)

  I <- rbind(cbind(Ibb,Ibs),cbind(Isb,Iss))

  k$I <- I
  
  # Inversa da matriz de Fisher
  vcov <- try(solve(as.matrix(k$I)))

  k$vcov <- vcov
  
  # Erro padrão
  stderror <- sqrt(diag(k$vcov))
  k$stderror <- stderror
  
  k
  
}

# Simulação ---------------------------------------------------------------
R = 10000 # Replicas de MC
vn =  c(20,30,50,80,100,200,500,1000) # Tamanhos da amostra
alfa = 0.01 #  nível de significância
set.seed(2021) # semente


# Parâmetros 
alpha = 2
beta = 2

pa = c(alpha,beta)

for(n in vn)
{
  
  
  # inicializations
  TC<-c() 
  coef_results<-c()
  results<-c()
  AM <- c()
  
  nc=0 # Contador de não convergência
  
  tempo.inicio = Sys.time() # Tempo inicial da simulação de MC
  
  i = 0 # Contador da simulação de MC
  while(i<R) # Inicio da simulação de MC
  {
    z <- runif(n) # Gerando uniformes
    y <- (1-(z^(1/beta)))^(1/alpha) # Gerando kumaswamy pela inversão
    
    fit <- Kuma.fit(y) # fit da distribuição
    
    if( fit$conv==0 )
    {
      results<-rbind(results,fit$pvalues)
      coef_results<-rbind(coef_results,fit$coef)

      LI<- fit$coef - qnorm(1-alfa/2)* fit$stderror
      LS<- fit$coef + qnorm(1-alfa/2)* fit$stderror

      AM <- rbind(AM, (LS - LI))

      TC <- rbind(TC, ((pa<LI) + (pa>LS)) )

      i = i+1
    }else{ # nonconvergence failures
      nc<-nc+1
      print(c("Nonconvergence",i,nc))
    }
  }
  
  m<-colMeans((coef_results),na.rm=T)
  va <- apply(coef_results,2,var)
  sd<-apply(coef_results,2,sd)
  assi<-apply(coef_results,2,moments::skewness)
  kurt<-apply(coef_results,2,moments::kurtosis)
  
  bias<-m-pa
  rb<- 100*bias/pa
  eqm <- var(m-pa)+ (bias^2)
  tc <- 1-colMeans((TC),na.rm=T) #  coverage rate
  am <- colMeans((AM),na.rm=T)
  
  tempo.fim = Sys.time()
  
  tempo.exec = tempo.fim - tempo.inicio
  
  M<-rbind(pa,m,va,sd,assi,kurt,bias,rb,eqm,tc,am,tempo.exec)
  row.names(M)<-c("Parameters","Mean","Var","SD","Assi","Curtose","Bias","RB","EQM","TC","AM","Time")
  
  print(c("n",n),quote=F)
  print(c("nc",nc-ni),quote=F)
  print(c("ni",ni),quote=F)
  print(round(M,3))
}







