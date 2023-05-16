//////////////////////////////////////////////////
// Programa: dist_kuma.ox
//
// Estimação por máxima verossimilhança os 
// os parâmetros da distribuição kumaraswamy
//
// Autor: Kleber Henrique
//
// Versão: 1.0.0
//
// Data da última modificação: 06/07/2022
//
//////////////////////////////////////////////////

#include <oxstd.h>
#include <oxprob.h>
#import <maximize>

/* Variáveis globais */
static decl y; // vetor da amostra aleatória

const decl nobs = 20; // número de observações
const decl nrep = 10000; // número de réplicas de Monte Carlo

/* Função Log-verossimilhança com gradiente analítico */
fLog(const vtheta, const adFunc, const avScore, const amHess)
{
      decl alpha = vtheta[0]; // alpha
      decl beta = vtheta[1]; // beta

      /* Função Log-verossimilhança */
      adFunc[0] = double(nobs*log(alpha*beta) + (alpha-1)*sumc(log(y)) + (beta-1)*sumc(log(1-(y.^alpha))));

     /* Gradiente analítico */
     if(avScore)
     {
         (avScore[0])[0] = (nobs/alpha) + sumc(log(y)) + (1-beta)*sumc(((y.^alpha)'*log(y))./(1-y.^alpha));
         (avScore[0])[1] = (nobs/beta) + sumc(log(1-(y.^alpha)));
     }

	 if( isnan(adFunc[0]) || isdotinf(adFunc[0]) )
		return 0;  // falha
	 else
		return 1; // sucesso
}

main()
{
     decl vtheta, dfunc; // vetores para armazenar os valores dos parâmetros e da função log-verossimilhança, respctivamente.
     decl alpha, beta; // parâmetros alpha e beta, respectivamente.
     decl opt, i; // vetores para a otimização e para o contador do laço, respectivamente.
     decl z; // vetor para o armazenamento de números gerados pela distribuição iniforme.
     decl cFailure = 0; // vetor de contagem de falhas;
     decl valpha, vbeta; // vetores para armazenamento dos valores estimados em cada interação de MC.
     decl qNorm; // quantil da distribuição normal.
     decl alpha2 = 0.1; // nível de significância adotado.
     decl InfICalpha, SupICalpha; // vetores dos limites inferior e superior para alpha, respectivamente.
     decl InfICbeta, SupICbeta; // vetores dos limites inferior e superior para beta, respectivamente.
	 decl AmICalpha, AmICbeta; // vetores da amplitude média dos IC para alpha e beta, respectivamente.
	 decl hessi, stderror; // vetores para a matriz da hessiana numérica e erro padrão, respectivamente.
     decl contalpha = 0, contbeta = 0; // contadores para alpha e beta, se pertencem ao IC.
	 decl exectime; 

	 /* Medidas calculadas */
	 decl mediaalpha, mediabeta, varalpha, varbeta, sdalpha, sdbeta, viesalpha, viesbeta;
	 decl viesrelativoalpha, viesrelativobeta, eqmalpha, eqmbeta, tcalpha, tcbeta;
	 decl amalpha, ambeta, LIalpha, LIbeta, LSalpha, LSbeta;

	 /* Iniciando o tempo para execução do programa */
	 exectime = timer();
	 
     /* Semente */
     ranseed("GM");

     /* Parâmetros verdadeiros */
     alpha = 2;
     beta = 2;

     /* Vetores de zeros para armazenar as estimativas dos parâmetros em cada interação de MC */
     valpha = zeros(nrep,1);
     vbeta = zeros(nrep,1);

	 /* Vetores de zeros para armazenar as estimativas dos IC para os parâmetros em cada interação de MC */
	 InfICalpha = zeros(nrep,1);
     SupICalpha = zeros(nrep,1);
	 InfICbeta = zeros(nrep,1);
     SupICbeta = zeros(nrep,1);
	 AmICalpha = zeros(nrep,1);
	 AmICbeta = zeros(nrep,1);
	 
    /* Controlador de iterações */
    MaxControl(50, -1);

    /* Inicio da simulação de MC */
    for(i = 0; i<nrep ; i++)
    {

      z = ranu(nobs, 1); // gerando uniforme
      y = (1-z.^(1/beta)).^(1/alpha); // transformando em kumaraswamy
      vtheta = <3; 3>; // chute inicual

	  /* Otimização não-linear via BFGS */
      opt = MaxBFGS(fLog, &vtheta, &dfunc, 0, TRUE);
	  
      if(opt == MAX_CONV || opt == MAX_WEAK_CONV)  // Verificando a convergência
      {
         valpha[i] = vtheta[0];
         vbeta[i] = vtheta[1];

		 Num2Derivative(fLog, vtheta, &hessi);
         stderror = sqrt(diagonal(invertsym(-hessi)));
         qNorm = quann(1-alpha2/2);

         InfICalpha[i] = valpha[i]-qNorm*stderror[0];
         SupICalpha[i] = valpha[i]+qNorm*stderror[0];

		 AmICalpha[i] = SupICalpha[i]-InfICalpha[i];
		 
		 if(alpha > InfICalpha[i] && alpha < SupICalpha[i])
            contalpha++;

         InfICbeta[i] = vbeta[i]-qNorm*stderror[1];
         SupICbeta[i] = vbeta[i]+qNorm*stderror[1];

		 AmICbeta[i] = SupICbeta[i]-InfICbeta[i];

         if(beta > InfICbeta[i] && beta < SupICbeta[i])
         contbeta++;
         }
         else
         {
           i--;
           cFailure++;
         }
}  /* Fim da simulação de MC */

	
	 /* Medidas calculadas */
	 mediaalpha = double(meanc(valpha));
	 mediabeta = double(meanc(vbeta));
	 varalpha = double(meanc(valpha));
	 varbeta = double(varc(vbeta));
	 sdalpha = sqrt(varalpha);
	 sdbeta = sqrt(varbeta);
	 viesalpha = double(fabs(meanc(valpha)-alpha));
	 viesbeta = double(fabs(meanc(vbeta)-beta));
	 viesrelativoalpha = fabs((double(meanc(valpha))-alpha))/alpha*100;
	 viesrelativobeta = fabs((double(meanc(vbeta))-beta))/beta*100;
	 eqmalpha = double(varc(valpha)+fabs(meanc(valpha)-alpha)^2);
	 eqmbeta = double(varc(vbeta)+fabs(meanc(vbeta)-beta)^2);
	 tcalpha = 100*contalpha/nrep;
	 tcbeta = 100*contbeta/nrep;
	 LIalpha = double(meanc(InfICalpha));
	 LIbeta = double(meanc(InfICbeta));	
	 LSalpha = double(meanc(SupICalpha));
	 LSbeta =  double(meanc(SupICbeta));
	 amalpha = double(meanc(AmICalpha));
	 ambeta = double(meanc(AmICbeta));
	 
	 

     println(" PRINCIPAIS RESULTADOS");
	 println("\n");
     println(" Número de observações : ", nobs);
     println(" Número de Réplicas de MC : ", nrep);
	 println(" Número de falhas : ", cFailure);
	 println(" Nível de significância : ", alpha2);
     println("\n");
	 println(" Status de convergência: ", MaxConvergenceMsg(opt));
	 println(" Maximização da log-verossimilhança: ", dfunc);
	 println("\n");
	 println(" Valor verdadeiro do parâmetro alpha = ", alpha);
	 println(" Valor verdadeiro do parâmetro beta = ", beta);
	 println("\n");
	 println(" Média dos estimadores de alpha: ", mediaalpha);
     println(" Média dos estimadores de beta : ", mediabeta);
     println("\n");
     println(" Variância dos estimadores de alpha: ", varalpha);
     println(" Variância dos estimadores de beta : ", varbeta);
     println("\n");
	 println(" Desvio dos estimadores de alpha: ", sdalpha);
     println(" Desvio dos estimadores de beta : ", sdbeta);
     println("\n");
	 println(" Assimetria dos estimadores de alpha: ", moments(valpha)[3]);
     println(" Assimetria dos estimadores de beta : ", moments(vbeta)[3]);
     println("\n");
	 println(" Curtose dos estimadores de alpha: ", moments(valpha)[4]);
     println(" Curtose dos estimadores de beta : ", moments(vbeta)[3]);
     println("\n");
	 println(" Viés para alpha: ", viesalpha);
     println(" Viés para beta : ", viesbeta);
     println("\n");
	 println(" Viés Relativo para alpha: (%) ", viesrelativoalpha);
     println(" Viés Relativo para beta : (%) ", viesrelativobeta);
     println("\n");
	 println(" EQM para alpha: ", eqmalpha);
     println(" EQM para beta : ", eqmbeta);
     println("\n");
	 println(" IC assintótico para, alpha 95% cobertura: ", "%6.3f",tcalpha);
     println(" IC assintótico beta, 95% cobertura: ", "%6.3f",tcbeta);
     println("\n");
	 println(" Amplitude média do IC assintótico para, alpha 95% cobertura: ", amalpha);
     println(" Amplitude média do IC assintótico para beta, 95% cobertura: ", ambeta);
     println("\n");
	 println("Tempo de execução: ", "%6.3f",timespan(exectime));
	 
}