#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define MAX 4000

int main (){

	int j,i,bla;
	float *ponteiro,lambdanm,vtolueno,k,c,nref,dndc,Na,lambda,teste;
	float normal[7],vsolvente[7],qquad[7],q[7],seno[7],coefang[MAX],coeflin[MAX],Mw[MAX],Rg[MAX],ri[MAX],ritemp[MAX],recipraizri[MAX];
	float MwGPC, A2[MAX], Lp[MAX], RgLpchute, Lpcerto[MAX], L;
	nref=1.335;// índice de refração
	dndc=0.125;// sem unidade
	Na=6.02214129e+23;// 1 mol
	lambda=0.00006328;// cm
	lambdanm=632.8;
	vtolueno=0.284386;// volts
	c=0.0005;// g/cm3
	normal[0]=0.765;
	normal[1]=0.8;
	normal[2]=1.033;
	normal[3]=1.;
	normal[4]=1.233;
	normal[5]=1.086;
	normal[6]=1.147;
	MwGPC=1250000.;//Mw determinado pelo GPC.
	L=((MwGPC * 1.04)/954.78526);
	//char linha[120];

for (i=0;i<MAX;i++){
	Rg[i]=140;
}

	for(i=0;i<MAX;i++){
		for(j=0;j<200;j++){
			RgLpchute = (((L*j)/3)-(j*j)+(2*(j*j*j)/L)-(2*((j*j*j*j)/(L*L)))*(1-exp(-L/j)));			
			//RgLpchute = (RgLpchute * (-1.));
			RgLpchute = sqrt(RgLpchute);
			//printf("%f \n",RgLpchute);
			if ((RgLpchute >= (Rg[i] - (Rg[i]*0.05))) && (RgLpchute <= (Rg[i] + (Rg[i]*0.05)))){
				Lp[i]=RgLpchute;
			}
		}
	}

for(i=0;i<MAX;i++){
printf("%f \n",Lp[i]);
}
return 0;
}
