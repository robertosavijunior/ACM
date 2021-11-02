#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define MAX 4500

struct estrutura{
float time,d8,d9,d10,d11,d12,d13,d14,ri,visco;		
};


int main (){

	int j,i;
	float *ponteiro,lambdanm,vtolueno,k,c,nref,dndc,Na,lambda,teste;
	float normal[7],vsolvente[7],qquad[7],q[7],seno[7],coefang[MAX],coeflin[MAX],Mw[MAX],Rg[MAX],ri[MAX],ritemp[MAX],recipraizri[MAX];
	float MwGPC, A2[MAX], Lp[MAX], RgLpchute, Lpcerto[MAX], L, Mw90[MAX],viscorel[MAX],viscoesp[MAX],viscointri[MAX],RgFlory[MAX];
	struct estrutura dados[MAX];
	struct estrutura rtheta[MAX];
	struct estrutura kcrtheta[MAX];
	

	nref = 1.335;// indice de refracao
	dndc = 0.125;// sem unidade
	Na = 6.02214129e+23;// 1 mol -> numero de avogadro
	lambda = 0.00006328;// cm
	lambdanm = 632.8;
	vtolueno = 0.284386;// volts
	c = 0.0005;// g/cm3
	normal[0] = 0.765;
	normal[1] = 0.8;
	normal[2] = 1.033;
	normal[3] = 1.;
	normal[4] = 1.233;
	normal[5] = 1.086;
	normal[6] = 1.147;
	MwGPC = 1192000.;//Mw determinado pelo GPC.
	L = ((MwGPC * 1.04)/954.78526);

	//calculo das constantes
	k=(4*(M_PI*M_PI)*(nref*nref)*(dndc*dndc))/(Na*(lambda*lambda*lambda*lambda));// constante optica

	seno[0]=sin(((60/2)*M_PI)/180);
	seno[1]=sin(((69/2)*M_PI)/180);
	seno[2]=sin(((80/2)*M_PI)/180);
	seno[3]=sin(((90/2)*M_PI)/180);
	seno[4]=sin(((100/2)*M_PI)/180);
	seno[5]=sin(((111/2)*M_PI)/180);
	seno[6]=sin(((121/2)*M_PI)/180);
//	printf("valores dos senos de theta/2%f %f %f %f %f %f %f \n\n\n\n\n\n",seno[0],seno[1],seno[2],seno[3],seno[4],seno[5],seno[6]);	
	
	q[0]=((4*M_PI*nref)/lambdanm)*seno[0];
	q[1]=((4*M_PI*nref)/lambdanm)*seno[1];	
	q[2]=((4*M_PI*nref)/lambdanm)*seno[2];
	q[3]=((4*M_PI*nref)/lambdanm)*seno[3];
	q[4]=((4*M_PI*nref)/lambdanm)*seno[4];
	q[5]=((4*M_PI*nref)/lambdanm)*seno[5];
	q[6]=((4*M_PI*nref)/lambdanm)*seno[6];
	
	qquad[0]=q[0]*q[0];
	qquad[1]=q[1]*q[1];
	qquad[2]=q[2]*q[2];
	qquad[3]=q[3]*q[3];
	qquad[4]=q[4]*q[4];
	qquad[5]=q[5]*q[5];
	qquad[6]=q[6]*q[6];

	printf("q  60°       - 69°       - 80°       - 90°       - 100°       - 111°       - 121°\n");
	printf("   %f  - %f  - %f  - %f  - %f   - %f   - %f\n",q[0],q[1],q[2],q[3],q[4],q[5],q[6]);
	printf("q² 60º       - 69°       - 80°       - 90°       - 100°       - 111°       - 121°\n");	
	printf("   %f  - %f  - %f  - %f  - %f   - %f   - %f\n",qquad[0],qquad[1],qquad[2],qquad[3],qquad[4],qquad[5],qquad[6]);
	printf("k = %.10f \n",k);
	printf("c = %f g/cm³\n",c);
	getchar();

	FILE *fp;
		fp = fopen ("ACM.csv","r");//ACMII.CSV OU ACM.csv
		if (fp == NULL) {
			printf ("Não foi possivel abrir o arquivo ACM-IIA.TXT\n");
			return 1;
		}

	for(i=0;i<MAX;i++){      
	fscanf(fp,"%f;%f;%f;%f;%f;%f;%f;%f;%f;%f", &dados[i].time, &dados[i].d8, &dados[i].d9, &dados[i].d10, &dados[i].d11, &dados[i].d12, &dados[i].d13, &dados[i].d14, &dados[i].ri, &dados[i].visco);
	//printf("%f;%f;%f;%f;%f;%f;%f;%f;%f;%f\n",dados[i].time, dados[i].d8, dados[i].d9, dados[i].d10, dados[i].d11, dados[i].d12, dados[i].d13, dados[i].d14, dados[i].ri, dados[i].visco);
	}
	vsolvente[0] = dados[50].d8;
	vsolvente[1] = dados[50].d9;
	vsolvente[2] = dados[50].d10;
	vsolvente[3] = dados[50].d11;
	vsolvente[4] = dados[50].d12;
	vsolvente[5] = dados[50].d13;
	vsolvente[6] = dados[50].d14;
		
			
	for(i=0;i<MAX;i++){
		
		//formula de rtheta		
		rtheta[i].d8=((dados[i].d8-vsolvente[0])/vtolueno)*0.944*(2258203.827/(lambdanm*lambdanm*lambdanm*lambdanm))*normal[0];
		rtheta[i].d9=((dados[i].d9-vsolvente[1])/vtolueno)*0.944*(2258203.827/(lambdanm*lambdanm*lambdanm*lambdanm))*normal[1];
		rtheta[i].d10=((dados[i].d10-vsolvente[2])/vtolueno)*0.944*(2258203.827/(lambdanm*lambdanm*lambdanm*lambdanm))*normal[2];
		rtheta[i].d11=((dados[i].d11-vsolvente[3])/vtolueno)*0.944*(2258203.827/(lambdanm*lambdanm*lambdanm*lambdanm))*normal[3];
		rtheta[i].d12=((dados[i].d12-vsolvente[4])/vtolueno)*0.944*(2258203.827/(lambdanm*lambdanm*lambdanm*lambdanm))*normal[4];
		rtheta[i].d13=((dados[i].d13-vsolvente[5])/vtolueno)*0.944*(2258203.827/(lambdanm*lambdanm*lambdanm*lambdanm))*normal[5];
		rtheta[i].d14=((dados[i].d14-vsolvente[6])/vtolueno)*0.944*(2258203.827/(lambdanm*lambdanm*lambdanm*lambdanm))*normal[6];
	}

	for(i=0;i<MAX;i++){
		//printf("%.10f ; %.10f ; %.10f ; %.10f ; %.10f ; %.10f ; %.10f\n",rtheta[i].d8,rtheta[i].d9,rtheta[i].d10,rtheta[i].d11,rtheta[i].d12,rtheta[i].d13,rtheta[i].d14);
		
		kcrtheta[i].d8  = (k * c) / rtheta[i].d8;
		kcrtheta[i].d9  = (k * c) / rtheta[i].d9;
		kcrtheta[i].d10 = (k * c) / rtheta[i].d10;
		kcrtheta[i].d11 = (k * c) / rtheta[i].d11;
		kcrtheta[i].d12 = (k * c) / rtheta[i].d12;
		kcrtheta[i].d13 = (k * c) / rtheta[i].d13;
		kcrtheta[i].d14 = (k * c) / rtheta[i].d14;
 
	/*  if (((i+1)%1000)==0){
			printf("Tecle para continuar mostrando\n");			
			getchar();
		}
	*/	
	}

	for(i=0;i<MAX;i++){
	
		coefang[i] = ( (7 * ((qquad[0] * kcrtheta[i].d8) + (qquad[1] * kcrtheta[i].d9) + (qquad[2] * kcrtheta[i].d10) + (qquad[3] * kcrtheta[i].d11)+ (qquad[4] * kcrtheta[i].d12) + (qquad[5] * kcrtheta[i].d13) + (qquad[6] * kcrtheta[i].d14)) - (qquad[0] + qquad[1] + qquad[2] + qquad[3] + qquad[4] + qquad[5] + qquad[6]) * (kcrtheta[i].d8 + kcrtheta[i].d9 + kcrtheta[i].d10 + kcrtheta[i].d11 + kcrtheta[i].d12 + kcrtheta[i].d13 + kcrtheta[i].d14)) / (7 * ((qquad[0] * qquad[0]) + (qquad[1] * qquad[1]) + (qquad[2] * qquad[2]) + (qquad[3] * qquad[3]) + (qquad[4] * qquad[4]) + (qquad[5] * qquad[5]) + (qquad[6] * qquad[6])) - ((qquad[0] + qquad[1] + qquad[2] + qquad[3] + qquad[4] + qquad[5] + qquad[6]) * (qquad[0] + qquad[1] + qquad[2] + qquad[3] + qquad[4] + qquad[5] + qquad[6]))) );

		coeflin[i] = ( ((kcrtheta[i].d8 + kcrtheta[i].d9 + kcrtheta[i].d10 + kcrtheta[i].d11 + kcrtheta[i].d12 + kcrtheta[i].d13 + kcrtheta[i].d14) / 7) - (coefang[i] * ((qquad[0] + qquad[1] + qquad[2] + qquad[3] + qquad[4] + qquad[5] + qquad[6]) / 7)));

		//printf("\n %d ---> a = %.10f ; b = %.10f \n",i+1,coefang[i], coeflin[i]);
	}

	for(i=0;i<MAX;i++){
	
		Mw[i] = (1/coeflin[i]);
		Rg[i] = (coefang[i] * Mw[i] * 3*(-1));
		Rg[i] = sqrt(Rg[i]);		
		printf("%d ---> Mw = %.10f",i+1,Mw[i]);
		printf(" - Rg² = %.10f \n", Rg[i]);

	/*	if (((i+1)%100)==0){
			printf("Tecle para continuar mostrando\n");			
			getchar();
		}
	*/
	}
 
	for(i=0;i<1668;i++){//esse numero Ã© onde a sensibilidade do detector de ri estoura      
		ritemp[i]=(dados[i].ri-dados[50].ri);
		ritemp[i]=((ritemp[i]*0.000292908)/0.180);//regra de 3 com o (FC 0.000585816 ACMII) (FC 0.000292908 ACM). 	
		ritemp[i]=(ritemp[i]*1000);//convertendo para mg/ml.
		ri[i]=(ritemp[i]/58.44);//convertendo para mol/L.
		recipraizri[i]=(1/(sqrt(ri[i])));//reciproco da raiz quadrada de ri.
	}
	for(i=1668;i<MAX;i++){//esse numero Ã© onde a sensibilidade do detector de ri estoura      
		ritemp[i]=(dados[i].ri-dados[50].ri);
		ritemp[i]=((ritemp[i]*0.000585816)/0.180);//regra de 3 com o (FC 0.001171632 ACMII) (FC 0.000585816 ACM). 
		ritemp[i]=(ritemp[i]*1000);//convertendo para mg/ml.
		ri[i]=(ritemp[i]/58.44);//convertendo para mol/L.
		recipraizri[i]=(1/(sqrt(ri[i])));//reciproco da raiz quadrada de ri.
	}
	for(i=0;i<MAX;i++){
	//	A2[i]=(coeflin[i]-(1/MwGPC))/2*c;
		A2[i]=(1/(2*c))*((kcrtheta[i].d11-(1/MwGPC)));
		printf("[NaCl] %f = A2 %.30f \n",ri[i],A2[i]);
	/*	if (((i+1)%100)==0){
			printf("Tecle para continuar mostrando\n");			
			getchar();
		}
	*/
	}

	for(i=0;i<MAX;i++){
		for(j=0;j<200;j++){
			RgLpchute = (((L*j)/3)-(j*j)+(2*(j*j*j)/L)-(2*((j*j*j*j)/(L*L)))*(1-exp(-L/j)));			
			RgLpchute = sqrt(RgLpchute);
			if ((RgLpchute >= (Rg[i] - (Rg[i]*0.05))) && (RgLpchute <= (Rg[i] + (Rg[i]*0.05)))){
				Lp[i]=RgLpchute;
			}
		}
	}
	
	for(i=0;i<MAX;i++){
		printf("[NaCl] %f = Lp %f \n",ri[i],Lp[i]);
		/*
		if (((i+1)%100)==0){
			printf("Tecle para continuar mostrando\n");			
			getchar();
		}
		*/
	}

	FILE *gravar;
		gravar = fopen ("RESULTADO.TXT","w");//RESULTADO.TXT E RESULTADO2.TXT
		if (gravar == NULL) {
			printf ("Nao foi possivel abrir o arquivo RESULTADO.TXT\n");
			return 1;
		}
	fprintf(gravar,"time ; Mw ; Rg ; ri ; 1/sqrt(ri) ; A2 ; Lp\n");
	for(i=0;i<MAX;i++){      
		if(Lp[i]<=10)		
			fprintf(gravar,"%f ; %f ; %f ; %f ; %f ; %.20f ;  \n",dados[i].time, Mw[i], Rg[i], ri[i], recipraizri[i], A2[i]);
		else
			fprintf(gravar,"%f ; %f ; %f ; %f ; %f ; %.20f ; %f \n",dados[i].time, Mw[i], Rg[i], ri[i], recipraizri[i], A2[i], Lp[i]);
		
	}
	
	for(i=0;i<MAX;i++){
		Mw90[i]= (1 / kcrtheta[i].d11);
		printf("Mw para 90 = %f \n",Mw90[i]);
	}
	FILE *gravamw90;
		gravamw90 = fopen ("Mw.txt","w");
		if (gravamw90 == NULL){
			printf("nao foi possivel criar o arquivo Mw.txt");
			return 1;
		}
	fprintf (gravamw90,"[NaCl] ; Mw90 \n");
	for (i=0;i<MAX;i++){
		fprintf(gravamw90,"%f ; %f \n",ri[i],Mw90[i]);
	}

	fclose (fp);
	fclose (gravar);
	fclose (gravamw90);

	FILE *pkcrtheta;
		pkcrtheta = fopen ("kcrtheta.txt","w");
		if (pkcrtheta == NULL){
			printf("nao foi possivel criar o arquivo kcrtheta.txt");
			return 1;
		}
	fprintf (pkcrtheta,"[NaCl] ; 60 ; 69 ; 80 ; 90 ; 100 ; 111 ; 121 \n");
	for (i=0;i<MAX;i++){
		fprintf(pkcrtheta,"%f ; %.20f  ; %.20f ; %.20f ; %.20f ; %.20f ; %.20f ; %.20f \n",ri[i],kcrtheta[i].d8,kcrtheta[i].d9,kcrtheta[i].d10,kcrtheta[i].d11,kcrtheta[i].d12,kcrtheta[i].d13,kcrtheta[i].d14);
	}

	for(i=0;i<MAX;i++){
		viscorel[i]=((dados[i].visco-dados[4500].visco)/(dados[50].visco-dados[4500].visco));
		viscoesp[i]=(viscorel[i]-1);
		viscointri[i]=(viscoesp[i]/c);
		RgFlory[i]=cbrt((viscointri[i]*MwGPC)/(cbrt(36)*2.56e+23));
		RgFlory[i]=RgFlory[i]*10000000;
		printf("%f \n",RgFlory[i]);
	}

	FILE *foxflory;
		foxflory = fopen ("foxflory.txt","w");
		if (foxflory == NULL){
			printf("nao foi possivel criar o arquivo foxflory.txt");
			return 1;
		}
	fprintf (foxflory,"ViscoIntrinseca ; RgFloryFox \n");
	for (i=0;i<MAX;i++){
		fprintf(foxflory,"%f ; %.10f  \n",viscointri[i],RgFlory[i]);
	}




	fclose (fp);
	fclose (gravar);
	fclose (gravamw90);
	fclose (pkcrtheta);
	fclose (foxflory);

/*
Voltagem Tolueno = "0.284386"	
d8	-	60Â°
d9	-	69Â°
d10	-	80Â°
d11	-	90Â°
d12	-	100Â°
d13	-	111Â°
d14	-	121Â°
*/

return 0;
}
