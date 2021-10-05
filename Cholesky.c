/*****************************************************************
* Projeto MS512 - Análise Numérica                               *
* Prof Maicon R. Correa                                          *
* Fatoração de Cholesky                                          *
*                                                                *
* Bruno S. T. Pagotto 154891                                     *
* João Pedro N. S. 176146                                        *
* Lucas G. Machado 172776                                        *
* Matheus L. Bernardi 184331                                     *
******************************************************************/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define N 12

void matriz_vazia(double A[][N]){
	int x,j;
	for(x=0;x<N;x++){
		for(j=0;j<N;j++)
			A[x][j]=0;
	}

	return;
}

void produto(double A[][N],double B[][N],double AB[][N]){//matriz AB deve ser nula, func�o altera AB para ser o produto A*B
	int i,j,k;
	matriz_vazia(AB);
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			for(k=0;k<N;k++)
				AB[i][j]+=A[i][k]*B[k][j];
		}
	}

	return;
}

void solveT_sup_sistem(double A[][N],double X[][N],double b[][N]){//funcao resolve um sistema triangular superior AX=b
	int l,j,k;
	double h=0;
	for(l=0;l<N;l++){
		for(j=N-1;j>=0;j--){
			for(k=j+1;k<N;k++)
				h+=A[j][k]*X[k][l];
			
			X[j][l]=(b[j][l]-h)/A[j][j];
			h=0;
		}
	}

	return;	
}

void solveT_inf_sistem(double A[][N],double X[][N],double b[][N]){//funcao resolve um sistema triangular inferior AX=b
	int l,j,k;
	double h=0;
	for(l=0;l<N;l++){
		for(j=0;j<N;j++){
			for(k=0;k<j;k++)
				h+=A[j][k]*X[k][l];
			
			X[j][l]=(b[j][l]-h)/A[j][j];
			h=0;
		}
	}

	return;
}

void transpor(double A[][N]){
	int i,j;
	double a;
	for(i=0;i<N;i++){
		for(j=i;j<N;j++){
			a=A[i][j];
			A[i][j]=A[j][i];
			A[j][i]=a;	
		}
	}	
	
	return;
}

int cholesky(double A[][N],double R[][N]){
	int i,j,k;
	double a=0,b=0,c=0;
	for(i=0;i<N;i++){
		for(j=i;j<N;j++){
			for(k=0;k<i;k++){
				a+=pow(R[k][i],2);
				b+=R[k][i]*R[k][j];
			}
			if(i==j){
			c=A[i][j]-a;
			if(c<=0){
				return -1;	
			}
			R[i][j]=sqrt(c);	
			}else if(j>i){
				R[i][j]=(A[i][j]-b)/R[i][i];
			}
			a=0;
			b=0;
			
		}
	}

	return 0;
}

double norma(double A[][N]){
	int i,j;
	double max=0,a=0,b=0;
	for(j=0;j<N;j++){
		for(i=0;i<N;i++){
			if(A[i][j]<0){
				b=-A[i][j];
			}else{
				b=A[i][j];
			}
			a+=b;
		}
		if(a>max){
			max=a;
		}
		a=0;
	}

	return max;
}

void imprimir(double C[][N]){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++)
			printf("%f ",C[i][j]);
		
		printf("\n");
	}

	return;
}

int main(){
	int i,j,k;
	double c,A[N][N],X[N][N],B[N][N],R[N][N],RT[N][N],C[N][N],Y[N][N];;
	
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			c=i+j+1;
			A[i][j]=1/c;
		}
	}
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			if(i==j){
				B[i][j]=1;
			}else{
				B[i][j]=0;
			}
		}
	}
	
	matriz_vazia(R);
	k=cholesky(A,R);
	if(k==-1){
		printf("nao e positiva definida");
		return 0;
	}
	
	matriz_vazia(RT);
	cholesky(A,RT);
	transpor(RT);


	matriz_vazia(Y);
	solveT_inf_sistem(RT,Y,B);

	matriz_vazia(X);
	solveT_sup_sistem(R,X,Y);

	matriz_vazia(C);
	produto(A,X,C);
	printf("matriz %dx%d:\n",N,N);
	printf("norma de AX=%.20f\n",norma(C));
	printf("norma de X=%.20f\n",norma(X));
	printf("norma de A=%.20f\n",norma(A));
	printf("numero de condicionamento de A=%.20f\n",norma(A)*norma(X));

	return 0;
}