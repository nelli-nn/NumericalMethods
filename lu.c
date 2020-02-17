#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <math.h>
#define I 5
#define J 5
void IdentityMatrix(double P[][J]){
	int i, j;
	for(i = 0; i < I; i++){
		for(j = 0; j < J; j++){
		  if(i == j)
				P[i][j] = 1;
			else 
				P[i][j] = 0;
		}
	}
}
void Multiply(double  C[][J], double A[][J], double B[][J]){
  int i, j, k, count = 0, p = 0;
	double sum = 0;
	for(i = 0; i < I; i++)
		for(j = 0; j < J; j++){
			C[i][j] = 0;
			for(k = 0; k < I; k++){
				C[i][j] += A[k][i] * B[j][k];
			}
		}
}
void GetMatrix(double matrixA[][J], FILE * file){
  int i = 0, j, l = 0, t = 0;
	char num[256];
	int symb = 0;
	char get, prev;
	for(j = 0; (symb = fgetc(file)) != EOF; j++){
		get = (char)symb;
		prev = symb;
		while((symb >= '0' && symb <= '9')  || symb == '-' || symb == '.'){
			num[l++] = symb;
			symb = fgetc(file);
			num[l] = '\0';
		}
		if(symb == '\n'){
			t = 0;
		  i++;
	  } else 
  		matrixA[i][t++] = atof(num);
		l = 0;
	}
}
void GetRight(double right[], FILE * file){
  int i = 0, j, l = 0, t = 0;
	char num[256];
	int symb = 0;
	char get, prev;
	for(j = 0; (symb = fgetc(file)) != EOF; j++){
		get = (char)symb;
		prev = symb;
		while((symb >= '0' && symb <= '9')  || symb == '-' || symb == '.'){
			num[l++] = symb;
			symb = fgetc(file);
			num[l] = '\0';
		}
		if(t < J)
  	  right[t++]= atof(num);
		l = 0;
	}
}
void MultiplyVec(double  C[J], double A[][J], double B[J]){
  int i, j, k, count = 0, p = 0;
	double sum = 0;
	for(i = 0; i < I; i++)
		for(j = 0; j < J; j++){
			C[j] = 0;
			for(k = 0; k < I; k++)
				C[j] += A[j][k] * B[k];
		}
}
void SolveLz(double z[], double matrixL[][J], double right[]){
  int i, j;
	double sum = 0;
	z[0] = right[0] / matrixL[0][0]; 
	for(i = 1; i < I; i++){
		for(j = 0; j < i; j++){
			sum += matrixL[i][j] * z[j];
		}
		z[i] = 1 / matrixL[i][i] * (right[i] - sum); 
	  sum = 0;
	}
}
void SwapRows(double P[][J], double matrixA[][J], int pivot , int i){
	int k;
	double tmp;
	for(k = 0; k < I; k++){
		tmp = P[i][k];
		P[i][k]=P[pivot][k];
		P[pivot][k]=tmp;
		tmp = matrixA[i][k];
		matrixA[i][k]= matrixA[pivot][k];
		matrixA[pivot][k]=tmp;
  }
}
void SolveUx(double z[], double matrixL[][J], double right[]){
  int i, j;
	double sum = 0;
	for(i = I - 1; i >= 0; i--){
		for(j = i + 1; j < J; j++){
			sum += matrixL[i][j] * z[j];
		}
		z[i] = (right[i] - sum) / matrixL[i][i];
	  sum = 0;
	}
}


void CreateLU(double m[][J], double matrixA[][J], double matrixU[][J], double matrixL[][J],  double P[][J]){
  int i, j, t = 0, row, pivot, k;
	double C[I][J],  n[J][I];
	double sum = 0, sum2 = 0, pivotValue;
	int exchanges = 0; 
	for(i = 0; i < I; i++){
		for(j = 0; j < J; j++){
  		C[i][j] = matrixA[i][j];
			m[i][j] = matrixA[i][j];
		}
	}
	for( i = 0; i < I; i++ ) {
    pivotValue = 0;
    pivot = -1;
    for(row = i; row < I; row++){
			if(fabs(C[row][i]) > pivotValue){
				pivotValue = fabs(C[row][i]);
				pivot = row;
			}
		}
		if(pivotValue != 0) {
			SwapRows(P, C, pivot, i);
			for(j = i+1; j < I; j++) {
				C[j][i] /= C[i][i];
				for(k = i+1; k < I; k++) {
					C[j][k] -= C[j][i] * C[i][k];
				}
 			}
		}
	}
	Multiply(n, matrixA, P);
	for(i = 0; i < I; i++){
		for(j = 0; j < J; j++){
       matrixA[i][j]= n[i][j];
		}
	}
	for(i = 0; i < I; i++){
		for(j = 0; j < J; j++){
       matrixA[i][j]= n[j][i];
		}
	}
	for(j = 0; j < J; j++){
		for(i = 0; i < I; i++){
			if(!j)
			matrixU[j][i] = matrixA[j][i];
			if(!j)
				matrixL[i][j] = matrixA[i][j] / matrixU[0][0];
     	for(t = 0; t <= j - 1; t++)
				sum += matrixL[j][t] * matrixU[t][i];
			for(t = 0; t <= i - 1; t++)
				sum2 += matrixL[j][t] * matrixU[t][i];
			if(i && i >= j)
        matrixU[j][i] = matrixA[j][i] - sum;
			else if(j){
				matrixL[i][j] = 0;
				matrixU[j][i] = 0;
			}
			if(i && j > i){
				sum2 = matrixA[j][i] - sum2;
				matrixL[j][i] = (sum2) / matrixU[i][i];
			}
		  else if(i == j){
				matrixL[j][i] = 1;
			} 
			sum = 0;
		  sum2 = 0;
		}
	} 
}


int main(void){
	int i, j;
	int count = 0;
	double y[I][J];
	double A[J];
	double matrixA[I][J];
	double m[I][J];
	double matrixU[I][J];
	double matrixL[I][J];
	double right[I];
	double z[I];
	double x[I];
	FILE * file, * myfile;
	errno_t err;

	if((err = fopen_s(&file, "xx.txt","r")) == 0){
	  GetMatrix(matrixA, file);
	}
	if((err = fopen_s(&myfile, "1x.txt", "r")) == 0){
	  GetRight(right, myfile);
	}
	fclose(file);
	fclose(myfile);
	IdentityMatrix(y);
	CreateLU(m, matrixA, matrixU, matrixL, y);
	for(i = 0; i < I; i++){
		for(j = 0; j < J; j++){
  		matrixA[i][j] =	m[i][j];
		}
	}

  MultiplyVec(A, y, right);
	SolveLz(z, matrixL, A);
	SolveUx(x, matrixU, z);
	return 0;
} 