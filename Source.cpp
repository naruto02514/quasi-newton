#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M 2
#define min 1.0e-5

double alpha;
double x1[M];
double x2[M];
double d[M];
double B1[M][M];
double B2[M][M];
double b[M];
double b1[M];
double H[M][M];
double H1[M][M];
int i, j, k, loop, flag, mode;

double fx1(double x, double y){
	return 2 * (1 - 200 * y)*x + 400 * x*x*x - 2;
}

double fx2(double x, double y){
	return 200 * (y - x*x);
}

double f(double x, double y){
	return 100 * ((y - x*x)*(y - x*x)) + ((1 - x)*(1 - x));
}

void initial(){
	loop = 0;
	flag = 0;
	x1[0] = 1.0;
	x1[1] = 0.8;
	B1[0][0] = 1.0;
	B1[0][1] = 0.0;
	B1[1][0] = 0.0;
	B1[1][1] = 1.0;
	b[0] = fx1(x1[0], x1[1]);
	b[1] = fx2(x1[0], x1[1]);
}


void hk(){
	if (loop == 0){
		double detA;
		detA = B1[0][0] * B1[1][1] - B1[0][1] * B1[1][0];
		H[0][0] = B1[1][1] / detA;
		H[0][1] = -B1[0][1] / detA;
		H[1][0] = -B1[1][0] / detA;
		H[1][1] = B1[0][0] / detA;
	}

	for (i = 0; i<M; i++){
		d[i] = 0.0;
		for (j = 0; j<M; j++){
			d[i] += -(H[i][j])*b[j];
		}
	}
}

void LU(){
	double u[M][M], l[M][M], c[M];
	for (i = 0; i<M; i++){  // L行列、U行列を1と0で初期化
		c[i] = 0.0;
		d[i] = 0.0;
		for (j = 0; j<M; j++){
			u[i][j] = 0.0;
			if (i == j)
				l[i][j] = 1.0;
			else
				l[i][j] = 0.0;
		}
	}
	for (i = 0; i<M; i++){
		for (j = i; j<M; j++){ // U行列の生成
			u[i][j] = B1[i][j];
			for (k = 0; k<i; k++){
				u[i][j] -= u[k][j] * l[i][k];
			}
		}
		for (j = i + 1; j<M; j++){ // L行列の生成 
			l[j][i] = B1[j][i];
			for (k = 0; k<i; k++){
				l[j][i] -= u[k][i] * l[j][k];
			}
			l[j][i] /= u[i][i];
		}
	}
	for (i = 0; i<M; i++){ // c行列の生成
		c[i] = -b[i];
		for (j = 0; j<i; j++){
			c[i] -= l[i][j] * c[j];
		}
	}
	for (i = M - 1; i >= 0; i--){ // x行列の生成
		d[i] = c[i];
		for (j = M - 1; j>i; j--){
			d[i] -= u[i][j] * d[j];
		}
		d[i] /= u[i][i];
	}
}

void alpha_update(){
	double FTD, ff;
	FTD = 0.0;
	alpha = 1.0;
	ff = f(x1[0], x1[1]); //f(x)
	for (i = 0; i<M; i++){
		FTD += b[i] * d[i]; //∇f(x)*dk
	}
	while (f(x1[0] + alpha*d[0], x1[1] + alpha*d[1]) >= ff + 0.5*alpha*FTD){
		alpha *= 0.5;
	}
}

void x2_update(){
	for (int i = 0; i<M; i++){
		x2[i] = x1[i] + alpha*d[i]; //x(k+1)=xk+α*dk
	}
}

void BGFS_HK(){
	double y[M], s[M], sTy, Hy[M], HysT[M][M], sHyT[M][M], yTHy, ssT[M][M];
	sTy = 0.0;
	yTHy = 0.0;
	b1[0] = fx1(x2[0], x2[1]);
	b1[1] = fx2(x2[0], x2[1]);
	for (i = 0; i<M; i++){
		Hy[i] = 0.0;
		s[i] = x2[i] - x1[i]; //sk=xk+1-xk
		y[i] = b1[i] - b[i]; //yk=∇f(x+1)-∇f(x)
	}
	for (i = 0; i<M; i++){
		sTy += s[i] * y[i];
		for (j = 0; j<M; j++){
			Hy[i] += H[i][j] * y[j];
			ssT[i][j] = s[i] * s[j];
		}
	}
	for (i = 0; i<M; i++){
		yTHy += y[i] * Hy[i];
		for (j = 0; j<M; j++){
			HysT[i][j] = Hy[i] * s[j];
			sHyT[i][j] = s[i] * Hy[j];
		}
	}

	for (i = 0; i<M; i++){
		for (j = 0; j<M; j++){
			H1[i][j] = H[i][j] - ((HysT[i][j] + sHyT[i][j]) / sTy) + (1 + (yTHy / sTy))*(ssT[i][j] / sTy);
		}
	}
}

void BFGS(){
	//データ初期化
	double sTy, sTBs, BS[M], BST[M], sT[M], y[M], s[M], yyT[M][M], BSBST[M][M];
	sTy = 0.0;
	sTBs = 0.0;
	b1[0] = fx1(x2[0], x2[1]);
	b1[1] = fx2(x2[0], x2[1]);
	for (i = 0; i<M; i++){ //データ初期化
		BS[i] = 0.0;
		BST[i] = 0.0;
		sT[i] = 0.0;
	}
	for (i = 0; i<M; i++){
		s[i] = x2[i] - x1[i]; //sk=xk+1-xk
		y[i] = b1[i] - b[i]; //yk=∇f(x+1)-∇f(x)
	}
	for (i = 0; i<M; i++){
		for (j = 0; j<M; j++){
			BS[i] += B1[i][j] * s[j]; //Bk*sk
			yyT[i][j] = y[i] * y[j]; //Yk*YkT
		}
	}
	for (i = 0; i<M; i++){
		sTBs += s[i] * BS[i];
		for (j = 0; j<M; j++){
			BSBST[i][j] = BS[i] * BS[j]; //BkSk*(BkSk)T
		}
		sTy += s[i] * y[i];  //(Sk)T*Yk
	}
	for (i = 0; i<M; i++){
		for (j = 0; j<M; j++){
			B2[i][j] = B1[i][j] - (BSBST[i][j] / sTBs) + (yyT[i][j] / sTy); //BFGS main function
		}
	}
}

void flagx(){
	if (fabs(b[0]) <= min && fabs(b[1]) <= min){
		flag = 1;
	}
}

void print_result(){
	printf("loop=%d\n", loop);
	printf("x = %f,y = %f\n", x2[0], x2[1]);
	printf("d[0] = %f,d[1] = %f\n", d[0], d[1]);
	for (i = 0; i<M; i++){
		x1[i] = x2[i];
		b[i] = b1[i];
		for (j = 0; j<M; j++){
			if (mode == 1){
				printf("B2[%d][%d]=%f\n", i, j, B2[i][j]);
				B1[i][j] = B2[i][j];
			}
			else{
				printf("H1[%d][%d]=%f\n", i, j, H1[i][j]);
				H[i][j] = H1[i][j];
			}
		}
	}
	printf("\n\n");
}

int main(){
	initial();
	while (mode != 1 || mode !=2){
		printf("1.BFGS(Bk)\n2.BFGS(Hk)\n");
		scanf_s("%s", &mode);
		if (mode == 0)
			printf("file.7z.001\n");
	}
	while (flag == 0){
		if (mode == 1)
			LU();
		else
			hk();
		alpha_update();
		x2_update();
		flagx();
		if (mode == 1)
			BFGS();
		else
			BGFS_HK();
		loop += 1;
		print_result();
	}
}