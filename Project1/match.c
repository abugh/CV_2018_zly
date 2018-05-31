#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdbool.h>
#include "readConfigure.h"

double calKPDis(double**dpt1,double** dpt2,int size,int num) {
	double res=0;
	for (int i = 0; i < num*num; i++) {
		for (int j = 0; j < 8; j++) {
			res += ((dpt1[i][j] - dpt2[i][j]) / (size*size)) * ((dpt1[i][j] - dpt2[i][j]) / (size*size));
		}
	}
	return sqrt(res);
}

int** match(int** pic1,int** pic2,int** keyPoints1,int** keyPoints2) {
	int kp_size1 = keyPoints1[0][0];
	int kp_size2 = keyPoints2[0][0];
	//max(kp_size1, kp_size2) * 4大小,4依次为x1,y1,x2,y2
	//res暂时没用到 =_=!!
	int ** res = (int**)malloc(sizeof(int*)*max(kp_size1, kp_size2));
	int size = 4; int num = 4;
	double*** dpt1 = descriptor(pic1, keyPoints1, size, num);
	double*** dpt2 = descriptor(pic2, keyPoints2, size, num);

	for (int i = 0; i < kp_size1; i++) {
		double min = 1000000000;
		int index = -1;
		for (int j = 0; j < kp_size2; j++) {
			double dis = calKPDis(dpt1[i], dpt2[j], size,num);
			if (dis < min) {
				min = dis;
				index = j;
			}
		}
		printf("图1中 %d ( %d , %d ) 的匹配点为:图2中 %d ( %d , %d ),MSE误差为 %lf\n",
			i, dpt1[i][0], dpt1[i][1], index, dpt2[index][0], dpt2[index][1], min);
	}

	free(dpt1);
	free(dpt2);
	return res;
}