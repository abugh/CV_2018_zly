#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdbool.h>
#include "readConfigure.h"

double* calGra(int i, int j, int** copy,int edge) {
	double * res = (double*)malloc(sizeof(double) * 2);
	for (int ii = max(0, i - 1); ii <= min(edge - 1, i + 1); ii++) {
		for (int jj = max(0, j - 1); jj <= min(edge - 1, j + 1); jj++) {
			if ((ii == i)&&(jj == j)) {
				continue;
			}
			double length = sqrt((ii - i)*(ii - i) + (jj - j)*(jj - j));
			//分别求在x,y两个方向上梯度
			double x = (double)(ii - i)*(double)(copy[ii][jj] - copy[i][j]) / length;
			double y = (double)(jj - j)*(double)(copy[ii][jj] - copy[i][j]) / length;
			
			res[0] = x/8; res[1] = y/8;
		}
	}

	return res;
}

double calDis(double a, double b) {
	return sqrt(a*a + b*b);
}

double myMod(double a, double b) {
	while (a < 0) {
		a += b;
	}
	return fmod(a, b);
}

//则整个描述符所占的邻域大小为(size * num)^2
//暂定size=4,num=4
double*** descriptor(int** pic, int** kp, int size, int num) {
	int kp_size = kp[0][0];
	printf("特征点数量 :  %d\n", kp_size);
	for (int i = 1; i <= kp_size; i++) {
		printf("特征点 :  %d  %d\n", kp[i][0], kp[i][1]);
	}

	double*** dpt = (double***)malloc(sizeof(double**)*(kp_size));
	for (int di = 0; di < kp_size; di++) {
		dpt[di] = (double**)malloc(sizeof(double*)*(num*num));
		for (int i = 0; i < num*num; i++) {
			dpt[di][i] = (double*)malloc(sizeof(double) * 8);
		}
	}

	//在整个邻域内((size * num)^2)使用Gauss，意思是，描述子中，离特征点越远的点影响越小
	double** GF = GaussFilter(size*num, sqrt(2));

	//直方图均衡化,减弱光照影响
	for (int ki = 1; ki <= kp_size; ki++) {
		int x = kp[ki][0]; int y = kp[ki][1];
		static const int a=1;
		int** copy = generateDynamic2IntArray(size*num, size*num);
		//把该特征点邻域内拷贝过来
		//debug log : 不能直接在原图上做，因为不同特征点邻域可能有重叠
		for (int i = 0; i < size*num; i++) {
			for (int j = 0; j < size*num; j++) {
				//pic第一行存着尺寸属性值
				if (x - size*num / 2 + 1 + i < 0 || x - size*num / 2 + 1 + i>=pic[0][0]) {
					printf("数组越界x %d  %d\n", x, x - size*num / 2 + 1 + i);
				}
				if (y - size*num / 2 + 1 + j < 0 || y - size*num / 2 + 1 + j>=pic[0][1]) {
					printf("数组越界y %d  %d\n", y, y - size*num / 2 + 1 + j);
				}
				//debug log : 横坐标值加2(比纵坐标多1)是因为横坐标从1开始
				copy[i][j] = pic[x - size*num / 2 + 2 + i][y - size*num / 2 + 1 + j];
			}
		}

		//对每个window内直方图均衡化 , 因为每个window就是一个描述子
		for (int i = 0; i < size*num; i+=size) {
			for (int j = 0; j < size*num; j+=size) {
				
				//求极大极小值
				int min = 255; int max = 0;
				for (int ii = i; ii < i + size; ii++) {
					for (int jj = j; jj < j + size; jj++) {
						if (copy[ii][jj] > max) {
							max = copy[ii][jj];
						}
						if (copy[ii][jj] < min) {
							min = copy[ii][jj];
						}
					}
				}

				//均衡化
				//printf("均衡化后像素值:\n");
				int scale = max - min;
				for (int ii = i; ii < i + size; ii++) {
					for (int jj = j; jj < j + size; jj++) {
						copy[ii][jj] = (double)(copy[ii][jj] - min) * (double)255 / (double)scale;
						//printf("%d  ", copy[ii][jj]);
					}
					//printf("\n");
				}

			}
		}

		//对邻域内每一个点计算梯度
		double sum_x = 0, sum_y = 0;
		//存储梯度长度和偏转角度
		double*** gradient = (double***)malloc(sizeof(double**)*(size*num));
		for (int i = 0; i < num*size; i++) {
			gradient[i] = (double**)malloc(sizeof(double*)*(size*num));
			for (int j = 0; j < num*size; j++) {
				gradient[i][j] = (double*)malloc(sizeof(double)*(2));
				double*gra = calGra(i,j,copy, size*num);
				gradient[i][j][0] = calDis(gra[0], gra[1]);
				gradient[i][j][1] = atan2(gra[1], gra[0]);//atan2函数y值在前
				free(gra);
				sum_x += gradient[i][j][0];
				sum_y += gradient[i][j][1];
			}
		}

		//主方向矫正
		//计算主方向偏转角度
		double theta = atan2(sum_y, sum_x);
		for (int i = 0; i < num*size; i++) {
			for (int j = 0; j < num*size; j++) {
				gradient[i][j][1] -= theta;
			}
		}

		//梯度高斯加权和
		for (int i = 0; i < num*size; i++) {
			for (int j = 0; j < num*size; j++) {
				gradient[i][j][0] = gradient[i][j][0] * GF[i][j] * (double)(num*size) * (double)(num*size);
				gradient[i][j][1] = gradient[i][j][1] * GF[i][j] * (double)(num*size) * (double)(num*size);
			}
		}

		//构造描述子
		//初始化描述子
		for (int i = 0; i < num*num; i++) {
			for (int j = 0; j < 8; j++) {
				//ki从1开始
				dpt[ki - 1][i][j] = 0;
			}
		}
		int dpt_index = 0;
		for (int i = 0; i < size*num; i += size) {
			for (int j = 0; j < size*num; j += size) {

				for (int ii = i; ii < i + size; ii++) {
					for (int jj = j; jj < j + size; jj++) {
						//将360度划分为8个区域，每个区域45度
						int direction = (myMod(gradient[ii][jj][0], (double)360)) / (double)45;
						dpt[ki - 1][dpt_index][direction] += gradient[ii][jj][0];
					}
				}
				

				dpt_index++;
			}
		}

		free(gradient);
		free(copy);
	}

	free(GF);

	return dpt;
}