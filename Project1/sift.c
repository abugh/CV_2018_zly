#include <stdbool.h>
#include "readConfigure.h"
#include <math.h>

double max(double a, double b) {
	if (a > b) {
		return a;
	}
	else {
		return b;
	}
}

double min(double a, double b) {
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

//debug log: sizeof double 或者int一定不要弄错,double 8 个字节，int 4 个字节
double** generateDynamic2DoubleArray(int width, int height) {
	double** arr = (double**)malloc(sizeof(int*)*width);
	for (int i = 0; i < width; i++) {
		arr[i] = (double*)malloc(sizeof(double)*height);
	}
	return arr;
}

int** generateDynamic2IntArray(int width, int height) {
	int** arr = (int**)malloc(sizeof(int*)*width);
	for (int i = 0; i < width; i++) {
		arr[i] = (int*)malloc(sizeof(int)*height);
	}
	return arr;
}

double** GaussFilter(int kernelSize, double sigma) {
	int n1, n2; double hg; double con;
	double** kernel = generateDynamic2DoubleArray(kernelSize, kernelSize);
	int center = kernelSize / 2; double sum = 0;
	for (int i = 0; i < kernelSize; i++) {
		for (int j = 0; j < kernelSize; j++) {
			n1 = abs(i - center);
			n2 = abs(j - center);
			con = -(double)(n1*n1 + n2*n2) / (double)(2 * sigma*sigma);
			//debug log:此处必须在开头#include <math.h>才能用真正的exp，要不然会使用一个奇怪的函数，编译器还不报错
			hg = exp(con);
			kernel[i][j] = hg;
			sum += hg;
		}
	}
	for (int i = 0; i < kernelSize; i++) {
		for (int j = 0; j < kernelSize; j++) {
			kernel[i][j] = kernel[i][j] / sum;
		}
	}
	return kernel;
}

bool isLocalMin(int fil, int i, int j,double*** pic_fil){
	for (int l = fil - 1; l <= fil + 1; l++) {
		for (int x = i - 1; x <= i + 1; x++) {
			for (int y = j - 1; y <= j + 1; y++) {
				if (l == fil && x == i && y == j) {
					continue;
				}
				//debug log ! < 和>=不等价？不知道为啥...
				//不能用>，这样的话一堆0都会被选进来
				if (pic_fil[fil][i][j] >= pic_fil[l][x][y]) {
					return false;
				}
			}
		}
	}

	return true;
}

int** SIFT(int** pic)
{
	int kernelSize = 3;  int numFilters = 5;
	int width = pic[0][0];
	int height = pic[0][1];
	int kp_size = sizeof(int*)*width*height / 100;
	int** kp = (int**)malloc(kp_size);//用来储存特征点
	int kp_index = 1;
	double*** filters = (double***)malloc(sizeof(double**) * numFilters);
	for (int sigma = 1; sigma <= numFilters; sigma += 1) {
		filters[sigma - 1] = GaussFilter(kernelSize, pow(2,(double)sigma)/(double)2);//方差是2的倍数(标准差是根号2的倍数)
	}

	/*printf("原图片\n");
	for (int x = 1; x <= width; x++) {
		for (int y = 0; y < height; y++) {
			printf("%d ", pic[x][y]);
		}
		printf("\n");
	}
	printf("\n");*/

	printf("Gauss Filter\n");
	for (int l = 0; l < numFilters; l++) {
		for (int x = 0; x < kernelSize; x++) {
			for (int y = 0; y < kernelSize; y++) {
				printf("%lf ", filters[l][x][y]);
			}
			printf("\n");
		}
		printf("\n");
	}

	//没做尺度不变(因为作业里就两张一样大小的图片)
	//生成卷积图

	//debug log:高斯卷积后的图像要用double类型，否则数值都十分接近

	//int*** pic_fil = (int***)malloc(sizeof(int**) * numFilters); 
	//int*** pic_fil_comp = (int***)malloc(sizeof(int**) * numFilters);
	double*** pic_fil = (int***)malloc(sizeof(int**) * numFilters); 
	double*** pic_fil_comp = (int***)malloc(sizeof(int**) * numFilters);//差分补图
	int center = kernelSize / 2;
	for (int fil = 0; fil<numFilters; fil++) {
		pic_fil[fil] = generateDynamic2DoubleArray(width, height);

		//用原图片初始化pic_fil
		//printf("原图像素值:\n");
		for (int i = 1; i <= width; i++) {
			for (int j = 0; j < height; j++) {
				pic_fil[fil][i - 1][j] = pic[i][j];
				//printf("%lf ", pic_fil[fil][i - 1][j]);
			}
		}
		//printf("滤波后像素值:\n");
		for (int i = kernelSize / 2; i < width - kernelSize / 2; i++) {
			for (int j = kernelSize / 2; j < height - kernelSize / 2; j++) {
				double sum = 0;
				for (int ii = -kernelSize / 2; ii <= kernelSize / 2; ii++) {
					for (int jj = -kernelSize / 2; jj <= kernelSize / 2; jj++) {
						//debug log:这里用一个double类型的sum来加和，不直接用pic_fil的int元素，防止精度丢失(已失效)
						sum += pic_fil[fil][i + ii][j + jj] * filters[fil][center + ii][center + jj];
					}
				}
				pic_fil[fil][i][j] = max(0,min(sum,255));
				//printf("%lf  ",pic_fil[fil][i][j]);
			}
		}
	}
	//取差分
	//直接在pic_fil原数组进行修改
	//printf("差分图像素值:\n");
	for (int fil = 0; fil < numFilters-1; fil++) {//大方差减小方差
		for (int i = kernelSize / 2; i < width - kernelSize / 2; i++) {
			for (int j = kernelSize / 2; j < height - kernelSize / 2; j++) {
				pic_fil[fil][i][j] = min(max(0, pic_fil[fil][i][j] - pic_fil[fil + 1][i][j]),255);
				//printf("%lf ", pic_fil[fil][i][j]);
			}
		}
	}

	//求补图的极小值相当于求极大值
	//求高斯差分的补图
	//printf("差分图像补图素值:\n");
	for (int fil = 0; fil < numFilters-1; fil++) {
		pic_fil_comp[fil] = generateDynamic2DoubleArray(width, height);
		for (int i = kernelSize / 2; i < width - kernelSize / 2; i++) {
			for (int j = kernelSize / 2; j < height - kernelSize / 2; j++) {
				pic_fil_comp[fil][i][j] = min(255,max(255 - pic_fil[fil][i][j],0));
				//printf("%lf ", pic_fil_comp[fil][i][j]);
			}
		}
	}


	//求极值，包括极大值和极小值
	//debug log : pic_fil中padding部分还是原来的像素值,而中间部分已经是差值，如果比较的话是错误的
	for (int i = kernelSize / 2+1; i < width - kernelSize / 2-1; i++) {
		for (int j = kernelSize / 2+1; j < height - kernelSize / 2-1; j++) {
			for (int fil = 1; fil < numFilters - 2; fil++) {//debug log:把对卷积层的循环放在最里面，因为要保证每一个点只在kp中保存一次
				if (isLocalMin(fil, i, j, pic_fil)) {
					kp[kp_index] = (int*)malloc(sizeof(int) * 2);
					kp[kp_index][0] = i; kp[kp_index][1] = j;
					kp_index++;
					break;
				}
				else if (isLocalMin(fil, i, j, pic_fil_comp)) {//在补图当中是极小值
					kp[kp_index] = (int*)malloc(sizeof(int) * 2);
					kp[kp_index][0] = i; kp[kp_index][1] = j;
					kp_index++;
					break;
				}
			}
		}
	}

	//printf("特征点数量: %d\n", kp_index-1);

	for (int i = 0; i < numFilters; i++) {
		free(filters[i]);
		free(pic_fil[i]);
		if (i < numFilters - 1) {//pic_fil_comp层数比其它少一层
			free(pic_fil_comp[i]);
		}
	}

	/*printf("SIFT 特征点:\n");
	for (int i = 1; i <= 5; i++) {
		printf("%d  %d\n", kp[i][0], kp[i][1]);
	}*/

	//kp从1开始计数,0用来存size
	kp[0] = (int*)malloc(sizeof(int));
	kp[0][0] = kp_index - 1;

	return kp;
}