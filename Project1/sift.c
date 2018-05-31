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

//debug log: sizeof double ����intһ����ҪŪ��,double 8 ���ֽڣ�int 4 ���ֽ�
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
			//debug log:�˴������ڿ�ͷ#include <math.h>������������exp��Ҫ��Ȼ��ʹ��һ����ֵĺ�������������������
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
				//debug log ! < ��>=���ȼۣ���֪��Ϊɶ...
				//������>�������Ļ�һ��0���ᱻѡ����
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
	int** kp = (int**)malloc(kp_size);//��������������
	int kp_index = 1;
	double*** filters = (double***)malloc(sizeof(double**) * numFilters);
	for (int sigma = 1; sigma <= numFilters; sigma += 1) {
		filters[sigma - 1] = GaussFilter(kernelSize, pow(2,(double)sigma)/(double)2);//������2�ı���(��׼���Ǹ���2�ı���)
	}

	/*printf("ԭͼƬ\n");
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

	//û���߶Ȳ���(��Ϊ��ҵ�������һ����С��ͼƬ)
	//���ɾ���ͼ

	//debug log:��˹�������ͼ��Ҫ��double���ͣ�������ֵ��ʮ�ֽӽ�

	//int*** pic_fil = (int***)malloc(sizeof(int**) * numFilters); 
	//int*** pic_fil_comp = (int***)malloc(sizeof(int**) * numFilters);
	double*** pic_fil = (int***)malloc(sizeof(int**) * numFilters); 
	double*** pic_fil_comp = (int***)malloc(sizeof(int**) * numFilters);//��ֲ�ͼ
	int center = kernelSize / 2;
	for (int fil = 0; fil<numFilters; fil++) {
		pic_fil[fil] = generateDynamic2DoubleArray(width, height);

		//��ԭͼƬ��ʼ��pic_fil
		//printf("ԭͼ����ֵ:\n");
		for (int i = 1; i <= width; i++) {
			for (int j = 0; j < height; j++) {
				pic_fil[fil][i - 1][j] = pic[i][j];
				//printf("%lf ", pic_fil[fil][i - 1][j]);
			}
		}
		//printf("�˲�������ֵ:\n");
		for (int i = kernelSize / 2; i < width - kernelSize / 2; i++) {
			for (int j = kernelSize / 2; j < height - kernelSize / 2; j++) {
				double sum = 0;
				for (int ii = -kernelSize / 2; ii <= kernelSize / 2; ii++) {
					for (int jj = -kernelSize / 2; jj <= kernelSize / 2; jj++) {
						//debug log:������һ��double���͵�sum���Ӻͣ���ֱ����pic_fil��intԪ�أ���ֹ���ȶ�ʧ(��ʧЧ)
						sum += pic_fil[fil][i + ii][j + jj] * filters[fil][center + ii][center + jj];
					}
				}
				pic_fil[fil][i][j] = max(0,min(sum,255));
				//printf("%lf  ",pic_fil[fil][i][j]);
			}
		}
	}
	//ȡ���
	//ֱ����pic_filԭ��������޸�
	//printf("���ͼ����ֵ:\n");
	for (int fil = 0; fil < numFilters-1; fil++) {//�󷽲��С����
		for (int i = kernelSize / 2; i < width - kernelSize / 2; i++) {
			for (int j = kernelSize / 2; j < height - kernelSize / 2; j++) {
				pic_fil[fil][i][j] = min(max(0, pic_fil[fil][i][j] - pic_fil[fil + 1][i][j]),255);
				//printf("%lf ", pic_fil[fil][i][j]);
			}
		}
	}

	//��ͼ�ļ�Сֵ�൱���󼫴�ֵ
	//���˹��ֵĲ�ͼ
	//printf("���ͼ��ͼ��ֵ:\n");
	for (int fil = 0; fil < numFilters-1; fil++) {
		pic_fil_comp[fil] = generateDynamic2DoubleArray(width, height);
		for (int i = kernelSize / 2; i < width - kernelSize / 2; i++) {
			for (int j = kernelSize / 2; j < height - kernelSize / 2; j++) {
				pic_fil_comp[fil][i][j] = min(255,max(255 - pic_fil[fil][i][j],0));
				//printf("%lf ", pic_fil_comp[fil][i][j]);
			}
		}
	}


	//��ֵ����������ֵ�ͼ�Сֵ
	//debug log : pic_fil��padding���ֻ���ԭ��������ֵ,���м䲿���Ѿ��ǲ�ֵ������ȽϵĻ��Ǵ����
	for (int i = kernelSize / 2+1; i < width - kernelSize / 2-1; i++) {
		for (int j = kernelSize / 2+1; j < height - kernelSize / 2-1; j++) {
			for (int fil = 1; fil < numFilters - 2; fil++) {//debug log:�ѶԾ������ѭ�����������棬��ΪҪ��֤ÿһ����ֻ��kp�б���һ��
				if (isLocalMin(fil, i, j, pic_fil)) {
					kp[kp_index] = (int*)malloc(sizeof(int) * 2);
					kp[kp_index][0] = i; kp[kp_index][1] = j;
					kp_index++;
					break;
				}
				else if (isLocalMin(fil, i, j, pic_fil_comp)) {//�ڲ�ͼ�����Ǽ�Сֵ
					kp[kp_index] = (int*)malloc(sizeof(int) * 2);
					kp[kp_index][0] = i; kp[kp_index][1] = j;
					kp_index++;
					break;
				}
			}
		}
	}

	//printf("����������: %d\n", kp_index-1);

	for (int i = 0; i < numFilters; i++) {
		free(filters[i]);
		free(pic_fil[i]);
		if (i < numFilters - 1) {//pic_fil_comp������������һ��
			free(pic_fil_comp[i]);
		}
	}

	/*printf("SIFT ������:\n");
	for (int i = 1; i <= 5; i++) {
		printf("%d  %d\n", kp[i][0], kp[i][1]);
	}*/

	//kp��1��ʼ����,0������size
	kp[0] = (int*)malloc(sizeof(int));
	kp[0][0] = kp_index - 1;

	return kp;
}