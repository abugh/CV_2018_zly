#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdbool.h>
#include "readConfigure.h"

void showBmpHead(BITMAPFILEHEADER* pBmpHead)
{
printf("λͼ�ļ�ͷ:\n");
printf("�ļ���С:%d\n",pBmpHead->bfSize);
printf("������:%d\n",pBmpHead->bfReserved1);
printf("������:%d\n",pBmpHead->bfReserved2);
printf("ʵ��λͼ���ݵ�ƫ���ֽ���:%d\n",pBmpHead->bfOffBits);

}


void showBmpInforHead(BITMAPINFOHEADER* pBmpInforHead)
{
printf("λͼ��Ϣͷ:\n");
printf("�ṹ��ĳ���:%d\n",pBmpInforHead->biSize);
printf("λͼ��:%d\n",pBmpInforHead->biWidth);
printf("λͼ��:%d\n",pBmpInforHead->biHeight);
printf("biPlanesƽ����:%d\n",pBmpInforHead->biPlanes);
printf("biBitCount������ɫλ��:%d\n",pBmpInforHead->biBitCount);
printf("ѹ����ʽ:%d\n",pBmpInforHead->biCompression);
printf("biSizeImageʵ��λͼ����ռ�õ��ֽ���:%d\n",pBmpInforHead->biSizeImage);
printf("X����ֱ���:%d\n",pBmpInforHead->biXPelsPerMeter);
printf("Y����ֱ���:%d\n",pBmpInforHead->biYPelsPerMeter);
printf("ʹ�õ���ɫ��:%d\n",pBmpInforHead->biClrUsed);
printf("��Ҫ��ɫ��:%d\n",pBmpInforHead->biClrImportant);
}

void showRgbQuan(RGBQUAD* pRGB)
{
printf("(%-3d,%-3d,%-3d)   ",pRGB->rgbRed,pRGB->rgbGreen,pRGB->rgbBlue);

}



int** readBMP(char* picName)
{

BITMAPFILEHEADER   bitHead;
BITMAPINFOHEADER bitInfoHead;
FILE* pfile;

char strFile[50];
strcpy(strFile, picName);
pfile = fopen(strFile, "rb");//���ļ�

if(pfile!=NULL)
{
   printf("file bkwood.bmp open success.\n");
   //��ȡλͼ�ļ�ͷ��Ϣ
   WORD fileType;
   fread(&fileType,1,sizeof(WORD),pfile);
   if(fileType != 0x4d42)
   {
    printf("file is not .bmp file!");
    return NULL;
   }
   //fseek(pfile,2,SEEK_CUR);   // "BM"
   fread(&bitHead,1,sizeof(BITMAPFILEHEADER),pfile);

   showBmpHead(&bitHead);
   printf("\n\n");

   //��ȡλͼ��Ϣͷ��Ϣ
   fread(&bitInfoHead,1,sizeof(BITMAPINFOHEADER),pfile);
   showBmpInforHead(&bitInfoHead);
   printf("\n");
}
else
{
   printf("file open fail!\n");
   return NULL;
}


RGBQUAD *pRgb ;

if(bitInfoHead.biBitCount < 24)//�е�ɫ��
{
   //��ȡ��ɫ�̽���Ϣ
   long nPlantNum = pow(2,(bitInfoHead.biBitCount));    //   Mix color Plant Number;
   pRgb=(RGBQUAD *)malloc(nPlantNum*sizeof(RGBQUAD));
   memset(pRgb,0,nPlantNum*sizeof(RGBQUAD));
   int num = fread(pRgb,4,nPlantNum,pfile);

   printf("Color Plate Number: %d\n",nPlantNum);

   printf("��ɫ����Ϣ:\n");
   for (int i =0; i<nPlantNum;i++)
   {
    if (i%5==0)
    {
     printf("\n");
    }
    showRgbQuan(&pRgb[i]);

   }

   printf("\n");

}


int width = bitInfoHead.biWidth;
int height = bitInfoHead.biHeight;
//�����ڴ�ռ��Դͼ�����ڴ�
int l_width   = WIDTHBYTES(width* bitInfoHead.biBitCount);//����λͼ��ʵ�ʿ�Ȳ�ȷ����Ϊ32�ı���
BYTE    *pColorData=(BYTE *)malloc(height*l_width);
memset(pColorData,0,height*l_width);
long nData = height*l_width;

//��λͼ������Ϣ����������
fread(pColorData,1,nData,pfile);



//��λͼ����ת��ΪRGB����
RGBQUAD* dataOfBmp;
dataOfBmp = (RGBQUAD *)malloc(width*height*sizeof(RGBQUAD));//���ڱ�������ض�Ӧ��RGB����
memset(dataOfBmp,0,width*height*sizeof(RGBQUAD));

if(bitInfoHead.biBitCount<24)//�е�ɫ�壬��λͼΪ�����ɫ
{
   int k;
   int index = 0;
   if (bitInfoHead.biBitCount == 1)
   {
    for(int i=0;i<height;i++)
     for(int j=0;j<width;j++)
     {
      BYTE mixIndex= 0;
      k = i*l_width + j/8;//k:ȡ�ø�������ɫ������ʵ�����������е����
      //j:��ȡ��ǰ���ص���ɫ�ľ���ֵ
      mixIndex = pColorData[k];
      switch(j%8)
      {
      case 0:
       mixIndex = mixIndex<<7;
       mixIndex = mixIndex>>7;
       break;
      case 1:
       mixIndex = mixIndex<<6;
       mixIndex = mixIndex>>7;
       break;
      case 2:
       mixIndex = mixIndex<<5;
       mixIndex = mixIndex>>7;
       break;

      case 3:
       mixIndex = mixIndex<<4;
       mixIndex = mixIndex>>7;
       break;
      case 4:
       mixIndex = mixIndex<<3;
       mixIndex = mixIndex>>7;
       break;

      case 5:
       mixIndex = mixIndex<<2;
       mixIndex = mixIndex>>7;
       break;
      case 6:
       mixIndex = mixIndex<<1;
       mixIndex = mixIndex>>7;
       break;

      case 7:
       mixIndex = mixIndex>>7;
       break;
      }

      //���������ݱ��浽�����ж�Ӧ��λ��
      dataOfBmp[index].rgbRed = pRgb[mixIndex].rgbRed;
      dataOfBmp[index].rgbGreen = pRgb[mixIndex].rgbGreen;
      dataOfBmp[index].rgbBlue = pRgb[mixIndex].rgbBlue;
      dataOfBmp[index].rgbReserved = pRgb[mixIndex].rgbReserved;
      index++;

    }
   }

   if(bitInfoHead.biBitCount==2)
   {
    for(int i=0;i<height;i++)
     for(int j=0;j<width;j++)
     {
      BYTE mixIndex= 0;
      k = i*l_width + j/4;//k:ȡ�ø�������ɫ������ʵ�����������е����
      //j:��ȡ��ǰ���ص���ɫ�ľ���ֵ
      mixIndex = pColorData[k];
      switch(j%4)
      {
      case 0:
       mixIndex = mixIndex<<6;
       mixIndex = mixIndex>>6;
       break;
      case 1:
       mixIndex = mixIndex<<4;
       mixIndex = mixIndex>>6;
       break;
      case 2:
       mixIndex = mixIndex<<2;
       mixIndex = mixIndex>>6;
       break;
      case 3:
       mixIndex = mixIndex>>6;
       break;
      }

      //���������ݱ��浽�����ж�Ӧ��λ��
      dataOfBmp[index].rgbRed = pRgb[mixIndex].rgbRed;
      dataOfBmp[index].rgbGreen = pRgb[mixIndex].rgbGreen;
      dataOfBmp[index].rgbBlue = pRgb[mixIndex].rgbBlue;
      dataOfBmp[index].rgbReserved = pRgb[mixIndex].rgbReserved;
      index++;


     }
   }
   if(bitInfoHead.biBitCount == 4)
   {
    for(int i=0;i<height;i++)
     for(int j=0;j<width;j++)
     {
      BYTE mixIndex= 0;
      k = i*l_width + j/2;
      mixIndex = pColorData[k];
      if(j%2==0)
      {//��
       mixIndex = mixIndex<<4;
       mixIndex = mixIndex>>4;
      }
      else
      {//��
       mixIndex = mixIndex>>4;
      }

      dataOfBmp[index].rgbRed = pRgb[mixIndex].rgbRed;
      dataOfBmp[index].rgbGreen = pRgb[mixIndex].rgbGreen;
      dataOfBmp[index].rgbBlue = pRgb[mixIndex].rgbBlue;
      dataOfBmp[index].rgbReserved = pRgb[mixIndex].rgbReserved;
      index++;

     }

   }
   if(bitInfoHead.biBitCount == 8)
   {
    for(int i=0;i<height;i++)
     for(int j=0;j<width;j++)
     {
      BYTE mixIndex= 0;

      k = i*l_width + j;

      mixIndex = pColorData[k];

      dataOfBmp[index].rgbRed = pRgb[mixIndex].rgbRed;
      dataOfBmp[index].rgbGreen = pRgb[mixIndex].rgbGreen;
      dataOfBmp[index].rgbBlue = pRgb[mixIndex].rgbBlue;
      dataOfBmp[index].rgbReserved = pRgb[mixIndex].rgbReserved;
      index++;



     }
   }
   if(bitInfoHead.biBitCount == 16)
   {
    for(int i=0;i<height;i++)
     for(int j=0;j<width;j++)
     {
      WORD mixIndex= 0;

      k = i*l_width + j*2;
      WORD shortTemp;
      shortTemp = pColorData[k+1];
      shortTemp = shortTemp<<8;

      mixIndex = pColorData[k] + shortTemp;

      dataOfBmp[index].rgbRed = pRgb[mixIndex].rgbRed;
      dataOfBmp[index].rgbGreen = pRgb[mixIndex].rgbGreen;
      dataOfBmp[index].rgbBlue = pRgb[mixIndex].rgbBlue;
      dataOfBmp[index].rgbReserved = pRgb[mixIndex].rgbReserved;
      index++;
     }
   }
}
else//λͼΪ24λ���ɫ
{
   int k;
   int index = 0;
   for(int i=0;i<height;i++)
    for(int j=0;j<width;j++)
    {
     k = i*l_width + j*3;
     dataOfBmp[index].rgbRed = pColorData[k+2];
     dataOfBmp[index].rgbGreen = pColorData[k+1];
     dataOfBmp[index].rgbBlue = pColorData[k];
     index++;
    }
}

//printf("����������Ϣ:\n");
//for (int i=0; i<width*height; i++)
//{
//   if (i%5==0)
//   {
//    printf("\n");
//   }
//   showRgbQuan(&dataOfBmp[i]);
//}

//bug log ��ע��Ҫ��width+1������������
int** pic = (int**)malloc(sizeof(int*)*(width+1));
for (int i = 0; i <= width; i++) {
	pic[i] = (int*)malloc(sizeof(int)*height);
}
//��һ����ΪͼƬͷ���洢һЩ��ͼƬ���͸�������Ϣ,pic[0][0]Ϊwidth��pic[0][1]Ϊheight
//����ͼƬ���ݺ���Ϣ
//ת���ɺڰ�ͼƬ
pic[0][0] = width; pic[0][1] = height;
int index = 0;
for (int i = 1; i <= width; i++) {
	for (int j = 0; j < height; j++) {
		pic[i][j] = (int)(0.299*dataOfBmp[index].rgbRed+0.587*dataOfBmp[index].rgbGreen+0.114*dataOfBmp[index].rgbBlue);
		index++;
	}
}

fclose(pfile);
if (bitInfoHead.biBitCount<24)
{
   free(pRgb);
}
free(dataOfBmp);
free(pColorData);
printf("\n");
/*printf("�������������\n");
char end[10];
scanf("%s", end);*/


return pic;
}




void main() 
{
	int** pic1 = readBMP("1.bmp");
	int** pic2 = readBMP("2.bmp");

	int** keyPoints1 = SIFT(pic1);
	int** keyPoints2 = SIFT(pic2);

	//ע��pic�ĵ�һ���ǿյģ�������ߴ���

	//ֱ��ͼ���⻯(����) + ������(��ת)
	match(pic1, pic2, keyPoints1, keyPoints2);


	free(keyPoints1);
	free(keyPoints2);
	free(pic1);
	free(pic2);
	while (1);
}
