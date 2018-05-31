#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdbool.h>
#include "readConfigure.h"

void showBmpHead(BITMAPFILEHEADER* pBmpHead)
{
printf("位图文件头:\n");
printf("文件大小:%d\n",pBmpHead->bfSize);
printf("保留字:%d\n",pBmpHead->bfReserved1);
printf("保留字:%d\n",pBmpHead->bfReserved2);
printf("实际位图数据的偏移字节数:%d\n",pBmpHead->bfOffBits);

}


void showBmpInforHead(BITMAPINFOHEADER* pBmpInforHead)
{
printf("位图信息头:\n");
printf("结构体的长度:%d\n",pBmpInforHead->biSize);
printf("位图宽:%d\n",pBmpInforHead->biWidth);
printf("位图高:%d\n",pBmpInforHead->biHeight);
printf("biPlanes平面数:%d\n",pBmpInforHead->biPlanes);
printf("biBitCount采用颜色位数:%d\n",pBmpInforHead->biBitCount);
printf("压缩方式:%d\n",pBmpInforHead->biCompression);
printf("biSizeImage实际位图数据占用的字节数:%d\n",pBmpInforHead->biSizeImage);
printf("X方向分辨率:%d\n",pBmpInforHead->biXPelsPerMeter);
printf("Y方向分辨率:%d\n",pBmpInforHead->biYPelsPerMeter);
printf("使用的颜色数:%d\n",pBmpInforHead->biClrUsed);
printf("重要颜色数:%d\n",pBmpInforHead->biClrImportant);
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
pfile = fopen(strFile, "rb");//打开文件

if(pfile!=NULL)
{
   printf("file bkwood.bmp open success.\n");
   //读取位图文件头信息
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

   //读取位图信息头信息
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

if(bitInfoHead.biBitCount < 24)//有调色板
{
   //读取调色盘结信息
   long nPlantNum = pow(2,(bitInfoHead.biBitCount));    //   Mix color Plant Number;
   pRgb=(RGBQUAD *)malloc(nPlantNum*sizeof(RGBQUAD));
   memset(pRgb,0,nPlantNum*sizeof(RGBQUAD));
   int num = fread(pRgb,4,nPlantNum,pfile);

   printf("Color Plate Number: %d\n",nPlantNum);

   printf("颜色板信息:\n");
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
//分配内存空间把源图存入内存
int l_width   = WIDTHBYTES(width* bitInfoHead.biBitCount);//计算位图的实际宽度并确保它为32的倍数
BYTE    *pColorData=(BYTE *)malloc(height*l_width);
memset(pColorData,0,height*l_width);
long nData = height*l_width;

//把位图数据信息读到数组里
fread(pColorData,1,nData,pfile);



//将位图数据转化为RGB数据
RGBQUAD* dataOfBmp;
dataOfBmp = (RGBQUAD *)malloc(width*height*sizeof(RGBQUAD));//用于保存各像素对应的RGB数据
memset(dataOfBmp,0,width*height*sizeof(RGBQUAD));

if(bitInfoHead.biBitCount<24)//有调色板，即位图为非真彩色
{
   int k;
   int index = 0;
   if (bitInfoHead.biBitCount == 1)
   {
    for(int i=0;i<height;i++)
     for(int j=0;j<width;j++)
     {
      BYTE mixIndex= 0;
      k = i*l_width + j/8;//k:取得该像素颜色数据在实际数据数组中的序号
      //j:提取当前像素的颜色的具体值
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

      //将像素数据保存到数组中对应的位置
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
      k = i*l_width + j/4;//k:取得该像素颜色数据在实际数据数组中的序号
      //j:提取当前像素的颜色的具体值
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

      //将像素数据保存到数组中对应的位置
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
      {//低
       mixIndex = mixIndex<<4;
       mixIndex = mixIndex>>4;
      }
      else
      {//高
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
else//位图为24位真彩色
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

//printf("像素数据信息:\n");
//for (int i=0; i<width*height; i++)
//{
//   if (i%5==0)
//   {
//    printf("\n");
//   }
//   showRgbQuan(&dataOfBmp[i]);
//}

//bug log ：注意要把width+1用括号括起来
int** pic = (int**)malloc(sizeof(int*)*(width+1));
for (int i = 0; i <= width; i++) {
	pic[i] = (int*)malloc(sizeof(int)*height);
}
//第一行作为图片头，存储一些像图片长和高这种信息,pic[0][0]为width，pic[0][1]为height
//拷贝图片内容和信息
//转化成黑白图片
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
/*printf("输入任意键结束\n");
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

	//注意pic的第一行是空的，用来存尺寸了

	//直方图均衡化(光照) + 主方向(旋转)
	match(pic1, pic2, keyPoints1, keyPoints2);


	free(keyPoints1);
	free(keyPoints2);
	free(pic1);
	free(pic2);
	while (1);
}
