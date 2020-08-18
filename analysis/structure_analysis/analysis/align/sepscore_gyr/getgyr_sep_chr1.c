#include <stdio.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define DIM 3





double getdist(double *coor1,double *coor2,int dim)
 {
 int i;
 double distsq;
 distsq = 0;

 for(i=0;i<dim;i++)
  distsq += (coor1[i]-coor2[i])*(coor1[i]-coor2[i]);

 return(sqrt(distsq));
 }




int getgyrbool(int **assign,int pointnum,int id,int interval)
 {
 int i,j,k,gyrbool;
 gyrbool = 1;

 if(id+interval > pointnum-1)
  {
  gyrbool = -1;
  }
 else
  {
  for(i=0;(i<interval)&&(gyrbool == 1);i++)
   {
   if((assign[id][1]!=assign[id+i][1])||(assign[id+i][2] == -1))
    gyrbool = -1;
   }
  }

 return(gyrbool);
 }






double getgyrvalue(double **coor,int pointnum,int id,int interval,double *center)
 {
 int i,j,k;
 double gyr;

 gyr = 0;
 for(j=0;j<DIM;j++)
  center[j] = 0;

 for(j=0;j<DIM;j++)
  {
  for(i=0;i<interval;i++)
   center[j] += coor[i+id][j];
  center[j] /= interval;
  }

 for(i=0;i<interval;i++)
  {
  for(j=0;j<DIM;j++)
   {
   gyr += (coor[id+i][j]-center[j])*(coor[id+i][j]-center[j]);
   }
  }

 gyr = sqrt(gyr/interval);

 return(gyr);

 }







double getgyrvalue_v2(double ***coor,int **assign,int pointnum,int id,int interval,double *center,int filenum)
 {
 int i,j,k,startid,endid;
 double gyr,tempgyr;

 startid = id-interval;
 endid = id+interval;

 if((startid < 0)||(endid >= pointnum))
  {
  gyr = 0;
  }
 else if((assign[startid][1] != assign[id][1])||(assign[endid][1] != assign[id][1]))
  {
  gyr = 0;
  }
 else
  {
  gyr = 0;
  for(i=0;i<filenum;i++)
   {

   for(j=0;j<DIM;j++)
    {
    center[j] = 0;

    for(k=startid;k<=endid;k++)
     center[j] += coor[i][k][j];

    center[j] /= 2*interval+1;
    }

   tempgyr = 0;
   for(k=startid;k<=endid;k++)
    {
    for(j=0;j<DIM;j++)
     tempgyr += (coor[i][k][j]-center[j])*(coor[i][k][j]-center[j]);
    }
   tempgyr = sqrt(tempgyr/(2*interval+1));
   gyr += tempgyr;
   }
  }

 gyr /= filenum;
 return(gyr);
 }







double getsepscore(double ***coor,int **assign,int pointnum,int id,int interval,int filenum)
 {
 int i,j,k,l,startid,endid;
 double score,tempscore;

 startid = id-interval;
 endid = id+interval;

 if((startid < 0)||(endid >= pointnum))
  {
  score = 0;
  }
 else if((assign[startid][1] != assign[id][1])||(assign[endid][1] != assign[id][1]))
  {
  score = 0;
  }
 else
  {
  score = 0;
  tempscore = 0;
  for(i=0;i<filenum;i++)
   for(k=startid;k<id;k++)
    for(l=id+1;l<=endid;l++)
     {
     tempscore += getdist(coor[i][k],coor[i][l],DIM);
     }
   tempscore /= filenum*interval*interval;

   }

 score = tempscore;

 return(score);
 }





void main()
{
char coorlistfilename[N],coorfilename[N],assignfilename[N],outputfilename[N];
int i,j,k,l,m,filenum,pointnum,**assign,*count,interval,*chrlen,chrnum,intervalnum,tempid,maxlen,endbool,chrid,startid,endid,*targetid,*gyrbool,tempstartid,tempendid;
double ***coor,*center,gyrvalue,sepscore;
FILE *coorlistfile,*coorfile,*assignfile,*outputfile;

printf("Input filename for original coordinate list:\n");
scanf("%s",coorlistfilename);

printf("Input filename for assignment:\n");
scanf("%s",assignfilename);

//printf("Input the chr id for target assignment:\n");
//scanf("%d",&chrid);
//
//printf("Input the interval for distance calculation:\n");
//scanf("%d",&interval);
//
printf("Input filename for output:\n");
scanf("%s",outputfilename);


filenum = getlinenum(coorlistfilename);
pointnum = getlinenum(assignfilename);

coor = doublematrixarray(filenum,pointnum,DIM);
assign = intmatrix(pointnum,3);
gyrbool = intarray(pointnum);
center = doublearray(3);
//chrlen = intarray(20);
chrid = 0;

//for(i=0;i<20;i++)
// chrlen[i] = 0;


maxlen = 0;
startid = pointnum-1;
endid = 0;
//printf("%s\n",assignfilename);
assignfile = openfile(assignfilename,'r');
for(i=0;i<pointnum;i++)
 {
 for(j=0;j<3;j++)
  fscanf(assignfile,"%d",&assign[i][j]);
// chrlen[assign[i][1]] ++;
 if(assign[i][1] == chrid)
  {
  if(i < startid)
   startid = i;
  if(i > endid)
   endid = i;
  }
// if(maxlen < chrlen[assign[i][1]])
//  maxlen = chrlen[assign[i][1]];
 }
fclose(assignfile);


coorlistfile = openfile(coorlistfilename,'r');
for(i=0;i<filenum;i++)
 {
 fscanf(coorlistfile,"%s",coorfilename);
 coorfile = openfile(coorfilename,'r');
 for(j=0;j<pointnum;j++)
  for(k=0;k<DIM;k++)
   fscanf(coorfile,"%lf",&coor[i][j][k]);
 fclose(coorfile);
 }
fclose(coorlistfile);




outputfile = openfile(outputfilename,'w');

interval = 10;
for(i=startid;i<=endid;i++)
 {
 gyrvalue = getgyrvalue_v2(coor,assign,pointnum,i,interval,center,filenum);
 sepscore = getsepscore(coor,assign,pointnum,i,50,filenum);
 fprintf(outputfile,"%d %lf %lf\n",i,gyrvalue,sepscore);
 }

fclose(outputfile);




//printf("test\n");

for(i=0;i<pointnum;i++)
 free(assign[i]);
free(assign);
//printf("test\n");
for(i=0;i<filenum;i++)
 {
 for(j=0;j<pointnum;j++)
  free(coor[i][j]);
 free(coor[i]);
 }
free(coor);

free(gyrbool);
free(center);

}



