#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define TRAPMAX 10
#define TRAPMAX_START 50
#define PI 3.14159265358979323846264338327950288419716939937510
#define NPN 6	//nucleosome point number
#define DIM 3
#define RANGE 0.1
#define EC 2.718281828459
#define TEMPPOW 2.0
#define TEMPPOW2 2.0
#define DISTCUTOFFSQ 4
#define FAC 1.5
//#define R_SPHERE_SQ 90000




void main()
{
char assignlistfilename[N],assignfilename[N],interlistfilename[N],interfilename[N],intralistfilename[N],intrafilename[N],outputassignfilename[N],outputtargetfilename[N];
int i,j,k,l,m,chrnum,*chrlen,internum,*interlen,*intralen,***intracontactlist,**interchridlist,***intercontactlist,***assignlist,**targetlist,pointnum,count,targetlen,tempcount;
FILE *assignlistfile,*assignfile,*interlistfile,*intrafile,*interfile,*intralistfile,*outputassignfile,*outputtargetfile;

printf("Input filename for assignment list:\n");
scanf("%s",assignlistfilename);

printf("Input filename for intra contact list:\n");
scanf("%s",intralistfilename);

printf("Input filename for inter contact list:\n");
scanf("%s",interlistfilename);

printf("Input filename for output assignment:\n");
scanf("%s",outputassignfilename);

printf("Input filename for output target contact:\n");
scanf("%s",outputtargetfilename);


chrnum = getlinenum(assignlistfilename);
internum = (chrnum-1)*(chrnum)/2;

chrlen = intarray(chrnum);
intralen = intarray(chrnum);
interlen = intarray(internum);
intracontactlist = intpointpointarray(chrnum);
interchridlist = intmatrix(internum,2);
intercontactlist = intpointpointarray(internum);
assignlist = intpointpointarray(chrnum);


//read file
assignlistfile = openfile(assignlistfilename,'r');
count = 0;
for(i=0;i<chrnum;i++)
 {
 fscanf(assignlistfile,"%s",assignfilename);
 chrlen[i] = getlinenum(assignfilename);
 assignlist[i] = intmatrix(chrlen[i],2);
 assignfile = openfile(assignfilename,'r');
 for(j=0;j<chrlen[i];j++)
  {
  fscanf(assignfile,"%d",&assignlist[i][j][0]);
  assignlist[i][j][1] = count;
  count ++;
  }
 fclose(assignfile);
 }
fclose(assignlistfile);

pointnum = count;

intralistfile = openfile(intralistfilename,'r');
targetlen = 0;
for(i=0;i<chrnum;i++)
 {
 fscanf(intralistfile,"%s",intrafilename);
 intralen[i] = getlinenum(intrafilename);
 targetlen += intralen[i];
 intracontactlist[i] = intmatrix(intralen[i],3);
 intrafile = openfile(intrafilename,'r');
 for(j=0;j<intralen[i];j++)
  {
  for(k=0;k<3;k++)
   fscanf(intrafile,"%d",&intracontactlist[i][j][k]);
  }
 fclose(intrafile);
 }
fclose(intralistfile);

interlistfile = openfile(interlistfilename,'r');
count = 0;
for(i=0;i<chrnum;i++)
 {
 for(j=i+1;j<chrnum;j++)
  {
//  count = 0;
  interchridlist[count][0] = i;
  interchridlist[count][1] = j;
  fscanf(interlistfile,"%s",interfilename);
  interlen[count] = getlinenum(interfilename);
  targetlen += interlen[count];
  intercontactlist[count] = intmatrix(interlen[count],3);
  interfile = openfile(interfilename,'r');
  for(k=0;k<interlen[count];k++)
   {
   for(l=0;l<3;l++)
    fscanf(interfile,"%d",&intercontactlist[count][k][l]);
   }
  fclose(interfile);

  count ++;
  }
 }

fclose(interlistfile);

targetlist = intmatrix(targetlen,3);

//get target list
count = 0;
for(i=0;i<chrnum;i++)
 {
 for(j=0;j<intralen[i];j++)
  {
  k = assignlist[i][intracontactlist[i][j][0]][1];
  l = assignlist[i][intracontactlist[i][j][1]][1];
//  printf("%d %d %d %d %d %d\n",targetlen,count,k,intracontactlist[i][j][0],l,intracontactlist[i][j][1]);
  targetlist[count][0] = k;
  targetlist[count][1] = l;
  targetlist[count][2] = intracontactlist[i][j][2];

  count ++;
  }
 }

tempcount = 0;
for(i=0;i<chrnum;i++)
 for(j=i+1;j<chrnum;j++)
  {
  for(m=0;m<interlen[tempcount];m++)
   {
//   printf("%d %d %d %d %d\n",i,j,tempcount,interlen[tempcount],m);
   k = assignlist[interchridlist[tempcount][0]][intercontactlist[tempcount][m][0]][1];
   l = assignlist[interchridlist[tempcount][1]][intercontactlist[tempcount][m][1]][1];
   targetlist[count][0] = k;
   targetlist[count][1] = l;
   targetlist[count][2] = intercontactlist[tempcount][m][2];
   count ++;
   }
  tempcount ++;
  }

outputassignfile = openfile(outputassignfilename,'w');
for(i=0;i<chrnum;i++)
 {
 for(j=0;j<chrlen[i];j++)
  fprintf(outputassignfile,"%d %d %d\n",assignlist[i][j][1],i,assignlist[i][j][0]);
 }
fclose(outputassignfile);


outputtargetfile = openfile(outputtargetfilename,'w');
for(i=0;i<targetlen;i++)
 {
 fprintf(outputtargetfile,"%d %d %d\n",targetlist[i][0],targetlist[i][1],targetlist[i][2]);
 }
fclose(outputtargetfile);





for(i=0;i<targetlen;i++)
 free(targetlist[i]);
free(targetlist);

for(i=0;i<chrnum;i++)
 {
 for(j=0;j<chrlen[i];j++)
  {
  free(assignlist[i][j]);
  }
 for(j=0;j<intralen[i];j++)
  free(intracontactlist[i][j]);
 free(assignlist[i]);
 free(intracontactlist[i]);
 }
free(assignlist);
free(intracontactlist);

for(i=0;i<internum;i++)
 {
 for(j=0;j<interlen[i];j++)
  free(intercontactlist[i][j]);
 free(intercontactlist[i]);
 free(interchridlist[i]);
 }
free(interchridlist);
free(intercontactlist);
free(interlen);
free(intralen);
free(chrlen);


}









