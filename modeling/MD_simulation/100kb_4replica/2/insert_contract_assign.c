#include <stdio.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define DIM 3
#include <math.h>



void main()
{
char inputfilename[N],assign1filename[N],assign2filename[N],outputfilename[N];
int i,j,k,statenum,insertnum,chrnum,*assign1,*chrlen1,totallen1,*assign2,*chrlen2,ratio,totallen2,tempid1,tempid2;
double **coor,**coor2;
FILE *inputfile,*assign1file,*assign2file,*outputfile;

printf("Input the filename for the conformation:\n");
scanf("%s",inputfilename);

printf("Input filename for first assignment:\n");
scanf("%s",assign1filename);

printf("Input filename for second assignment:\n");
scanf("%s",assign2filename);

printf("Input the filename for output:\n");
scanf("%s",outputfilename);


totallen1 = getlinenum(assign1filename);
totallen2 = getlinenum(assign2filename);
statenum = getlinenum(inputfilename);
chrnum = 20;

insertnum = (int)(totallen2/((double)totallen1)+0.5);
printf("%d\n",insertnum);


if(totallen1 != statenum)
 {
 printf("ERROR for the input\n");
 exit(0);
 }

assign1 = intarray(totallen1);
assign2 = intarray(totallen2);
chrlen1 = intarray(chrnum);
chrlen2 = intarray(chrnum);

coor = doublematrix(statenum,DIM);
coor2 = doublematrix(totallen2,DIM);

inputfile = openfile(inputfilename,'r');
for(i=0;i<statenum;i++)
 for(j=0;j<DIM;j++)
  fscanf(inputfile,"%lf",&coor[i][j]);
fclose(inputfile);

for(i=0;i<chrnum;i++)
 {
 chrlen1[i] = 0;
 chrlen2[i] = 0;
 }

assign1file = openfile(assign1filename,'r');
for(i=0;i<totallen1;i++)
 {
 fscanf(assign1file,"%d",&j);
 fscanf(assign1file,"%d",&assign1[i]);
 fscanf(assign1file,"%d",&j);
 chrlen1[assign1[i]] ++;
 }
fclose(assign1file);

assign2file = openfile(assign2filename,'r');
for(i=0;i<totallen2;i++)
 {
 fscanf(assign2file,"%d",&j);
 fscanf(assign2file,"%d",&assign2[i]);
 fscanf(assign2file,"%d",&j);
 chrlen2[assign2[i]] ++;
 }
fclose(assign2file);



tempid1 = 0;
tempid2 = 0;
for(i=0;i<chrnum;i++)
 {
 for(j=0;j<chrlen2[i];j++)
  {
  for(k=0;k<DIM;k++)
   {
   if(j/insertnum < chrlen1[i]-1)
    coor2[tempid2+j][k] = (coor[tempid1+j/insertnum][k]*(insertnum+(j/insertnum)*insertnum-j)+coor[tempid1+j/insertnum+1][k]*(j-(j/insertnum)*insertnum))/insertnum;
   else
    coor2[tempid2+j][k] = (coor[tempid1+j/insertnum-1][k]*((j/insertnum)*insertnum-j)+coor[tempid1+j/insertnum][k]*(insertnum+j-(j/insertnum)*insertnum))/insertnum;

   }
  }
 tempid1 += chrlen1[i];
 tempid2 += chrlen2[i];
 }





outputfile = openfile(outputfilename,'w');
for(i=0;i<totallen2;i++)
 {
 for(j=0;j<DIM;j++)
  fprintf(outputfile,"%lf ",coor2[i][j]);
 fprintf(outputfile,"\n");
 }


fclose(outputfile);

for(i=0;i<statenum;i++)
 {
 free(coor[i]);
 }
free(coor);

for(i=0;i<totallen2;i++)
 {
 free(coor2[i]);
 }
free(coor2);
free(assign1);
free(assign2);
free(chrlen1);
free(chrlen2);





}







