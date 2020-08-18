#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/mldaetlib/inlist.c"
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
char contactfilename[N],assignfilename[N],outputfilename[N];
int i,j,k,l,m,statenum,**assign,contactnum,**contact,**contactmatrix,chrnum,id1,id2;
FILE *contactfile,*assignfile,*outputfile;

printf("Input filename for original contact:\n");
scanf("%s",contactfilename);

printf("Input filename for assignment:\n");
scanf("%s",assignfilename);

printf("Input filename for output:\n");
scanf("%s",outputfilename);


statenum = getlinenum(assignfilename);
contactnum = getlinenum(contactfilename);

assign = intmatrix(statenum,3);
contact = intmatrix(contactnum,3);

chrnum = 0;
assignfile = openfile(assignfilename,'r');
for(i=0;i<statenum;i++)
 {
 for(j=0;j<3;j++)
  fscanf(assignfile,"%d",&assign[i][j]);
 if(chrnum < assign[i][1])
  chrnum = assign[i][1];
 }
fclose(assignfile);
chrnum ++;

contactfile = openfile(contactfilename,'r');
for(i=0;i<contactnum;i++)
 {
 for(j=0;j<3;j++)
  {
  fscanf(contactfile,"%d",&contact[i][j]);
  }
 }
fclose(contactfile);


contactmatrix = intmatrix(chrnum,chrnum);
for(i=0;i<chrnum;i++)
 for(j=0;j<chrnum;j++)
  contactmatrix[i][j] = 0;

for(i=0;i<contactnum;i++)
 {
// printf("%d %d %d %d %d %d\n",i,contact[i][0],assign[contact[i][0]][1],contact[i][1],assign[contact[i][1]][1],chrnum);
// printf("%d %d\n",assign[contact[i][0]][1],assign[contact[i][1]][1]);
// id1=assign[contact[i][0]][1];
// id2=assign[contactmatrix[i][1]][1];
// printf("%d %d\n",id1,id2);
// contactmatrix[id1][id2] = 1;
// printf("%d %d\n",assign[contact[i][0]][1],assign[contactmatrix[i][1]][1]);
// contactmatrix[assign[contact[i][0]][1]][assign[contactmatrix[i][1]][1]] = 1;
 contactmatrix[assign[contact[i][0]][1]][assign[contact[i][1]][1]] = 1;
// printf("%d %d %d %d %d %d\n",i,contact[i][0],assign[contact[i][0]][1],contact[i][1],assign[contact[i][1]][1],chrnum);
 }

outputfile = openfile(outputfilename,'w');
for(i=0;i<chrnum;i++)
 for(j=i+1;j<chrnum;j++)
  {
  if(contactmatrix[i][j] != 0)
   fprintf(outputfile,"%d %d %d\n",i,j,1);
  }
fclose(outputfile);

for(i=0;i<chrnum;i++)
 free(contactmatrix[i]);
free(contactmatrix);

for(i=0;i<statenum;i++)
 free(assign[i]);
free(assign);

for(i=0;i<contactnum;i++)
 free(contact[i]);
free(contact);

}
