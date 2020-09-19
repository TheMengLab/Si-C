#include <stdio.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200
#define DIM 3
#include <math.h>




int getnextid(int *assign,int statenum,int currentid)
 {
 int i,j,nextid,endbool;

 nextid = currentid;
 endbool = 0;
 for(i=currentid+1;(i<statenum)&&(endbool == 0);)
  {
  if(assign[i] == -1)
   {
   i ++;
   nextid = i;
   }
  else
   endbool = 1;
  }
 
 return(nextid);
 }







void modifycoor(double **coor,int statenum,int currentid,int nextid)
 {
 int i,j;
 double range;

 range = nextid-currentid;

 for(i=currentid+1;(i<nextid)&&(i<statenum);i++)
  {
  for(j=0;j<DIM;j++)
   coor[i][j] = ((nextid-i)*coor[currentid][j]+(i-currentid)*coor[nextid][j])/range;
  }

 }



//void modifycoor_chr(double **coor,int statenum,int currentid,int nextid,int *chrassign)



void getsmoothcoor(double **coor,int statenum,int currentid,int len,double *smoothcoor)
 {
 int i,j,k,pointnum;

 pointnum = 2*len+1;
 for(j=0;j<DIM;j++)
  smoothcoor[j] = 0;

 for(i=currentid-len;i<=currentid+len;i++)
  for(j=0;j<DIM;j++)
   smoothcoor[j] += coor[i][j];

 for(j=0;j<DIM;j++)
  smoothcoor[j] /= pointnum;

 }







void main()
{
char input1filename[N],input2filename[N],outputfilename[N],assignfilename[N],chrassignfilename[N];
int i,j,k,statenum,maxlen,*assign,templen,nextid,*chrassign;
double **coor,**tempcoor;
FILE *input1file,*input2file,*outputfile,*assignfile,*chrassignfile;

printf("Input the filename for the conformation:\n");
scanf("%s",input1filename);

printf("Input the maximum window for smooth:\n");
scanf("%d",&maxlen);

printf("Input the filename for assignment:\n");
scanf("%s",assignfilename);

printf("Input filename for chr assignfilename:\n");
scanf("%s",chrassignfilename);

printf("Input the filename for output:\n");
scanf("%s",outputfilename);

statenum = getlinenum(input1filename);
//insertnum = 20;

coor = doublematrix(statenum,DIM);
tempcoor = doublematrix(statenum,DIM);
assign = intarray(statenum);
chrassign = intarray(statenum);
//tempcoor = doublematrix(statenum-insertnum,DIM);


input1file = openfile(input1filename,'r');
for(i=0;i<statenum;i++)
 for(j=0;j<DIM;j++)
  {
  fscanf(input1file,"%lf",&coor[i][j]);
  }
fclose(input1file);


//read assignment
assignfile = openfile(assignfilename,'r');
for(i=0;i<statenum;i++)
 fscanf(assignfile,"%d",&assign[i]);
fclose(assignfile);

chrassignfile = openfile(chrassignfilename,'r');
for(i=0;i<statenum;i++)
 {
 fscanf(chrassignfile,"%d",&chrassign[i]);
 }
fclose(chrassignfile);




nextid = 0;
if(assign[0] == -1)
 {
 nextid = getnextid(assign,statenum,0);

for(i=0;i<nextid;i++)
 for(j=0;j<DIM;j++)
  coor[i][j] = coor[nextid][j];
 }
//modify coor
for(i=nextid;i<statenum;i++)
 {
 nextid = getnextid(assign,statenum,i);
 if(nextid >= statenum)
  nextid = statenum-1;
// printf("%d %d\n",i,nextid);
 if(nextid != i)
  {
  if(chrassign[i] == chrassign[nextid])
   modifycoor(coor,statenum,i,nextid);
  else
   {
   for(j=i+1;j<nextid;j++)
    {
    if(chrassign[j] == chrassign[i])
     for(k=0;k<DIM;k++)
      coor[j][k] = coor[i][k];
    else if(chrassign[j] == chrassign[nextid])
     for(k=0;k<DIM;k++)
      coor[j][k] = coor[nextid][k];
    else
     {
     printf("ERROR for the input:\n");
     exit(0);
     }
    }
   }
  }
//  modifycoor_chr(coor,statenum,i,nextid,chrassign);
 }

//smooth the coor
for(i=0;i<statenum;i++)
 {
 templen = maxlen;
 if(i < maxlen)
  templen = i;
 if(i > statenum-1-maxlen)
  templen = statenum-1-i;
 
 getsmoothcoor(coor,statenum,i,templen,tempcoor[i]);
 }


outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum;i++)
 {
 for(j=0;j<DIM;j++)
  fprintf(outputfile,"%lf ",tempcoor[i][j]/1);
 fprintf(outputfile,"\n");
 }
fclose(outputfile);



for(i=0;i<statenum;i++)
 {
 free(coor[i]);
 }
free(coor);

for(i=0;i<statenum;i++)
 free(tempcoor[i]);
free(tempcoor);

free(assign);
free(chrassign);






}







