#include <stdio.h>
#include <stdlib.h>
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/mldaetlib/allocate.c"
#define N 200
#define DIM 3

void main()
 {
 char inputfilename[N],outputfilename[N],assignfilename[N],res;
 int i,j,eof1,eof2,beginindex,index,temp,linenum,*assign,id;
 double coor[DIM];
 FILE *inputfile,*outputfile,*assignfile;

 printf("Input the filename of the cavity coordinates:\n");
 scanf("%s",inputfilename);

 printf("Input the beginning index for the cavities(The number of atoms+1):\n");
 scanf("%d",&beginindex);

 printf("Input the filename for assignment:\n");
 scanf("%s",assignfilename);

 printf("Input the filename for the output pdb filename:\n");
 scanf("%s",outputfilename);

 linenum = getlinenum(inputfilename);
 printf("%d\n",linenum);

 assign = intarray(linenum);

 assignfile = openfile(assignfilename,'r');
 for(i=0;i<linenum;i++)
  {
  fscanf(assignfile,"%d",&assign[i]);
  }
 fclose(assignfile);

 inputfile = openfile(inputfilename,'r');
 outputfile = openfile(outputfilename,'w');

 id = 1;

 res = 'a';
// printf("%c\n",res);
 for(i=0;i<linenum;i++)
  {
//  printf("%d\n",i);
  index = i+beginindex;
//  fscanf(inputfile,"%d",&temp);
//  printf("%d\n",temp);
  for(j=0;j<DIM;j++)
   {
   fscanf(inputfile,"%lf",&coor[j]);
//   printf("%lf ",coor[j]);
   }
//  printf("\n");
  if(i > 0)
   {
   if(assign[i] != assign[i-1])
    {
    id = 1;
    res ++;
    }
   else
    id ++;
   }
//   printf("%c\n",res);
   fprintf(outputfile,"ATOM %6d  C     %c%6d    %8.3lf%8.3lf%8.3lf\n",index%1000000,res,id,coor[0],coor[1],coor[2]);
//  fprintf(outputfile,"ATOM %6d  C   XXX     3    %8.3lf%8.3lf%8.3lf\n",index%1000000,coor[0]/10,coor[1]/10,coor[2]/10);
  }
 

 fclose(inputfile);
 fclose(outputfile);

 free(assign);
 }



