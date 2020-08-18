#include <stdio.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/mldaetlib/inlist.c"
#define N 200
#define DIM 3
#include <math.h>





double getdist(double *coor1,double *coor2,int dim)
 {
 int i;
 double distsq;
 distsq = 0;

 for(i=0;i<dim;i++)
  distsq += (coor1[i]-coor2[i])*(coor1[i]-coor2[i]);

 return(sqrt(distsq));
 }






int getdensity(double *distlist,int statenum,int currentid,double densitydist,int *assign)
 {
 int i,j,startid,endid;
 
 for(i=currentid;i>=0;i--)
  {
  if(i==0)
   {
   startid = i;
   break;
   }
  else if((distlist[i] > densitydist)||(assign[i] < 0))
   {
   startid = i+1;
   break;
   }
  }


 for(i=currentid;i<statenum;i++)
  {
  if(i==statenum-1)
   {
   endid = i;
   break;
   }
  else if((distlist[i] > densitydist)||(assign[i] < 0))
   {
   endid = i-1;
   break;
   }
  }

 return(endid-startid+1);

 }







 
int getneighborpeak(int *density,int statenum,int currentid,int genedist)
 {
 int i,maxvalue,maxid,startid,endid;
 
 startid = currentid-genedist;
 endid = currentid+genedist;

 if(startid < 0)
  startid = 0;
 if(endid > statenum-1)
  endid = statenum-1;

 maxvalue = 0;
 maxid = currentid;
 for(i=startid;i<=endid;i++)
  {
  if(density[i] >= maxvalue)
   {
   maxvalue = density[i];
   maxid = i;
   }
  }

 return(maxid);
 }







int getneighborhighpeak(double *density,int statenum,int currentid,int genedist,int *assign)
 {
 int i,minid,startid,endid;
 double minvalue;

 startid = currentid-genedist;
 endid = currentid+genedist;

 if(startid < 0)
  startid = 0;
 if(endid > statenum-1)
  endid = statenum-1;

 minvalue = density[currentid];
// minvalue = 1000;
 minid = currentid;
// printf("%d\n",currentid);
 for(i=startid;i<=endid;i++)
  {
  if((density[i] > minvalue)&&(assign[i] != -1))
   {
   minvalue = density[i];
   minid = i;
   }
  }

 return(minid);
 }








int getneighborlowpeak(double *density,int statenum,int currentid,int genedist,int *assign)
 {
 int i,minid,startid,endid;
 double minvalue;

 startid = currentid-genedist;
 endid = currentid+genedist;

 if(startid < 0)
  startid = 0;
 if(endid > statenum-1)
  endid = statenum-1;

 minvalue = density[currentid];
// minvalue = 1000;
 minid = currentid;
// printf("%d\n",currentid);
 for(i=startid;i<=endid;i++)
  {
  if((density[i] < minvalue)&&(assign[i] != -1))
   {
   minvalue = density[i];
   minid = i;
   }
  }

 return(minid);
 }








int getfinalassign(int *peakassign,int statenum,int currentid,int *assignbool)
 {
 int finalid;
 if(peakassign[currentid] == currentid)
  finalid = currentid;
 else
  finalid = getfinalassign(peakassign,statenum,peakassign[currentid],assignbool);

 assignbool[currentid] = 1;
 return(finalid);

 }






void main()
{
char inputfilename[N],outputfilename[N],assignfilename[N],peakassignfilename[N];
int i,j,statenum,*peakidlist,peaknum,genedist,*peakassign,*assign,*assignbool;
double *avelist,*sigmalist;
FILE *inputfile,*outputfile,*assignfile,*peakassignfile;

printf("Input filename for original data:\n");
scanf("%s",inputfilename);

//printf("Input filename for point assignment:\n");
//scanf("%s",assignfilename);

printf("Input the genome dist to consider peaks:\n");
scanf("%d",&genedist);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

printf("Input filename for peakassign:\n");
scanf("%s",peakassignfilename);

statenum = getlinenum(inputfilename);

avelist = doublearray(statenum);
sigmalist = doublearray(statenum);
peakidlist = intarray(statenum);
peakassign = intarray(statenum);
assign = intarray(statenum);
assignbool = intarray(statenum);


inputfile = openfile(inputfilename,'r');
for(i=0;i<statenum;i++)
 {
 fscanf(inputfile,"%lf",&avelist[i]);
 sigmalist[i] = 0;
// fscanf(inputfile,"%lf",&sigmalist[i]);
 }
fclose(inputfile);

//assignfile = openfile(assignfilename,'r');
for(i=0;i<statenum;i++)
 assign[i] = 1;
// fscanf(assignfile,"%d",&assign[i]);
//fclose(assignfile);

//get peak assign
for(i=0;i<statenum;i++)
 {
 if(assign[i] < 0)
  peakassign[i] = -1;
 else
  peakassign[i] = getneighborhighpeak(avelist,statenum,i,genedist,assign);
//  peakassign[i] = getneighborlowpeak(avelist,statenum,i,genedist,assign);
// printf("%d\n",peakassign[i]);
 }


//overall assign
for(i=0;i<statenum;i++)
 assignbool[i] = 0;

for(i=0;i<statenum;i++)
 {
// printf("%d\n",i);
 if(assign[i] < 0)
  assignbool[i] = 1;
 else
  peakassign[i] = getfinalassign(peakassign,statenum,i,assignbool);
// printf("%d\n",peakassign[i]);
 }
//printf("%d\n",peakassign[2920]);
//printf("%d\n",peakassign[2921]);

peaknum = 0;
for(i=0;i<statenum;i++)
 {
 if((peakassign[i] >= 0)&&(peaknum == 0))
  {
  peakidlist[peaknum] = peakassign[i];
  peaknum ++;
  }
 else if((peakassign[i] >=0)&&(peaknum != 0))
  {
  if(inintlist(peakidlist,peaknum,peakassign[i]) < 0)
   {
   peakidlist[peaknum] = peakassign[i];
//   printf("%d %d\n",peaknum,peakidlist[peaknum]);
   peaknum ++;
   }
  }
 }

//printf("%d\n",peaknum);

outputfile = openfile(outputfilename,'w');
for(i=0;i<peaknum;i++)
 fprintf(outputfile,"%d %lf %lf\n",peakidlist[i],avelist[peakidlist[i]],sigmalist[peakidlist[i]]);
fclose(outputfile);

peakassignfile = openfile(peakassignfilename,'w');
for(i=0;i<statenum;i++)
 fprintf(peakassignfile,"%d\n",peakassign[i]);
fclose(peakassignfile);




free(avelist);
free(sigmalist);
free(peakidlist);
free(peakassign);
free(assign);
free(assignbool);


}












