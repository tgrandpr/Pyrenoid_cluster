/* This code was used to generate the unique molecular bonds and total 
bonds on Rubisco for simulations of an algal pyrenoid model involving 
EPYC1 and Rubisco in the paper T. GrandPre et al., Impact of Linker Length on Biomolecular Condensate Formation, PRX Life.
To run the code, a lammps file is needed to be read in. 
*/



#include <stdio.h>
#include <time.h>
#include <string.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <errno.h>
using namespace std;

/*define variables*/
int qt_types=6;
int ttdense;
int ttdenseR;
int mol;
int N= 18000;
double Lx=250.0;
double Ly=63;
double Lz=63;
int  frame;
   vector<int> cluster_number(N+1); 
    vector<int> bondr(N+1);
    vector<int> stickerbond(N+1);
void clusterHist();
int clustercheck(int, int, int);
 vector<int> check(N+1); 
 vector<int> nn_type(N+1); 
    double ave_size;
    double ave_total;
 int n_clusters=0;
  vector< vector<int> > ClusterSize(qt_types, vector<int>(N+1));
  vector<vector<double> > vvlist(N+1, vector<double>(N+1));//molecular neighborlist
    int cluster_id;
      double large_cluster;
    int large_cluster_id;


int main(int argc, char *argv[]){

frame=0;
// These are two files that are printed.
  FILE* Fout3;
      Fout3 = fopen("clustersize.txt", "w");

      FILE *  cluster4;
    cluster4 = fopen("clusterdist.txt", "w");


// here is the number of time slices in the simulation. 
int n_steps=501;

// these are the x, y, and z positions. 
vector<vector<double> > xx(n_steps+1, vector<double>(N+1));
vector<vector<double> > yy(n_steps+1, vector<double>(N+1));
vector<vector<double> > zz(n_steps+1, vector<double>(N+1));

// this is the molecular id used to calculate unique molecular bonds. 
vector<vector<int> > molid(N+1, vector<int>(3));

// This is a molecular bond indicator function. 
vector<vector<int> > bind(N+1, vector<int>(N+1));


//These are other variables to be used. 
vector<int> idd(N);
vector<int> bonds(N);
vector<int> bondsm(N);
vector<int> tbonds(N);
vector<vector<int> >  indexr(N+1, vector<int>(N+1));

int bound=0;
int bound1=0;
int tbound=0;
double rcut=1; // the cutoff of the interaction between type 1 (EPYC1) and type 2 (binding sites of RUBISCO) 

//These are other files that will be printed out such as complexes in the dense and dilute phases. 

 FILE * output2;
 output2=fopen("prob2.txt","w");


   FILE * x_out;
    x_out = fopen("dimerg.txt","w");

     FILE * x_out2;
    x_out2 = fopen("rubg.txt","w");

      FILE * x_out3;
    x_out3= fopen("out.lammpstrj","w");

      FILE * x_out4;
    x_out4 = fopen("dimergdense.txt","w");

     FILE * x_out5;
    x_out5 = fopen("rubgdense.txt","w");


/*Read in file*/  

  FILE *data;
  int t,j,i,b;
  char line;
  //THis is the simulation file that is read in to compute the complexes. 
  data=fopen("shortLinker_A7_ns0.4_Anneal_Rep1.lammpstrj","r");




  for (t=0; t<n_steps; t++) {

    // b is the number of word segments at the top of the lammps file/
      for (b=0; b<27; b++) {
      fscanf(data, "%s", &line);
      }
    
      for(int i=1;i<=N;i++){
      
       int k,type;
       double x,y,z;
      
       fscanf(data,"%d %d %lf %lf %lf",&k,&type,&x,&y,&z);
//These are defined on - LX to Lx and so forth so I change the domaim to be positive. 
       xx[t][k]=x+Lx;
       yy[t][k]=y+Ly;
       zz[t][k]=z+Lz;
       if(t==0){
       idd[k]=type;
       }

       
      
    }
    
    
  }
   fclose(data);

//Assign the molecule IDs
double rmol=0;
double emol=0;
   int mcnt=0;//counter
    mol=1;//molecular id for both
     int moln=1;//molecular id EPYC1 specific 
     int molnr=1;//Rubisco specific

    for(int i=1;i<=N;i++){//assign a molecular id

    if(idd[i]==1){
    mcnt++;

    molid[i][idd[i]]=mol;
 
        if(mcnt%5==0){
          emol=moln;
          moln++;// this will provide the number of EPYC1 molecules
          mol++;
          mcnt=0;
        }
    }
    if (idd[i]==2 || idd[i]==3){
      mcnt++;
   
     molid[i][idd[i]]=mol;

    if(mcnt%9==0){
          rmol=molnr;
          molnr++;// number of rubisco molecules.
          mol++;
          mcnt=0;
        }

    }
    }
  
             
ttdense=0;
      ttdenseR=0;
     for (t=0; t<n_steps; t++){//loop over time here
        
 int cnt=0;
  int tcnt=0;
  int cnt1=0;
         for(i=1;i<=N;i++){//zero every time step.

            tbonds[i]=0;//total bonds
             bondr[i]=0;//bonds molecule EPYC1 or Rubisco
        
               stickerbond[i]=0;

             for (j=1;j<=N;j++){
  
             bind[i][j]=0;// the indicator function between molecules. 
      

    }
         }
         

         for (i=1; i<=N-1; i++) {//loop over N's to construct distancs
			for (j=i+1; j<=N; j++) {

				double xr = xx[t][i]-xx[t][j];
// these are periodic boundaries.
				if(xr>2*Lx*0.5) {
					 xr-=2*Lx;}
				else if(xr<(-2*Lx*0.5)) {
					 xr=2*Lx+xr;}

				double yr = yy[t][i]-yy[t][j];

				if(yr>2*Ly*0.5) {
					 yr-=2*Ly;}
				else if(yr<(-2*Ly*0.5)) {
					 yr=2*Ly+yr;}

           double zr = zz[t][i]-zz[t][j];

				if(zr>2*Lz*0.5) {
					 zr-=2*Lz;}
				else if(zr<(-2*Lz*0.5)) {
					 zr=2*Lz+zr;}

				double r2 = xr*xr+yr*yr+zr*zr,r = sqrt(r2);

        //Here we will count all of the times type 2 and type 1 are within the cutoff.
        //We will also count if two there are many to one bonds. 

    // check if there is a bond between rubisco stickers and EPYC1 stickers

				if (idd[i]==1 & idd[j]==2 & r<rcut ) {// EPYC1 bead(1) and RUBISCO binding site (2) and between different bonds.
            tbonds[i]++;//total bonds
            tbonds[j]++;
            tbound++;
            tcnt++;

            //sticker bonds per molecule
            stickerbond[molid[j][idd[j]]]++;
             stickerbond[molid[i][idd[i]]]++;
            // Now compute molecualr bonds and unique bonds. 
            if(bind[molid[j][idd[j]]][molid[i][idd[i]]]!=1 ){//if the indicator function is not one already then this is a new bond. 

              bind[molid[j][idd[j]]][molid[i][idd[i]]]=1; //set the bond vector between a rubisco molecule and a EPYC1 to 1.
               bind[molid[i][idd[i]]][molid[j][idd[j]]]=1;// two molecules are bound

                 bondr[molid[j][idd[j]]]++;// the number of bonds that a rubisco molecule has at a time without double counting the same molecules twice.
              bondr[molid[i][idd[i]]]++;//the number of bonds for EPYC1

              vvlist[molid[j][idd[j]]][bondr[molid[j][idd[j]]]]=molid[i][idd[i]];//the molecular ID's for each unique bond on Rubisco or EPYC1
              vvlist[molid[i][idd[i]]][bondr[molid[i][idd[i]]]]=molid[j][idd[j]];

              
              

            }

 
				}
       
//things for cluster code.
         if(bondr[molid[j][idd[j]]]>1){nn_type[molid[j][idd[j]]]=1;}
        
		}

        if(bondr[molid[i][idd[i]]]>1) {nn_type[molid[i][idd[i]]]=1;}
	}//loop over N

// now compute the total number of bonds. 
  int ttbonds=0;

for (j=1; j<=N; j++) {

ttbonds+=tbonds[j];


}


double btot=0;
double etot=0;
double cbtot=0;
double cetot=0;


for (int k=1; k<(mol);k++){//loop over the molecules and plot the bonds and the 
 for (int p=1;p<=bondr[k];p++){

  if(p==bondr[k]){

if(k>emol){
  if(bondr[k]>0){
btot+=bondr[k];
cbtot++;
  }
}
if(k<emol){
  if(bondr[k]>0){
etot+=bondr[k];
cetot++;
  }
}
  }
 }
 }



                       //Step1: determine the instantaneous interface
                       double nslice=100;
                       double nslice1=100;
                        double dx=((2*Lx)/nslice);//The box goes from negatuve Lx to Lx
                        vector<double> h(nslice);

                         double dx1=((2*Lx)/nslice1);
                        vector<double> h1(nslice1);
                    
                        double ad=4*Lz*Ly*dx;
                           double ad1=4*Lz*Ly*dx1;
                        double ldx;
                        double rdx;

                          double ldx1;
                        double rdx1;


                        for (int j=0;j<nslice;j++){
                            h[j]=0;
                            
                          }
                           for (int j=0;j<nslice1;j++){
                            h1[j]=0;
                            
                          }
// hard coded interface here. for a check from a specific simulation. 
ldx1=-168.00+Lx;

rdx1=156.00+Lx;
           // //compute clusters
             clusterHist();
          

              frame++;
			
// print a file that shows the cluster numbers which molecules are in each cluster. 

              fprintf(x_out3,"%s\n","ITEM: TIMESTEP");
        fprintf(x_out3,"%d\n",(int)(t+1));
        fprintf(x_out3,"%s\n","ITEM: NUMBER OF ATOMS");
        fprintf(x_out3,"%d\n",N);
        fprintf(x_out3,"%s\n","ITEM: BOX BOUNDS pp pp pp");
        fprintf(x_out3,"%lf %lf\n",-Lx,Lx);
        fprintf(x_out3,"%lf %lf\n",-Ly,Ly);
        fprintf(x_out3,"%lf %lf\n",-Lz,Lz);
        fprintf(x_out3,"%s\n","ITEM: ATOMS id type x y z vx vy");



                       
                for (j=0; j<N; j++) {
                    
                  
   
             fprintf(x_out3,"%d %d %lf %lf %lf %d %d\n",j,idd[j],xx[t][j]-Lx,yy[t][j]-Ly,zz[t][j]-Lz,molid[j][idd[j]],cluster_number[molid[j][idd[j]]]);
     }




                    //now compute the molecular bonds in the dilute phase
                  int core=1;
                for (j=1; j<=N; j++) {


                

                //next we want to ask if this particle is outside of the 
                    if(xx[t][j]< ldx1 & idd[j]==1 & molid[j][idd[j]]>=core){
                        fprintf(x_out,"%d %d %d %lf %lf %lf %d %d %d %d\n",t,j,idd[j], xx[t][j]-Lx,ldx-Lx,rdx-Lx,molid[j][idd[j]],bondr[molid[j][idd[j]]],cluster_number[molid[j][idd[j]]], stickerbond[molid[j][idd[j]]]);
                        core=molid[j][idd[j]]+1;
                        
                    }

                    if(xx[t][j]> rdx1 & idd[j]==1 & molid[j][idd[j]]>=core){

                      fprintf(x_out,"%d %d %d %lf %lf %lf %d %d %d %d\n",t,j,idd[j], xx[t][j]-Lx,ldx-Lx,rdx-Lx,molid[j][idd[j]],bondr[molid[j][idd[j]]],cluster_number[molid[j][idd[j]]], stickerbond[molid[j][idd[j]]]);

                     core=molid[j][idd[j]]+1;
                   
                    }

                    if(xx[t][j]< ldx1 & idd[j]==3 & molid[j][idd[j]]>=core){
                        fprintf(x_out2,"%d %d %d %lf %lf %lf %d %d %d %d\n",t,j,idd[j], xx[t][j]-Lx,ldx1-Lx,rdx1-Lx,molid[j][idd[j]],bondr[molid[j][idd[j]]],cluster_number[molid[j][idd[j]]], stickerbond[molid[j][idd[j]]]);
                        core=molid[j][idd[j]]+1;
                        
                    }

                    if(xx[t][j]> rdx1 & idd[j]==3 & molid[j][idd[j]]>=core){

                      fprintf(x_out2,"%d %d %d %lf %lf %lf %d %d %d %d\n",t,j,idd[j], xx[t][j]-Lx,ldx1-Lx,rdx1-Lx,molid[j][idd[j]],bondr[molid[j][idd[j]]],cluster_number[molid[j][idd[j]]], stickerbond[molid[j][idd[j]]]);

                     core=molid[j][idd[j]]+1;
                    
                    }
                
                  
            
     }




     core=1;

                for (j=1; j<=N; j++) {

                    //next we want to ask if this sticker of EPYC1 is in the slab.
                    if(xx[t][j]> ldx & xx[t][j]< rdx & idd[j]==1 & molid[j][idd[j]]>=core){
                        fprintf(x_out4,"%d %d %d %lf %lf %lf %d %d %d %d\n",t,j,idd[j], xx[t][j]-Lx,ldx-Lx,rdx-Lx,molid[j][idd[j]],bondr[molid[j][idd[j]]],cluster_number[molid[j][idd[j]]], stickerbond[molid[j][idd[j]]]);
                    
                        //
                        core=molid[j][idd[j]]+1;
                        
                    }

            

                    if(xx[t][j]>ldx1 & xx[t][j]<rdx1 & idd[j]==3 & molid[j][idd[j]]>=core){
                        fprintf(x_out5,"%d %d %d %lf %lf %lf %d %d %d %d\n",t,j,idd[j], xx[t][j]-Lx,ldx1-Lx,rdx1-Lx,molid[j][idd[j]],bondr[molid[j][idd[j]]],cluster_number[molid[j][idd[j]]], stickerbond[molid[j][idd[j]]]);
                        core=molid[j][idd[j]]+1;
                     
                    }
                     if(xx[t][j]> ldx & xx[t][j]< rdx & idd[j]==1){

                         ttdense+= tbonds[j];

                     }

                       if(xx[t][j]> ldx & xx[t][j]< rdx & idd[j]==2){

                         ttdenseR+= tbonds[j];

                     }





                }
                //print the total number of bonds for Rubisco and EPYC1

                fprintf(output2,"%d %d %d\n",t+1,ttdense/(t+1),ttdenseR/(t+1));





}//time

    for(int atom=1; atom<=mol; atom++){
		if(ClusterSize[1][atom]!=0){
			fprintf(cluster4,"%d %e\n",atom, (double)(ClusterSize[1][atom])/(frame));
		}
	}
	

}//main
// Everything below is for the cluster code. 
void clusterHist(){
	int ncluster;
	int type,k,atom;
	
	cluster_id=0;

	for(k=1; k<=mol; k++) check[k]=0;//loop over molecules
    type=1;
	
    for(atom=1; atom<=mol; atom++){
 
        ncluster=0;
                    
        if(nn_type[atom]==type){
            if(check[atom]!=1) {
                cluster_id=(cluster_id+1);// this is its cluster ID. 
                ncluster=clustercheck(atom,ncluster,type);//this is the cluster size.

            }	
        }
        
        ClusterSize[type][ncluster]++;//this might potentiall help me see the structures I want.
        
        if(type==1){
            ave_size+=ncluster;




            if(large_cluster<ncluster) {//computes the largest cluster. 
                 large_cluster=ncluster;
                 if(type ==1) large_cluster_id=cluster_id;
             }
        }
        
        if(ncluster!=0) n_clusters++;
    }
    ave_size=ave_size/n_clusters;
	
}

int clustercheck(int atom, int ncluster, int type){
	int j,jj;
	
	if(check[atom]!=1){
		ncluster++;
		check[atom]=1;
		cluster_number[atom]=cluster_id;
		
		for(j=1; j<=bondr[atom]; j++){
			jj=vvlist[atom][j];
			if(nn_type[jj]==type) ncluster=clustercheck(jj,ncluster,type);
		}	
	}	
	
	return ncluster;
}
