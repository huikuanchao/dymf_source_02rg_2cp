#include "globals.h"
void torque_2d(void){

	int i,j,m,k,l, gind, t1,t2, center_ind,ind;
	double *tmp_trp,mt, mdr2, mdrp2, mdr, dr[Dim], drp[Dim],delr ;
     int sites_per_gnp = 1 +ngb_per_partic * ( 1 + Ngb ); +ng_per_partic * ( 1 + Ng );
// set all the torque to zero//

#pragma omp parallel for private(j)
   for(i=0;i<nP;i++){
	for(j = 0; j<4 ;j++)
		trq[i][j] = 0;
	for(j = 0; j<3; j++)
		real_trq[i][j] = 0;
   }


if(sites_per_gnp >0){
#pragma omp parallel for private(dr,drp,m,j,k,l,ind,gind,center_ind,mdrp2,mdr2)\
reduction(+:Ubond)
	for(i=0;i<nP;i++){
	   center_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ;	   
	   //real torque on each nP
	   for( m=0 ; m<ng_per_partic ; m++){	
	      ind = center_ind + m * ( Ng + 1 ) + 1 ;
	     // gind = ng_per_partic*i + m;
	      mdr2 = pbc_mdr2( x[ind+1] , x[ind] , dr ) ;//dist between r' and r 
	      mdrp2 = pbc_mdr2(x[ind],x[center_ind],drp) ; //dist between rcm and r
	      if(tp[ind] != -1) { cout<<"wrong trq calc !!"<<endl; exit(1);} 

		Ubond += 1.5*mdr2;	   
	      for(j=0;j<1;j++){
            	        //for(k=0 ; k<Dim ; k++){
		//	for(l=0; l<Dim ; l++){
		   		real_trq[i][j] += 3.0*drp[0]*dr[1] - 3.0*drp[1]*dr[0];	   
		//	}//l
	      	

		
		//}//k
		
	      }//j, torque dim


	    

	   }//m=ng_per_oartic
	   
	  for( m=0 ; m<ngb_per_partic ; m++){	
	      ind = center_ind + ng_per_partic*(Ng+1) +m * ( Ngb + 1 ) + 1 ;
	     // gind = ng_per_partic*i + m;
	      mdr2 = pbc_mdr2( x[ind+1] , x[ind] , dr ) ;//dist between r' and r 
	      mdrp2 = pbc_mdr2(x[ind],x[center_ind],drp) ; //dist between rcm and r
	      if(tp[ind] != -1) { cout<<"wrong trq calc !!"<<endl; exit(1);} 

		Ubond += 1.5*mdr2;	   
	      for(j=0;j<1;j++){

	        //for(k=0 ; k<Dim ; k++){
		//	for(l=0; l<Dim ; l++){
	           real_trq[i][j] += 3.0*drp[0]*dr[1] - 3.0*drp[1]*dr[0];	   
		//	}//l
	      	

		
		//}//k
		
	      }//j, torque dim


	    

	   }//m=ng_per_oartic
	   


	    


	   //get the quaternions' langevin force term, trq//

	}//nP
}//site_np>0?


}
