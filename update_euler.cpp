#include "globals.h"
void calc_A();
void rot_noise(double**);
double crct_q(double* , double* );


void update_euler_2d(){

         int i, j,k, cent_ind,gind, ind;
	         int sites_per_gnp = 1 + ngb_per_partic*(Ngb+1)+ng_per_partic * ( 1 + Ng );
		    double vect[Dim],tmp_gn[nP] , de_theta,**tmp_q,lamb,tmp_x,tmp_y;

       for(i=0 ; i<nP ; i++){
	

            tmp_gn[i] = sqrt(2.0*Diff_rot*delt)* gasdev2();               

       }

#pragma omp parallel for private(k,vect,tmp_x,tmp_y,de_theta,j,ind,cent_ind)
      for(i=0 ; i<nP ; i++){
         cent_ind =nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ; 
	 de_theta =Diff_rot*delt*real_trq[i][0]   + tmp_gn[i]; 
         
	 for( j =0 ; j<(ngb_per_partic+ ng_per_partic) ; j++){
             ind = cent_ind  + 1+ (j >= ng_per_partic ? ng_per_partic*(Ng+ 1) + (j-ng_per_partic)*(Ngb +1)    : j*(Ng+ 1) ) ;
           
	    if(tp[ind] != -1 ){
		cout<<"wrong xind in updating euler rot"<<endl;
		exit(1);
	    }

	    for(k=0 ; k<Dim ; k++){
	       vect[k] = (x[ind][k]-x_bac[cent_ind][k]);
               if(vect[k] > L[k]/2.0){
		  vect[k] -= L[k];
	       }
	       else if(vect[k] < -L[k]/2.0){
		 vect[k] += L[k]; 
	       }

             }//k

            tmp_x = cos(de_theta)*vect[0] - sin(de_theta)*vect[1];//
	    tmp_y = sin(de_theta)*vect[0] + cos(de_theta)*vect[1];

	    for(k =0 ;k < Dim ; k++){

	       if(k == 0)	{x[ind][k] = x[cent_ind][k] + tmp_x;
	       }
	       else {x[ind][k] = x[cent_ind][k] + tmp_y;
	       }
	       if(x[ind][k] < 0.0 )
	             x[ind][k] += L[k];
	       else if(x[ind][k] > L[k]) 
	             x[ind][k] -= L[k];
	       else{

	       }
	    }//k
	
	 }//j

      }//i

}

void update_euler(){

	int i, j,k, cent_ind,gind, ind;
	int sites_per_gnp = 1 + ng_per_partic * ( 1 + Ng );
	double **tmp_q,lamb;
    
        tmp_q = (double**) calloc(nP, sizeof(double*));
	for(i=0 ; i<nP ; i++)
		tmp_q[i] = (double*) calloc(4, sizeof(double));

     rot_noise(q_noise);

#pragma omp parallel for private(j,lamb)
	for(i=0 ; i<nP ; i++){
	  for(j=0;j<4;j++){

           	//get the predicted quaterns, q~//	
		tmp_q[i][j] =  euler_q[i][j] + Diff_rot*delt*trq[i][j] + q_noise[i][j];	
		
	      
	  }//j
	  lamb = crct_q(euler_q[i],tmp_q[i]); // Lagrange multiplier , lambda
	 
	  
	  for(j=0;j<4;j++){

	  	euler_q[i][j] = tmp_q[i][j] + lamb*euler_q[i][j]; //unit length correction 
	  }//j

	}//i=nP

       calc_A();


          

#pragma omp parallel for private(k,j,cent_ind,ind,gind)
       for(i = 0 ; i<nP ;i++){
           cent_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ;
	    for(j = 0; j<ng_per_partic; j++){
		ind = cent_ind + j * ( Ng + 1 ) + 1 ; 
		gind = ng_per_partic*i + j;
		matrix_trs_vect(euler_A[i],grf_bf_x[gind],x[ind]);

       		if(tp[ind] != -1) { cout<<"wrong em update !!"<<endl; exit(1);} 
		
		for(k=0 ;k<Dim; k++){
			x[ind][k] += x[cent_ind][k] ;
		      if ( x[ind][k] > L[k] )
		            x[ind][k] -= L[k] ;
			else if ( x[ind][k] < 0.0 )
			    x[ind][k] += L[k] ;


		}
	   }//ng_per_partic


       }//nP

       for(i=0 ;i<nP;i++)
       	 free(tmp_q[i]);
       free(tmp_q);

}


double crct_q(double* q, double* q_tld){

	double delta,q_qtld=0, qtld2=0;
	int i;
        for(i=0 ; i<4 ;i++){

		q_qtld += q[i]*q_tld[i];
		qtld2 += q_tld[i]*q_tld[i];

	}

	delta = q_qtld*q_qtld - qtld2 +1;

        if(delta <0 ){
	  cout<<"can not find lamb_a !"<<endl;
	  exit(1);
	}
  
  return (-q_qtld+sqrt(delta));
}
