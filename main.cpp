#define MAIN
#include "globals.h"
void update_positions( void ) ;
void update_positions_init( void ) ;
void update_euler(void);
void update_euler_2d();
void initialize( void ) ;
void write_np(void) ;
void write_gro( void ) ;
void write_rst_gro(void );
void write_quaternions(void) ;
void write_grid( void ) ;
void forces( void ) ;
void torque(void);
void torque_2d();
double integrate( double* ) ;
void write_stress( void ) ;
void bond_stress( void ) ;
void calc_Unb( void ) ;
void calc_stress();

int main( int argc , char** argv ) {

  int i,j,k,l;

  if ( argc == 3 && !strcmp( "-nt" , argv[1] ) ) {
    nthreads = atoi( argv[2] ) ;
    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
  }
  else {
    nthreads = 1 ;
    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
  }


  initialize() ;

  write_gro( ) ;
  write_quaternions( ) ;
  // Save chi for pre-equilibration steps //
  double tmp_ang,chi_bkp = chiAB ;

  FILE *otp ;
  otp = fopen( "data.dat" , "w" ) ;

  printf("Entering main loop!\n") ; fflush( stdout ) ;
  for ( step = 0 ; step < nsteps ; step++ )  {

    if ( step < pre_equil_steps )
      chiAB = 0.0 ;
    else
      chiAB = chi_bkp ;

   

    forces() ;
    if(sigma>0){
    if(Dim ==3 )torque();
     else torque_2d();
    }
  //  exit(1);

    if(step >0 || rst_para == 1)
   	update_positions() ;
    else
    	update_positions_init() ;
    
    if(sigma>0){
       if(Dim == 3 ) update_euler();
        else update_euler_2d();
   }  


    if ( stress_freq > 0 && step % stress_freq == 0 ) {
      calc_stress() ;

      for ( j=0 ; j<Dim ; j++ ){ 
        for ( k=0 ; k<Dim ; k++ ) {
          sts_buf[buff_ind][j][k]= Rg3*Ptens[j][k];//( j<Dim ? Stress_bonds[j][k]:0.0) ;
           sts_buf_pp[buff_ind][j][k] = Rg3*Stress_PP[j][k];
	   sts_buf_ng[buff_ind][j][k] = Rg3*Stress_Ng[j][k];

	/*for ( k=0 ; k<nP;  k++ ){
		//euler_adot[k][j] = euler_q[k][j];	
		sts_buf[buff_ind][j][k+Dim] = euler_q[k][j];
	}*/
	}
      }
      buff_ind++ ;
    }
   
    if(step % sample_freq == 0){
	write_np();
    }  

    
    if ( step > sample_wait && step % sample_freq == 0 ) {
      fftw_fwd( rho[0] , ktmp ) ;
      for ( i=0 ; i<M ; i++ ) {
        avg_sk[0][i] += ktmp[i] * conj(ktmp[i]) ;
      }

      if ( nP > 0 ) {

        fftw_fwd( rho[2] , ktmp ) ;
        for ( i=0 ; i<M ; i++ ) {
          avg_sk[2][i] += ktmp[i] * conj( ktmp[i] ) ;
        }
      }

      num_averages += 1.0 ;
    }


    if ( step % print_freq == 0 || step == nsteps-1 ) {
      printf("step %d of %d  Ubond: %lf\n" , step , nsteps , Ubond ) ;
      fflush( stdout ) ;
      write_gro() ;
      write_rst_gro();
      write_quaternions();


      if ( stress_freq > 0 )
        write_stress() ;


      write_grid_data( "rhoda.dat" , rhoda ) ;
      write_grid_data( "rhodb.dat" , rhodb ) ;
      
      if ( nA > 0.0 )
        write_grid_data( "rhoha.dat" , rhoha ) ;

      if ( nB > 0.0 )
        write_grid_data( "rhohb.dat" , rhohb ) ;

      if ( nP > 0.0 ) 
        write_grid_data( "rhop.dat" , rhop ) ;

      if ( step > sample_wait ) {
        for ( i=0 ; i<M ; i++ ) 
          ktmp2[i] = avg_sk[0][i] / num_averages ;
        write_kspace_data( "avg_sk_A.dat" , ktmp2 ) ;
        
        if ( nP > 0 ) {
          for ( i=0 ; i<M ; i++ )
            ktmp2[i] = avg_sk[2][i] / num_averages ;
          write_kspace_data( "avg_sk_np.dat" , ktmp2 ) ;
        }
      }

      calc_Unb() ;

      fprintf( otp , "%d %lf %lf %lf %lf %lf %lf %lf\n" , step , Ubond , U_chi_gg, U_kappa_gg,U_chi_pg,U_kappa_pg,U_kappa_pp , Utt) ;
      fflush( otp ) ;

    }// if step % print_Freq == 0

  }

  fclose( otp ) ;

  return 0 ;

}
