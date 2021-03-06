#include "globals.h"
void lagrange_get_weights( double , double , double* ) ;
void spline_get_weights( double , double , double* ) ;
void add_segment( int ) ;

///////////////////////////////////
// Add all particles to the grid //
///////////////////////////////////
void charge_grid( ) {

  int i, j;
  
#pragma omp parallel for
  for ( i=0 ; i<M ; i++ ) 
    for ( j=0 ; j<ntypes ; j++ ) 
      rho[j][i] = 0.0 ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ ) {
    rhogb[i] = rhoda[i] = rhoha[i] = rhodb[i] = rhohb[i] = rhop[i] = smrhop[i] = rhoga[i] = 0.0 ;
    for ( j=0 ; j<nthreads ; j++ ) 
      rhogb_t[j][i] = rhoda_t[j][i] = rhoha_t[j][i] = rhodb_t[j][i] = rhohb_t[j][i] 
        = rhoga_t[j][i] = rhop_t[j][i] = 0.0 ;
  }

  ////////////////////////////////////////////////////
  // Add segments to each processors density fields //
  ////////////////////////////////////////////////////
#pragma omp parallel for
  for ( i=0 ; i<nstot ; i++ ) {
    if ( tp[i] != -1 ){
     add_segment( i ) ;
     }
  }

#pragma omp parallel for private(j)
  for ( i=0 ; i<M ; i++ ) {
    for ( j=0 ; j<nthreads ; j++ ) {
      rhoda[i] += rhoda_t[j][i] ;
      rhodb[i] += rhodb_t[j][i] ;
      rhoha[i] += rhoha_t[j][i] ;
      rhohb[i] += rhohb_t[j][i] ;
      rhop[i] += rhop_t[j][i] ;
      rhoga[i] += rhoga_t[j][i] ;
      rhogb[i] += rhogb_t[j][i] ;
    }

    rhot[i] = rhogb[i] + rhoda[i] + rhodb[i] + rhoha[i] + rhohb[i]  + rhoga[i];

    rho[0][i] = rhoda[i] + rhoha[i] + rhoga[i] ;
    rho[1][i] = rhodb[i] + rhohb[i] + rhogb[i];
    rho[2][i] = rhop[i] ;
  }

}




//////////////////////////////////////////////
// Adds the segment associated with particle //
// "id" to the PME grid using Lagrange      //
// interpolation scheme. From JCP V103 3668 //
//////////////////////////////////////////////
void add_segment( int id ) {

  int tid = omp_get_thread_num() ;

  if ( spline_weights + lagrange_weights == 2 ||
       spline_weights + lagrange_weights == 0 )
    die("Invalid spline/lagrange weight function") ;

  int j, g_ind[Dim] , ix, iy, iz, nn[Dim] , Mindex, grid_ct; 
  //double **W , gdx , W3;
  
  double gdx , W3;

  // W = ( double** ) calloc( Dim , sizeof( double* ) );



  ///////////////////////////////////////////////
  // First, determine the relevant weights for //
  // all grid points in all directions.        //
  ///////////////////////////////////////////////
  for ( j=0 ; j<Dim ; j++ ) {

   // W[j] = ( double* ) calloc( pmeorder+1 , sizeof( double ) );


    // Distance to nearest grid point if even //
    if ( pmeorder % 2 == 0 ) {
      g_ind[j] = int( ( x[id][j] + 0.5 * dx[j] ) / dx[j] ) ;

      gdx = x[id][j] - double( g_ind[j] ) * dx[j] ;
    }
 

    // Distance to nearest mid-point between grid points if odd //
    else {
      g_ind[j] = int( ( x[id][j]  ) / dx[j] ) ;
      if ( g_ind[j] >= Nx[j] ){
       //cout<<id<<" "<<g_ind[j]<<endl;
//       printf("Using the boundary grid point!\n") ;

       g_ind[j] = Nx[j]-1;
      }
      gdx = x[id][j]  - ( double( g_ind[j] ) + 0.5 ) * dx[j] ;
    }


    /////////////////////////////////////////
    // Get the weights for each grid point //
    /////////////////////////////////////////
    if ( lagrange_weights ) 
      lagrange_get_weights( gdx , dx[j] , W_tsn[id][j] );

    else if ( spline_weights ) 
      spline_get_weights( gdx , dx[j] , W_tsn[id][j] );

  }//for ( j=0 ; j<3...



  ////////////////////////////////////////////////////
  // Assign the weights to all relevant grid points //
  ////////////////////////////////////////////////////
  grid_ct = 0 ;
  
  ///////////////////////////////////////////
  // 3D version of particle-to-mesh scheme //
  ///////////////////////////////////////////
  if ( Dim == 3 ) {
    for ( ix = 0 ; ix < pmeorder+1 ; ix++ ) {
      
      nn[0] = g_ind[0] + ix - ( pmeorder/2 + pmeorder % 2 );
      
      if ( nn[0] < 0 ) nn[0] += Nx[0] ;
      else if ( nn[0] >= Nx[0] ) nn[0] -= Nx[0] ;
  
      for ( iy = 0 ; iy < pmeorder+1 ; iy++ ) {
  
        nn[1] = g_ind[1] + iy - ( pmeorder/2 + pmeorder % 2 ) ;
        
        if ( nn[1] < 0 ) nn[1] += Nx[1] ;
        else if ( nn[1] >= Nx[1] ) nn[1] -= Nx[1] ;
  
        for ( iz = 0 ; iz < pmeorder+1 ; iz++ ) {
  
          nn[2] = g_ind[2] + iz - ( pmeorder/2 + pmeorder % 2 ) ;
          
          if ( nn[2] < 0 ) nn[2] += Nx[2] ;
          else if ( nn[2] >= Nx[2] ) nn[2] -= Nx[2] ;
  
  
          // stack() returns index in [0,M]
          Mindex = stack( nn ) ;
  
          if ( Mindex >= M ) {
            char nm[40] ;
            sprintf(nm, "%d Index = %d out of range, particle %d %lf %lf %lf\n" , 
                step, Mindex, id , x[id][0], x[id][1], x[id][2] ) ;

            die(nm) ;
          }
  
          W3 = W_tsn[id][0][ix] * W_tsn[id][1][iy] * W_tsn[id][2][iz] / gvol ;
  
          if ( id < nsD && tp[id] == 0 )
            rhoda_t[tid][ Mindex ] += W3 / CG_ratio ;
          else if ( id < nsD && tp[id] == 1 )
            rhodb_t[tid][ Mindex ] += W3 / CG_ratio ;
          else if ( id < ( nsD + nsA ) )
            rhoha_t[tid][ Mindex ] += W3 / CG_ratio ;
          else if ( id < ( nsD + nsA + nsB ) )
            rhohb_t[tid][ Mindex ] += W3 / CG_ratio ;
          else if ( tp[id] == 2 )
            rhop_t[tid][ Mindex ] += W3 / CG_ratio ;
          else if ( tp[id] == 0 )
            rhoga_t[tid][ Mindex ] += W3 / CG_ratio ;
          else if( tp[id] == 1)
	    rhogb_t[tid][ Mindex ] += W3 / CG_ratio ;
	  else {
            char nm[40] ;
            sprintf(nm, "Invalid partic type. tp[%d] = %d\n" , id, tp[id] ) ;
            die(nm);
          }
  
          grid_inds[ id ][ grid_ct ] = Mindex ;
          grid_W[ id ][ grid_ct ] = W3 ;
  
          grid_ct++ ;
  
        }
      }
    }
  }

  ///////////////////////////////////////////
  // 2D version of particle-to-mesh scheme //
  ///////////////////////////////////////////
  else if ( Dim == 2 ) {
    for ( ix = 0 ; ix < pmeorder+1 ; ix++ ) {
      
      nn[0] = g_ind[0] + ix - ( pmeorder/2 + pmeorder % 2 );
      
      if ( nn[0] < 0 ) nn[0] += Nx[0] ;
      else if ( nn[0] >= Nx[0] ) nn[0] -= Nx[0] ;
  
      for ( iy = 0 ; iy < pmeorder+1 ; iy++ ) {
  
        nn[1] = g_ind[1] + iy - ( pmeorder/2 + pmeorder % 2 ) ;
        
        if ( nn[1] < 0 ) nn[1] += Nx[1] ;
        else if ( nn[1] >= Nx[1] ) nn[1] -= Nx[1] ;
  
        Mindex = stack( nn ) ;
  
        if ( Mindex >= M ) {
          char nm[40] ;
          sprintf(nm, "%d Index = %d out of range, particle %d %lf %lf %lf\n" , 
              step, Mindex, id , x[id][0], x[id][1], x[id][2] ) ;

          die(nm) ;
        }
  
        W3 = W_tsn[id][0][ix] * W_tsn[id][1][iy] / gvol ;
          
        if ( id < nsD && tp[id] == 0 )
          rhoda_t[tid][ Mindex ] += W3 / CG_ratio ;
        else if ( id < nsD && tp[id] == 1 )
          rhodb_t[tid][ Mindex ] += W3 / CG_ratio ;
        else if ( id < ( nsD + nsA ) )
          rhoha_t[tid][ Mindex ] += W3 / CG_ratio ;
        else if ( id < ( nsD + nsA + nsB ) )
          rhohb_t[tid][ Mindex ] += W3 / CG_ratio ;
        else if ( tp[id] == 2 )
          rhop_t[tid][ Mindex ] += W3 / CG_ratio ;
        else if ( tp[id] == 0 )
          rhoga_t[tid][ Mindex ] += W3 / CG_ratio ;
        else if (tp[id] == 1 )
	   rhogb_t[tid][ Mindex ] += W3 / CG_ratio ;
	else {
          char nm[40] ;
          sprintf(nm, "Invalid partic type. tp[%d] = %d\n" , id, tp[id] ) ;
          die(nm);
        }
  
        grid_inds[ id ][ grid_ct ] = Mindex ;
        grid_W[ id ][ grid_ct ] = W3 ;
  
        grid_ct++ ;
  
        
      }
    }
  }

  ///////////////////////////
  // Free allocated memory //
  ///////////////////////////
//  for ( j=0 ; j<Dim ; j++ ) 
  //  free( W[j] ) ;

 // free(W) ;

}// End lagrange add charge







///////////////////////////////////////////////////////////////////
// Taken directly from appendix in Petersen JCP V103 3668 (1995) //
///////////////////////////////////////////////////////////////////
void lagrange_get_weights( double dx , double H , double *W ) {
  double norm, dx2, dx4, H2, H4;

  if ( pmeorder == 0 ) 
    W[0] = 1. ;

  else if ( pmeorder == 2 ) {
    norm = 2.*H*H ;
    dx2 = dx * dx ;
    W[0] = (dx2 - H * dx ) / norm ;
    W[1] = ( -2.*dx2 + 2*H*H ) / norm ;
    W[2] = (dx2 + H * dx ) / norm ;
  }

  else if ( pmeorder == 3 ) {
    norm = 48.*H*H*H ;
    dx2 = dx * dx ;
    W[0] = ( -8.*dx2*dx + 12.*H*dx2 + 2.*H*H*dx - 3.*H*H*H ) / norm ;
    W[1] = ( 24.*dx2*dx - 12.*H*dx2 - 54.*H*H*dx + 27.*H*H*H ) / norm ;
    W[2] = ( -24.*dx2*dx - 12.*H*dx2 + 54.*H*H*dx + 27.*H*H*H ) / norm ;
    W[3] = ( 8.*dx2*dx + 12.*H*dx2 - 2.*H*H*dx - 3.*H*H*H ) / norm ;
  }

  else if ( pmeorder == 4 ) {
    H2 = H * H;
    H4 = H2 * H2 ;
    dx2 = dx*dx;
    dx4 = dx2*dx2;

    norm = 24.0 * H4;

    W[0] = ( dx4 - 2.*H*dx2*dx - H2*dx2 + 2*H2*H*dx ) / norm ;
    W[1] = ( -4.*dx4 + 4.*H*dx2*dx + 16.*H2*dx2 - 16.*H2*H*dx ) / norm ;
    W[2] = ( 6.*dx4 - 30.*H2*dx2 + 24.*H4 ) / norm ;
    W[3] = ( -4.*dx4 - 4.*H*dx2*dx + 16.*H2*dx2 + 16.*H2*H*dx ) / norm ;
    W[4] = ( dx4 + 2.*H*dx2*dx - H2*dx2 - 2*H2*H*dx ) / norm ;

  }

  else if ( pmeorder == 5 ) {
    // I think there is a bug in here! //
    cout << "Possible bug in lagrange interpolation with pmeorder==5\n" ;
    H2 = H * H ;
    H4 = H2 * H2 ;

    dx2 = dx * dx;
    dx4 = dx2 * dx2 ;

    norm = 3840. * H4 * H;

    W[0] = ( -32.*dx4*dx + 80.*H*dx4 + 80.*H2*dx2*dx - 200.*H2*H*dx2
             - 18.*H4*dx + 45.*H4*H ) / norm ;
    W[1] = ( 160.*dx4* - 240.*H*dx4 - 1040.*H2*dx2*dx + 1560.*H2*H*dx2
           + 250.*H4*dx - 375.*H4*H ) / norm ;
    W[2] = ( -320.*dx4*dx + 160.*H*dx4 + 2720.*H2*dx2*dx - 1360.*H2*H*dx2
           - 4500.*H4*dx + 2250*H4*H ) / norm ;

    W[2] = ( 320.*dx4*dx + 160.*H*dx4 - 2720.*H2*dx2*dx - 1360.*H2*H*dx2
           + 4500.*H4*dx + 2250*H4*H ) / norm ;
    W[1] = ( -160.*dx4* - 240.*H*dx4 + 1040.*H2*dx2*dx + 1560.*H2*H*dx2
           - 250.*H4*dx - 375.*H4*H ) / norm ;
    W[0] = ( 32.*dx4*dx + 80.*H*dx4 - 80.*H2*dx2*dx - 200.*H2*H*dx2
             + 18.*H4*dx + 45.*H4*H ) / norm ;
  }


  else {
    cout << "PME order is " << pmeorder << endl;
    die("PME not set up for this interpolation order!\n");
  }


}


void spline_get_weights( double dx , double H , double *W ) {
  double sx = dx / H ;
  
  double sx2, sx3, sx4, sx5;
  double scale = double( M ) / V ;

  if ( pmeorder == 0 ) 
    W[0] = 1. * scale ;

  else if ( pmeorder == 1 ) {
    W[0] = ( 1. - 2. * sx ) / 2.0 ;
    W[1] = ( 1. + 2. * sx ) / 2.0 ;
  }

  else if ( pmeorder == 2 ) {
    sx2 = sx * sx ;

    W[0] = (1. - 4. * sx + 4.*sx2 ) / 8. ;
    W[1] = (3. - 4. * sx2 ) / 4. ;
    W[2] = (1. + 4. * sx + 4.*sx2 ) / 8. ;

  }

  else if ( pmeorder == 3 ) {
    sx2 = sx * sx ;
    sx3 = sx2 * sx ;

    W[0] = ( 1. - 6.*sx + 12.*sx2 - 8.*sx3 ) / 48. ;
    W[1] = ( 23. - 30.*sx - 12.*sx2 + 24.*sx3 ) / 48. ;
    W[2] = ( 23. + 30.*sx - 12.*sx2 - 24.*sx3 ) / 48. ;
    W[3] = ( 1. + 6.*sx + 12.*sx2 + 8.*sx3 ) / 48. ;
  }

  else if ( pmeorder == 4 ) {
    sx2 = sx * sx ;
    sx3 = sx2 * sx2 ;
    sx4 = sx2 * sx2 ;

    W[0] = (1. - 8.*sx + 24.*sx2 - 32.*sx3 + 16.*sx4 ) / 384.;
    W[1] = (19. - 44.*sx + 24.*sx2 + 16.*sx3 - 16.*sx4 ) / 96. ;
    W[2] = (115. - 120.*sx2 + 48.*sx4 ) / 192. ;
    W[3] = (19. + 44.*sx + 24.*sx2 - 16.*sx3 - 16.*sx4 ) / 96. ;
    W[4] = (1. + 8.*sx + 24.*sx2 + 32.*sx3 + 16.*sx4 ) / 384.;
  }

  else if ( pmeorder == 5 ) {
    sx2 = sx * sx ;
    sx3 = sx2 * sx2 ;
    sx4 = sx2 * sx2 ;
    sx5 = sx4 * sx ;

    W[0] = (1. - 10.*sx + 40.*sx2 - 80.*sx3 + 80.*sx4 - 32.*sx5 ) / 3840.;
    W[1] = (237. - 750.*sx + 840.*sx2 - 240.*sx3 - 240.*sx4 + 160.*sx5 ) / 3840.;
    W[2] = (841. - 770.*sx - 440.*sx2 + 560.*sx3 + 80.*sx4 - 160.*sx5 ) / 1920. ;
    W[3] = (841. + 770.*sx - 440.*sx2 - 560.*sx3 + 80.*sx4 + 160.*sx5 ) / 1920. ;
    W[4] = (237. + 750.*sx + 840.*sx2 + 240.*sx3 - 240.*sx4 - 160.*sx5 ) / 3840.;
    W[5] = (1. + 10.*sx + 40.*sx2 + 80.*sx3 + 80.*sx4 + 32.*sx5 ) / 3840.;

  }


  else
    die("P3M not set up for this interpolation order!\n");

}

