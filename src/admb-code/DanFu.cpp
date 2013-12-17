    
    // some parameters  for with logistnormal
    init_bounded_number LN_sigmaCSLF(lbd,ubd,phase)         
    init_bounded_number LN_phi1CSLF(lbd,ubd,phase)            
    init_bounded_number LN_taoCSLF(lbd,ubd,phase)            
    number LN_phi2CSLF;

    
    // some variables/structures associated with logistnormal
    vector LN_sigmasCSLF(1,nBinCSLFs)     //nBinCSLFs and nCSLFs: number of length bins and  number of years of obserations
    matrix LN_covmatCSLF(1,nBinCSLFs,1,nBinCSLFs)
    matrix LN_VmatCSLF(1,nBinCSLFs-1,1,nBinCSLFs-1)
    matrix LN_VinvCSLF(1,nBinCSLFs-1,1,nBinCSLFs-1)
    matrix LN_wwCSLF(1,nCSLFs,1,nBinCSLFs-1)
  
  
    // nll_logistnorm gets called depending on whether the user defined LN1, LN2 and LN3
  
               case 7: // logistnormal LN1
               {
                  LN_VmatCSLF.initialize();
  	        LN_VinvCSLF.initialize();
  	        LN_wwCSLF.initialize();
  	        
  	        LN_sigmasCSLF = LN_sigmaCSLF; 
  	        LN_covmatCSLF = covmat_logistnorm(LN_sigmasCSLF,nBinCSLFs);
  	        likeCSLF = nll_logistnorm(CSLFtruncated,predCSLFtruncated,LN_wtsCSLF,LN_covmatCSLF,LN_VmatCSLF,LN_VinvCSLF,LN_wwCSLF);
  	        
  	        
               }
               break;
               case 8: // logistnormal LN2
               {
                  LN_VmatCSLF.initialize();
  	        LN_VinvCSLF.initialize();
  	        LN_wwCSLF.initialize();
  	     
  	        LN_sigmasCSLF = LN_sigmaCSLF; 
  	        
  	        LN_covmatCSLF = covmat_logistnorm(LN_sigmasCSLF,LN_phi1CSLF,nBinCSLFs);
  	        likeCSLF = nll_logistnorm(CSLFtruncated,predCSLFtruncated,LN_wtsCSLF,LN_covmatCSLF,LN_VmatCSLF,LN_VinvCSLF,LN_wwCSLF);
               
               }
               break;
               case 9: // logistnormal LN3
               {
                  LN_VmatCSLF.initialize();
  	        LN_VinvCSLF.initialize();
  	        LN_wwCSLF.initialize();
  	        LN_phi2CSLF = -1+(2-fabs(LN_phi1CSLF))*LN_taoCSLF;
  	        LN_sigmasCSLF = LN_sigmaCSLF;  
  	        LN_covmatCSLF = covmat_logistnorm(LN_sigmasCSLF,LN_phi1CSLF,LN_phi2CSLF,nBinCSLFs);
  	    
  	        likeCSLF = nll_logistnorm(CSLFtruncated,predCSLFtruncated,LN_wtsCSLF,LN_covmatCSLF,LN_VmatCSLF,LN_VinvCSLF,LN_wwCSLF);
               
               }
               break;






 // functions that implements various bits of logist_nommal 

  void suppresszeroes(dmatrix& compdat,_CONST double eps_const)
  {
    compdat = compdat+eps_const;
    for(i=compdat.rowmin();i<=compdat.rowmax();i++) 
      compdat(i) = compdat(i)/sum(compdat(i));
	
  }
  	
 dvariable zerofun(dvariable x, _CONST double r_const) 
   {
       RETURN_ARRAYS_INCREMENT();
       dvariable z;
       if(x>=r_const)
         z=x;
       else
         z=r_const/(2-(x/r_const));
       RETURN_ARRAYS_DECREMENT();
       return (z);
   }
 	
  
  dvar_vector getrho(dvariable& phi1,int kk)
  {
    //calculation of AR(1) acf for  LN2 
    RETURN_ARRAYS_INCREMENT();
    dvar_vector rho(1,kk-1);
    for(int i=1;i<=kk-1;i++)
      rho(i)= pow(phi1,i);
    RETURN_ARRAYS_DECREMENT();   
    return(rho);
  }

  dvar_vector getrho(dvariable& phi1,dvariable& phi2,int kk)
  {
    //calculation of AR(2) acf for  LN3 
    RETURN_ARRAYS_INCREMENT();
    
    if((phi2<-1) || (1-fabs(phi1) < phi2)){
      cout<<"Unstable phi paramete\n";
      exit(1);
    }
    dvar_vector acvect(0,kk);
    acvect=1;
    acvect(1) = phi1/(1-phi2);
    for(int i =2;i<=kk;i++){
      acvect[i] = phi1*acvect(i-1)+phi2*acvect(i-2);
    }
    dvar_vector rho(1,kk);
    for(int i=1;i<=kk;i++){
      rho[i] =acvect(i);
    }
    
    RETURN_ARRAYS_DECREMENT();   
    return(rho);
  }
  
  // covariance matrix for LN1
  dvar_matrix covmat_logistnorm(dvar_vector& sigma, int nBin)
  {
    //calculation of covriance matrix for LN2 
    RETURN_ARRAYS_INCREMENT();
    if(sigma.size() != nBin) {
      cout<<"covmat_LN1: sigma.size() != nBin   \n";
      exit(1);
    }
    dvar_matrix covmat = identity_matrix(1, nBin);
    for(int i=1;i<=nBin;i++)
        covmat(i,i)=sigma(i)*sigma(i);
    RETURN_ARRAYS_DECREMENT();
    return (covmat);
  }
 
 // covariance matrix for LN2
  dvar_matrix covmat_logistnorm(dvar_vector& sigma,dvariable& phi1, int nBin)
  {
    //calculation of covriance matrix for LN2 
    RETURN_ARRAYS_INCREMENT();
    if(sigma.size() != nBin) {
      cout<<"covmat_logistnorm: sigma.size() != nBin   \n";
      exit(1);
    }
    dvar_matrix covmat = identity_matrix(1, nBin);
    dvar_vector rhovec = getrho(phi1,nBin);
    for(int i=1;i<=nBin;i++){
      for(int j=1;j<=nBin;j++){
        if(i!=j) covmat(i,j)=rhovec(abs(i-j));
      }
    } 
    for (int i=1;i<=nBin;i++)
      covmat.rowfill(i,elem_prod(extract_row(covmat,i),sigma));
    for (int j=1;j<=nBin;j++)
      covmat.colfill(j,elem_prod(extract_column(covmat,j),sigma));
    
    RETURN_ARRAYS_DECREMENT();
    return (covmat);
  }
  
  // covariance matrix for LN3
    dvar_matrix covmat_logistnorm(dvar_vector& sigma,dvariable& phi1, dvariable& phi2,int nBin)
    {
      //calculation of covriance matrix for LN2 
      RETURN_ARRAYS_INCREMENT();
      if(sigma.size() != nBin) {
        cout<<"covmat_logistnorm: sigma.size() != nBin   \n";
        exit(1);
      }
      dvar_matrix covmat = identity_matrix(1, nBin);
      dvar_vector rhovec = getrho(phi1,phi2,nBin);
      for(int i=1;i<=nBin;i++){
        for(int j=1;j<=nBin;j++){
          if(i!=j) covmat(i,j)=rhovec(abs(i-j));
        }
      } 
      for (int i=1;i<=nBin;i++)
        covmat.rowfill(i,elem_prod(extract_row(covmat,i),sigma));
      for (int j=1;j<=nBin;j++)
        covmat.colfill(j,elem_prod(extract_column(covmat,j),sigma));
      
      RETURN_ARRAYS_DECREMENT();
      return (covmat);
    }

  dvar_vector nll_logistnorm( dmatrix&  obs,  dvar_matrix& exp,dvector& wts, dvar_matrix& covmat,dvar_matrix& Vmat,dvar_matrix& Vinv,dvar_matrix& ww)
  {
    //calculate the negtive log likelihood of logistic-normal distribution
    RETURN_ARRAYS_INCREMENT();
  
      int nBin = obs.colsize();
      int nYear = obs.rowsize();
      dvar_vector negloglik(1,nYear);
      dmatrix O(1,nYear,1,nBin);
      dvar_matrix E(1,nYear,1,nBin);
      dmatrix Kmat(1,nBin-1,1,nBin);
      //dvar_matrix Vmat(1,nBin-1,1,nBin-1);
      //dvar_matrix Vinv(1,nBin-1,1,nBin-1);
      //dvar_matrix ww(1,nYear,1,nBin);
     
      for(int i=1,j=obs.colmin();i<=nBin;i++,j++){
        O.colfill(i,extract_column(obs,j));
        E.colfill(i,extract_column(exp,j));
       }
     
      for(int j=1;j<=nBin-1;j++){
      	ww.colfill(j,(log(extract_column(O,j))-log(extract_column(O,nBin)))-(log(extract_column(E,j))-log(extract_column(E,nBin))));
      }
      dmatrix temp_m = identity_matrix(1,nBin-1);
      dvector temp_v(1,nBin-1);
      temp_v = -1;
      for(int j=1;j<=nBin-1;j++) {
      	Kmat.colfill(j, extract_column(temp_m,j));
      }
      Kmat.colfill(nBin,temp_v);
      Vmat = Kmat * covmat * trans(Kmat);
      Vinv = inv(Vmat);

      negloglik = 0.5*(nBin-1)*log(2*pi)+rowsum(log(O))+0.5*log(det(Vmat))+(nBin-1)*log(wts);
      for(int i=1;i<=nYear;i++)
      	negloglik(i) = negloglik(i)+(0.5/(wts(i)*wts(i)))*ww(i) * Vinv * ww(i);
      
      RETURN_ARRAYS_DECREMENT();
      return (negloglik);
  } 
   
   

