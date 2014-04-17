  #define HOME_VERSION
  #include "fvar.hpp"
  #include "multinomial.h"
  
  void dfcholeski_solve(void);
  
  dvar_vector choleski_solve(_CONST dvar_matrix& MM,const dvar_vector& vv,
    const prevariable& det,const int& sgn)
  {
    // kludge to deal with constantness
    if (MM.colsize() != MM.rowsize())
    {
      cerr << "Error in chol_decomp. Matrix not square" << endl;
      ad_exit(1);
    }
    dmatrix M=value(MM);
    dvector v=value(vv);
    // int rowsave=M.rowmin();
    // int colsave=M.colmin();
    M.rowshift(1);
    M.colshift(1);
    int n=M.rowmax();
  
    dmatrix L(1,n,1,n);
  #ifndef SAFE_INITIALIZE
      L.initialize();
  #endif
  
    int i,j,k;
    double tmp;
      if (M(1,1)<=0)
      {
        cerr << "Error matrix not positive definite in choleski_decomp"
          <<endl;
        ad_exit(1);
      }
    L(1,1)=sqrt(M(1,1));
    for (i=2;i<=n;i++)
    {
      L(i,1)=M(i,1)/L(1,1);
    }
  
    for (i=2;i<=n;i++)
    {
      for (j=2;j<=i-1;j++)
      {
        tmp=M(i,j);
        for (k=1;k<=j-1;k++)
        {
          tmp-=L(i,k)*L(j,k);
        }
        L(i,j)=tmp/L(j,j);
      }
      tmp=M(i,i);
      for (k=1;k<=i-1;k++)
      {
        tmp-=L(i,k)*L(i,k);
      }
      if (tmp<=0)
      {
        cerr << "Error matrix not positive definite in choleski_decomp"
          <<endl;
        ad_exit(1);
      }
      L(i,i)=sqrt(tmp);
    }

    double cdet=0.0;
    for (int i=1;i<=n;i++)
    {
      cdet+=log(L(i,i));
    }
  
    v.shift(1);
    dvector x(1,n);
    x(1)=v(1)/L(1,1);
    for (int i=2;i<=n;i++)
    {
      double ssum=0.0;
      for (int j=1;j<=i-1;j++)
      {
        ssum+=L(i,j)*x(j);
      }
      x(i)=(v(i)-ssum)/L(i,i);
    }
  
    save_identifier_string("PO");
    value(det)=cdet;
    det.save_prevariable_position();
    save_identifier_string("SY");
    v.save_dvector_value();
    v.save_dvector_position();
    save_identifier_string("XY");
    dvar_vector vx=nograd_assign(x);
    save_identifier_string("TY");
    vx.save_dvar_vector_position();
    //save_identifier_string("rs");
    //vc.save_dvar_matrix_position();
    save_identifier_string("rt");
    MM.save_dvar_matrix_value();
    save_identifier_string("rl");
    MM.save_dvar_matrix_position();
    save_identifier_string("RO");
    vv.save_dvar_vector_position();
    save_identifier_string("QQ");
    gradient_structure::GRAD_STACK1->
        set_gradient_stack(dfcholeski_solve);
    return vx;
  }
  
  void dfcholeski_solve(void)
  {
    verify_identifier_string("QQ");
    dvar_vector_position vvpos=restore_dvar_vector_position();
    verify_identifier_string("RO");
    dvar_matrix_position MMpos=restore_dvar_matrix_position();
    verify_identifier_string("rl");
    dmatrix M=restore_dvar_matrix_value(MMpos);
    verify_identifier_string("rt");
    dvar_vector_position vxpos=restore_dvar_vector_position();
    verify_identifier_string("TY");
    dvector dfx=restore_dvar_vector_derivatives(vxpos);
    verify_identifier_string("XY");
    dvector_position vpos=restore_dvector_position();
    dvector v=restore_dvector_value(vpos);
    verify_identifier_string("SY");
    prevariable_position detpos=restore_prevariable_position();
    double dfdet=restore_prevariable_derivative(detpos);
    verify_identifier_string("PO");
  
    if (M.colsize() != M.rowsize())
    {
      cerr << "Error in chol_decomp. Matrix not square" << endl;
      ad_exit(1);
    }
    int rowsave=M.rowmin();
    int colsave=M.colmin();
    M.rowshift(1);
    M.colshift(1);
    int n=M.rowmax();
  
    dmatrix L(1,n,1,n);
    dmatrix dfL(1,n,1,n);
    dvector tmp(1,n);
    dmatrix tmp1(1,n,1,n);
    dmatrix dftmp1(1,n,1,n);
    dmatrix dfM(1,n,1,n);
    dvector dftmp(1,n);
    dvector dfv(1,n);
    tmp.initialize();
    tmp1.initialize();
    dftmp.initialize();
    dftmp1.initialize();
    dfM.initialize();
    dfL.initialize();
    dfv.initialize();
  #ifndef SAFE_INITIALIZE
      L.initialize();
  #endif
  
    int i,j,k;
    if (M(1,1)<=0)
    {
      cerr << "Error matrix not positive definite in choleski_decomp"
        <<endl;
      ad_exit(1);
    }
    L(1,1)=sqrt(M(1,1));
    for (i=2;i<=n;i++)
    {
      L(i,1)=M(i,1)/L(1,1);
    }
  
    for (i=2;i<=n;i++)
    {
      for (j=2;j<=i-1;j++)
      {
        tmp1(i,j)=M(i,j);
        for (k=1;k<=j-1;k++)
        {
          tmp1(i,j)-=L(i,k)*L(j,k);
        }
        L(i,j)=tmp1(i,j)/L(j,j);
      }
      tmp(i)=M(i,i);
      for (k=1;k<=i-1;k++)
      {
        tmp(i)-=L(i,k)*L(i,k);
      }
      if (tmp(i)<=0)
      {
        cerr << "Error matrix not positive definite in choleski_decomp"
          <<endl;
        ad_exit(1);
      }
      L(i,i)=sqrt(tmp(i));
    }

    double cdet=0.0;
    for (int i=1;i<=n;i++)
    {
      cdet+=log(L(i,i));
    }
  
    v.shift(1);
    dvector vsum(2,n);
    dvector x(1,n);
    x(1)=v(1)/L(1,1);
    for (int i=2;i<=n;i++)
    {
      double ssum=0.0;
      for (int j=1;j<=i-1;j++)
      {
        ssum+=L(i,j)*x(j);
      }
      vsum(i)=ssum;
      x(i)=(v(i)-vsum(i))/L(i,i);
    }
    dvector dfvsum(2,n);
    dfvsum.initialize();
  
    //for (int i=2;i<=n;i++)
    for (int i=n;i>=2;i--)
    {
      //x(i)=(v(i)-vsum(i))/L(i,i);
      dfv(i)+=dfx(i)/L(i,i);
      dfvsum(i)-=dfx(i)/L(i,i);
      dfL(i,i)-=dfx(i)*(v(i)-vsum(i))/square(L(i,i));
      dfx(i)=0.0;
      //for (int j=1;j<=i-1;j++)
      for (int j=i-1;j>=1;j--)
      {
        //vsum(i)+=L(i,j)*x(j);
        dfL(i,j)+=dfvsum(i)*x(j);
        dfx(j)+=dfvsum(i)*L(i,j);
      }
      dfvsum(i)=0;
    }
    //x(1)=v(1)/L(1,1);
    dfv(1)+=dfx(1)/L(1,1);
    dfL(1,1)-=dfx(1)*v(1)/square(L(1,1));
  
    
    for (int i=1;i<=n;i++)
    {
      //cdet+=log(L(i,i));
      dfL(i,i)+=dfdet/L(i,i);
    }
  
  
    for (i=n;i>=2;i--)
    {
      //L(i,i)=sqrt(tmp(i));
      dftmp(i)+=dfL(i,i)/(2.0*L(i,i));
      dfL(i,i)=0.0;
      for (k=i-1;k>=1;k--)
      {
        //tmp(i)-=L(i,k)*L(i,k);
        dfL(i,k)-=2.*dftmp(i)*L(i,k);
      }
      //tmp(i)=M(i,i);
      dfM(i,i)+=dftmp(i);
      dftmp(i)=0.0;
      for (j=i-1;j>=2;j--)
      {
        //L(i,j)=tmp1(i,j)/L(j,j);
        double linv=1./L(j,j);
        dftmp1(i,j)+=dfL(i,j)*linv;
        dfL(j,j)-=dfL(i,j)*tmp1(i,j)*linv*linv;
        dfL(i,j)=0.0;
        for (k=j-1;k>=1;k--)
        {
          //tmp(i,j)-=L(i,k)*L(j,k);
          dfL(i,k)-=dftmp1(i,j)*L(j,k);
          dfL(j,k)-=dftmp1(i,j)*L(i,k);
        }
        //tmp(i,j)=M(i,j);
        dfM(i,j)+=dftmp1(i,j);
        dftmp1(i,j)=0.0;
      }
    }
    double linv=1./L(1,1);
    for (i=n;i>=2;i--)
    {
      //L(i,1)=M(i,1)/L(1,1);
      dfM(i,1)+=dfL(i,1)*linv;
      dfL(1,1)-=dfL(i,1)*M(i,1)*linv*linv;
      dfL(i,1)=0.0;
    }
    //L(1,1)=sqrt(M(1,1));
    dfM(1,1)+=dfL(1,1)/(2.*L(1,1));
  
    dfM.rowshift(rowsave);
    dfM.colshift(colsave);
  
    save_double_derivative(dfdet,detpos);
    dfM.save_dmatrix_derivatives(MMpos);
    dfv.save_dvector_derivatives(vvpos);
  }
  
  #undef HOME_VERSION
