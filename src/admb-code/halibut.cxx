/** get_ft
  Solving the baranov catch equation using Newtons method
  -this implementation is for multiple fleets where the catch
   is based on weight.
*/
dmatrix get_ft_bysex(dvector& ct,const dvector& m, const d3_array& V,const dmatrix& ba)
{
	/*  ARGUMENTS:
	   ct is the observed catch for each fleet summed over sex
	   m is the natural mortality rate by sex
	   va is the vulnerability by sex, row fleet, col age va(1,ngear,1,nage)
	   ba is the start of year biomass at age by sex
	*/
	//cout<<"entering get_ft_bysex"<<endl;
	int h,i,a,A;
	double minsurv = 0.02;
	int ng=size_count(ct);	//number of gears
	int nh=size_count(m);	//number of sexs
	a=ba.colmin(); A=ba.colmax();
	dmatrix ft(1,nh,1,ng);
	dmatrix ctmp(1,ng,1,nh);
	dmatrix ptmp(1,ng,1,nh);			//proportion of predicted catch 
	dmatrix btmp(1,ng,1,nh);
	dmatrix ct_hat(1,nh,1,ng);	//predicted catch
	dmatrix dct_hat(1,nh,1,ng);	//derivative
	dmatrix zage(1,nh,a,A);
	dmatrix sage(1,nh,a,A);
	dmatrix ominus(1,nh,a,A);
	d3_array F(1,nh,1,ng,a,A);
	
	
	
	//initial guess for ft=0.5ct/(0.98 Bt); assuming 50:50 female:male catch by weight.
	//ctmp=ct;
	
	for(i=1;i<=ng;i++)
	{
		for(h=1;h<=nh;h++)
		{
			btmp(i,h) = ba(h)*V(h)(i)*exp(-0.5*m(h));
		}
		ptmp(i)  = btmp(i) / sum(btmp(i));
		ctmp(i)  = ptmp(i) *ct(i);
	}
	//cout<<"ctmp\n"<<ctmp<<endl;
	
	for(h=1;h<=nh;h++)
	{
		for(i=1;i<=ng;i++)
		{   
			ft(h,i) = ctmp(i,h)/(0.98*btmp(i,h));

			if(1.-ft(h)(i)<minsurv)
			{
				//cout<<"minsurv"<<endl;
				ft(h)(i)=1.-minsurv;
				ctmp(i,h)=ft(h)(i)*0.98*btmp(i,h);
			}
		}
	}
	//ct=ctmp;	//don't do this for the differentiable version.
	
	
	// March 26, Still working here trying to solve the catch equation for sex and gears.
	// Need to calculate total catch by sex to get newton step right.
	
	//now solve baranov catch equation iteratively.
	for(int iter=1; iter<=100; iter++)
	{ 
		
		for(h=1;h<=nh;h++)
		{
			for(i=1;i<=ng;i++)
			{
				//cout<<ft(h)(i)<<endl;
				F(h)(i)=ft(h)(i)*V(h)(i);
			}
			
			zage(h)=m(h)+colsum(F(h));
			sage(h)=mfexp(-zage(h));
			ominus(h)=(1.-sage(h));
		}
		
		for(i=1;i<=ng;i++)
		{   
			for(h=1;h<=nh;h++)
			{
				dvector t1 = elem_div(F(h)(i),zage(h));
				dvector t2 = elem_prod(t1,ominus(h));
				dvector t3 = elem_prod(ominus(h),ba(h));
			
				ct_hat(h)(i) = t2*ba(h);
				
				dct_hat(h)(i)= sum(
							elem_div(t3,zage(h))
							- elem_div(elem_prod(F(h)(i),t3),square(zage(h)))
							+ elem_prod(elem_prod(t1,sage(h)),ba(h)));
							
				ft(h)(i) -= (ct_hat(h)(i)-ctmp(i,h))/dct_hat(h)(i);
			}
		}
		//cout<<iter<<"\t"<<colsum(ct_hat)-ct<<endl;//"\t"<<ptmp<<endl;
	}
	//cout<<ft<<"\t\t"<<ct<<"\t\t"<<ctmp<<endl;
	
	return(ft);
	
}  

