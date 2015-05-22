// milka.cpp
/**
 * @Milka Source code for Operating Model
 * @author Steven Martell & Catarina Wor
 * @details The default constuctort uses model_data as the base
 * class for the OperatingModel class.  
 * 
 * The OperatingModel class has a major function called 
 * runScenario:
 * 
 * ———————————————————————————————————————————————— *
 * STATUS LEGEND
 *    : not implemented yet.
 *  - : partially implemented
 *  + : implemented & testing
 *   : Good to go! 
 * ———————————————————————————————————————————————— * 
 * runScenario:                                STATUS
 *      |- readMSEcontrols                     [-]
 *      |- initParameters                      [-]
 *          |- surveyQ                         [ ]
 *          |- stock-recruitment parameters    [ ]
 *      |- initMemberVariables                 [-]
 *      |- conditionReferenceModel             [-]
 *      |- setRandomVariables                  [-]
 *      |- | getReferencePointsAndStockStatus  [-]
 *         | calculateTAC                      [-]
 *         | allocateTAC                       [-]
 *         | implementFisheries                [-]
 *              |- calcSelectivity             [ ]
 *              |- calcRetentionDiscards       [ ]
 *              |- calcTotalMortality          [-]
 *         | calcRelativeAbundance             [-]
 *         | calcCompositionData               [-]
 *         | calcEmpiricalWeightAtAge          [-]
 *         | updateReferenceModel              [-]
 *         | writeDataFile                     [-]
 *         | runStockAssessment                [-]
 *      |- |            
 *      |- writeSimulationVariables            [-]
 *      |- calculatePerformanceMetrics         [ ]
 * ———————————————————————————————————————————————— *
 */

#include <admodel.h>
#include "milka.h"
#include "include/lib_iscam.h"

#undef COUT
#define COUT(object) cout<<#object"\n"<<object<<endl;
// Destructor
OperatingModel::~OperatingModel(){}

// Constructor
OperatingModel::OperatingModel(ModelVariables _mv,int argc,char * argv[])
:model_data(argc,argv), mv(_mv)
{
    // cout<<"Inheritance version using model_data as base class"<<endl;
    // cout<<"Ngroup "<<ngroup<<endl;
    // cout<<"Catch Data\n"<<dCatchData<<endl;
    // cout<<"d3 Survey Data\n"<<d3_survey_data<<endl;
    // cout<<"eof "<<eof<<endl;

}

/**
 * @brief Verify Equilibrium MSY calculations numerically.
 * @details This runs the model for upto 100 years to verify
 * that the equilibrium MSY-based reference points are calcualted correctly.
 * 
 */
void OperatingModel::checkMSYcalcs()
{
    dvector tmp_tau = m_dTau;
    m_dTau = 0;


    readMSEcontrols();

    initParameters();

    initMemberVariables();

    conditionReferenceModel();

    calcMSY();

    for(int i = nyr+1; i <= m_nPyr; i++ )
    {
        calculateTAC();
        
        calcTotalMortality(i);

        implementFisheries(i);

        updateReferenceModel(i);
    }
    // COUT(m_sbt);
    // COUT(m_d3_Ct);


    m_dTau = tmp_tau;
}

void OperatingModel::runScenario(const int &seed)
{
    readMSEcontrols();

    initParameters();

    initMemberVariables();

    conditionReferenceModel();

    calcMSY();

    setRandomVariables(seed);

    for(int i = nyr+1; i <= m_nPyr; i++ )
    {
        getReferencePointsAndStockStatus(i);
        if(verbose) cout<<"getReferencePointsAndStockStatus OK"<<endl;

        calculateTAC();
        if(verbose) cout<<"calculateTAC OK"<<endl;

        allocateTAC(i);
        if(verbose) cout<<"allocateTAC OK"<<endl;

        implementFisheries(i);
        if(verbose) cout<<"implementFisheries OK"<<endl;

        calcTotalMortality(i);
        if(verbose) cout<<"calcTotalMortality OK"<<endl;

        calcRelativeAbundance(i);
        if(verbose) cout<<"calcRelativeAbundance OK"<<endl;

        calcCompositionData(i);
        if(verbose) cout<<"calcCompositionData OK"<<endl;

        calcEmpiricalWeightAtAge(i);
        if(verbose) cout<<"calcEmpiricalWeightAtAge OK"<<endl;

        updateReferenceModel(i);
        if(verbose) cout<<"updateReferenceModel OK"<<endl;

        writeDataFile(i);
        if(verbose) cout<<"writeDataFile OK"<<endl;

        writeParameterFile(i);
        if(verbose) cout<<"writeParameterFile OK"<<endl;

        // implement perfect info option
        // if flag is 0 - write .res with true params
        runStockAssessment();
    
    }

    writeSimulationVariables();
}


/**
 * @brief Calculate MSY-based reference points for reference model.
 * @details This routine calculates the MSY-based reference points based 
 * on the true parameter values and selectivity coefficients used in 
 * the population reference model.  Needed for perfect information 
 * scenarios.
 * 
 * Requires msy.hpp
 * rfp::msy(ro,h,rho,ma,wa,fa,_V)
 * 
 * rho is the fraction of mortality that occurs prior to spawning.
 */
void OperatingModel::calcMSY()
{
    m_msy.allocate(1,ngroup,1,nfleet);
    m_bmsy.allocate(1,ngroup);
    m_fmsy.allocate(1,ngroup,1,nfleet);
    m_msy.initialize();
    m_fmsy.initialize();
    m_bmsy.initialize();


    /* Fecundity at age and natural mortality*/
    dmatrix fa_bar(1,n_ags,sage,nage);
    dmatrix  M_bar(1,n_ags,sage,nage);
    for(int ig=1;ig<=n_ags;ig++)
    {
        fa_bar(ig) = elem_prod(dWt_bar(ig),ma(ig));
        M_bar(ig)  = colsum(value(m_M(ig).sub(pf_cntrl(3),pf_cntrl(4))));
        M_bar(ig) /= pf_cntrl(4)-pf_cntrl(3)+1; 
    }

    // | (1) : Matrix of selectivities for directed fisheries.
    // |     : log_sel(gear)(n_ags)(year)(age)
    // |     : ensure dAllocation sums to 1.
    dvector d_ak(1,nfleet);
    d3_array  d_V(1,n_ags,1,nfleet,sage,nage);
    
    for(int k=1;k<=nfleet;k++)
    {
        int kk  = nFleetIndex(k);
        d_ak(k) = dAllocation(kk);
        for(int ig=1;ig<=n_ags;ig++)
        {
            d_V(ig)(k) = ( exp(d4_logSel(kk)(ig)(nyr)) );
        }
    }
    d_ak /= sum(d_ak);

    
    // initial guess for fmsy
    dvector dftry(1,nfleet);
    dftry  = 0.6/nfleet * mean(M_bar);
    
    
    double rho = d_iscamCntrl(13);
    for(int g=1; g<=ngroup; g++)
    {
        rfp::msy<double,dvector,dmatrix,d3_array>
        cMSY(m_dRo(g),m_dSteepness(g),rho,M_bar,dWt_bar,fa_bar,d_V);
        
        m_fmsy(g) = cMSY.getFmsy(dftry,d_ak);
        m_bmsy(g) = cMSY.getBmsy();
        m_msy(g)  = cMSY.getMsy();
        cMSY.print();
    }
}


/**
 * @brief Read control file for Management Strategy Evaluation.
 * @details Use cifstream to read in controls for MSE related options.
 * 
 */
void OperatingModel::readMSEcontrols()
{
    if(verbose) cout<<"MSE Control file\n"<<ProcedureControlFile<<endl;
    if(verbose) cout<<"MSE Scenario file\n"<<ScenarioControlFile<<endl;

    cifstream ifs_mpc(ProcedureControlFile);
    ifs_mpc >> m_nPyr;

    //assessment option
    ifs_mpc >> m_nAssessOpt;

    // Control file.
    ifs_mpc >> m_controlFile;
    adstring flg = "none";
    if(m_controlFile == flg)
    {
        m_controlFile = ControlFile;
    }
    //harvest control rule
    ifs_mpc >> m_nHCR;

    ifs_mpc >> m_dBthreshold;
    ifs_mpc >> m_dBlimit;
    ifs_mpc >> m_maxf;

    m_nGearIndex.allocate(1,ngear);
    m_nCSex.allocate(1,ngear);
    m_nASex.allocate(1,ngear);
    m_nAGopen.allocate(1,ngear,1,narea);
    
    // Controls for sexing catch and comps and fishing in given areas.
    dmatrix tmp(1,ngear,-7,narea);
    ifs_mpc >> tmp;
    m_nGearIndex = ivector(column(tmp,-7));
    m_nCSex      = ivector(column(tmp,-6));
    m_nASex      = ivector(column(tmp,-5));
    m_nATau      = column(tmp,-4);
    m_nWSex      = ivector(column(tmp,-3));
    m_dLslim     = column(tmp,-2);
    m_dUslim     = column(tmp,-1);
    m_dDiscMortRate = column(tmp,0);

    for( k = 1; k <= ngear; k++ )
    {
        m_nAGopen(k) = ivector(tmp(k)(1,narea));
    }


    int eof=0;
    ifs_mpc >> eof;
    // cout<<"End of MPC file "<<eof<<endl;
    if(eof != 999)
    {
        cout<<"Error reading Management Procedure Control File"<<endl;
        cout<<eof<<endl;
        ad_exit(1);
    }
    //cout<<"finished MSE controls"<<endl;






    //
    // READ SCENARIO CONTROL FILE
    // 
    cifstream ifs_scn(ScenarioControlFile);
    
    // 
    // Controls for recruitment options
    // 
    ifs_scn >> m_nRecType;

    // 
    // Dispersel kernel for new recruits.
    // 
    m_dispersal.allocate(1,narea,1,narea); m_dispersal.initialize();
    ifs_scn >> m_dispersal; 

    // 
    // Autocorrelation coefficient in recruitment.
    // 
    ifs_scn >> m_gamma_r;

    // 
    // Recruitment Regime.
    // 
    ifs_scn >> m_PDO_phase;

    // 
    // Size-at-age (increase or decrease)
    // 
    ifs_scn >> m_SAA_flag;
    


    // End of file
    int eof_scn=0;
    ifs_scn >> eof_scn;
    if(eof_scn != 999)
    {
        cout<<"Error reading Scenario file"<<endl;
        cout<<eof_scn<<endl;
        ad_exit(1);
    }

}

/**
 * @brief Initialize model parameters based on model variable struct.
 * @details [long description]
 */
void OperatingModel::initParameters()
{

    // Initializing data members
    m_nNyr = nyr; // needs to be updated for each year inside the mse loop do we need this here??
    m_irow = nCtNobs; // counter for current number of rows in the catch table.
    m_nyrs = m_nPyr - m_nNyr;
    
    m_yr.allocate(1,m_nPyr-syr+1);
    m_yr.fill_seqadd(syr,1);

    // needs to be updated for each year in the mse loop

    // m_nn is a counter for the number of rows of catch data that will be
    // added to the data file each year.
    m_nn = 0;
    // for( k = 1; k <= ngear; k++ )
    for( k = 1; k <= nfleet; k++ )
    {
        m_nn += sum(m_nAGopen(k));
        m_nn += m_nCSex(k)*m_nn;
    }
    m_nCtNobs = nCtNobs + m_nyrs*m_nn;
    
    int ncol = dCatchData.colmax();
    m_dCatchData.allocate(1,m_nCtNobs,1,ncol);
    m_dCatchData.initialize();
    m_dCatchData.sub(1,nCtNobs) = dCatchData;
    

    m_dSubLegalData.allocate(nCtNobs+1,m_nCtNobs,1,8);
    m_dSubLegalData.initialize();

     // Fishing mortality rate parameters
    // How many more Ft parameters are needed?
    m_log_ft_pars.allocate(1,m_nCtNobs);
    m_log_ft_pars.initialize();
    m_log_ft_pars.sub(1,ft_count) = mv.log_ft_pars(1,ft_count);

    // Below could be deprecated
   // m_d3_Ct.allocate(1,n_ags,syr,m_nPyr,1,ngear);
   // m_d3_Ct.initialize();
    //for(int ig = 1; ig <= n_ags; ig++ )
    //{
    //    m_d3_Ct(ig).sub(syr,nyr) = d3_Ct(ig);
    //}
    

    // allocate & initialize survey data arrays
    m_n_it_nobs.allocate(1,nItNobs);
    m_n_it_nobs.initialize();
    m_n_it_nobs = n_it_nobs + m_nyrs;
    
    m_d3SurveyData.allocate(1,nItNobs,1,m_n_it_nobs,1,9);
    m_d3SurveyData.initialize();
    for(int k=1;k<=nItNobs;k++)
    {
        m_d3SurveyData(k).sub(1,n_it_nobs(k)) = d3_survey_data(k);  
    }
    // Age-composition arrays   
    m_A_irow.allocate(1,nAgears);
    m_A_irow.initialize(); 

    m_n_A_nobs.allocate(1,nAgears);
    m_n_A_nobs.initialize();
    for( k = 1; k <= nAgears; k++ )
    {
        if(n_A_nobs(k) > 0)
        {
            m_n_A_nobs(k) = n_A_nobs(k) + m_nyrs + m_nyrs * m_nASex(k);
        }
    }
    
    m_d3_A.allocate(1,nAgears,1,m_n_A_nobs,n_A_sage-6,n_A_nage);
    m_d3_A.initialize();
    // cout<<"Ok today"<<endl;
    
    for(int k=1;k<=nAgears;k++)
    {
        m_d3_A(k).sub(1,n_A_nobs(k)) = d3_A(k); 
    }
        
    // Weight-at-age array
    m_W_irow.allocate(1,nWtTab);
    m_W_irow.initialize(); 

    m_nWtNobs.allocate(1,nWtTab);
    m_nWtNobs.initialize();
    for( k = 1; k <= nWtTab; k++ )
    {
        if( nWtNobs(k) > 0 )
        {
            m_nWtNobs(k) = nWtNobs(k) + m_nyrs + m_nyrs * sum(m_nWSex);
        }
    }

    m_d3_inp_wt_avg.allocate(1,nWtTab,1,m_nWtNobs,sage-5,nage);
    m_d3_inp_wt_avg.initialize();
    for(int k=1;k<=nWtTab;k++)
    {
        m_d3_inp_wt_avg(k).sub(1,nWtNobs(k)) = d3_inp_wt_avg(k);
    }
    
    // TODO: allow user to specify observation error
    // initializing population parameters
    m_dRo        = exp(mv.log_ro);
    m_dBo        = mv.sbo;
    m_dSteepness = mv.steepness;
    m_dM         = exp(mv.m);
    m_dRho       = mv.rho;    // now autocorrelation
    m_dVarphi    = mv.varphi; //sqrt(1.0/mv.varphi);
    m_dSigma     = 0.2 * mv.varphi;//elem_prod(sqrt(m_dRho) , m_dVarphi);
    m_dTau       = m_dVarphi; //elem_prod(sqrt(1.0-m_dRho) , m_dVarphi);

    m_dRbar.allocate(1,n_ag);
    m_dRinit.allocate(1,n_ag);
    for(int ih = 1; ih <= n_ag; ih++ )
    {
        m_dRbar  = exp(mv.log_rbar(ih));
        m_dRinit = exp(mv.log_rinit(ih));
    }

    switch(int(d_iscamCntrl(2)))
    {
        case 1:
            //Beverton-Holt model
            m_dKappa = elem_div(4.*m_dSteepness,(1.-m_dSteepness));
            break;
        case 2:
            //Ricker model
            m_dKappa = pow((5.*m_dSteepness),1.25);
        break;
    }


    
    // cout<<"finished init parameters"<<endl;
}


void OperatingModel::initMemberVariables()
{
    m_N.allocate(1,n_ags,syr,m_nPyr+1,sage,nage); m_N.initialize();
    m_M.allocate(1,n_ags,syr,m_nPyr,sage,nage); m_M.initialize();
    m_F.allocate(1,n_ags,syr,m_nPyr,sage,nage); m_F.initialize();
    m_Z.allocate(1,n_ags,syr,m_nPyr,sage,nage); m_Z.initialize();
    m_S.allocate(1,n_ags,syr,m_nPyr,sage,nage); m_S.initialize();
    m_ft.allocate(1,n_ags,1,ngear,syr,m_nPyr);  m_ft.initialize();
    m_d3_wt_avg.allocate(1,n_ags,syr,m_nPyr+1,sage,nage); m_d3_wt_avg.initialize();
    m_d3_wt_mat.allocate(1,n_ags,syr,m_nPyr+1,sage,nage); m_d3_wt_mat.initialize();

    m_log_rt.allocate(1,n_ag,syr-nage+sage,nyr); m_log_rt.initialize();
    
    m_est_bo.allocate(1,ngroup);
    m_est_bmsy.allocate(1,ngroup);
    m_est_sbtt.allocate(1,ngroup);
    m_est_btt.allocate(1,ngroup);
    m_est_fmsy.allocate(1,ngroup,1,nfleet);
    m_est_msy.allocate(1,ngroup,1,nfleet);
    m_est_N.allocate(1,n_ags,sage,nage);
    m_est_M.allocate(1,n_ags,sage,nage);
    m_est_wa.allocate(1,n_ags,sage,nage);
    m_est_log_sel.allocate(1,n_ags,sage,nage);

    //Spawning stock biomass
    m_sbt.allocate(syr,m_nPyr,1,ngroup);m_sbt.initialize();
    m_sbt.sub(syr,nyr)=(trans(mv.sbt)).sub(syr,nyr);
    //total biomass
    m_bt.allocate(syr,m_nPyr,1,ngroup);m_bt.initialize();
    m_bt.sub(syr,nyr)=(trans(mv.bt)).sub(syr,nyr);

    m_dbeta.allocate(1,ngroup);m_dbeta.initialize();

    m_dTAC.allocate(1,ngroup,1,nfleet);

    m_q = mv.q;

    // Initialize Mortality arrays from ModelVariables (mv)
    // cohort-specific weight-at-age deviate
    dvector wa_dev(nyr-sage,m_nPyr);
    wa_dev.initialize();
    wa_dev = 0;
    switch( m_SAA_flag )
    {
        case 1:
            wa_dev = 0.3;
        break;
        case -1:
            wa_dev = -0.3;
        break;
        default:
            wa_dev = 0;
        break;
    }

    for(int ig = 1; ig <= n_ags; ig++ )
    {
        m_M(ig).sub(syr,nyr) = (*mv.d3_M)(ig);
        m_F(ig).sub(syr,nyr) = (*mv.d3_F)(ig);
        m_Z(ig).sub(syr,nyr) = m_M(ig).sub(syr,nyr) + m_F(ig).sub(syr,nyr);
        m_S(ig).sub(syr,nyr) = exp(-m_Z(ig).sub(syr,nyr));
        m_d3_wt_avg(ig).sub(syr,nyr+1) = d3_wt_avg(ig).sub(syr,nyr+1);
        m_d3_wt_mat(ig).sub(syr,nyr+1) = d3_wt_mat(ig).sub(syr,nyr+1);

        // Temporary extend natural mortality out to m_nPyr
        // Modify m_d3_wt_avg & m_d3_wt_mat to accomodate changes
        // in size-at-age in the future.  The idea would be that
        // each cohort would be given a multiplier specific to that 
        // cohort. 
        // m_SAA_flag = 1, size at age increases
        // m_SAA_flag =-1, size-at-age decreases
        for( i = nyr+1; i <= m_nPyr; i++ )
        {
            m_M(ig)(i) = m_M(ig)(nyr);
            m_d3_wt_avg(ig)(i+1) = d3_wt_avg(ig)(nyr+1);
            m_d3_wt_mat(ig)(i+1) = d3_wt_mat(ig)(nyr+1);

            for( j = sage; j <= nage; j++)
            {
                int cohort = i - j;
                if(cohort <= nyr+1) continue;
                
                m_d3_wt_avg(ig)(i+1)(j) = d3_wt_avg(ig)(nyr+1)(j)*exp(wa_dev(cohort)); 
                m_d3_wt_mat(ig)(i+1)(j) = d3_wt_mat(ig)(nyr+1)(j)*exp(wa_dev(cohort)); 
            }
        }
    }
    // COUT(m_d3_wt_avg(1));
    // exit(1);
    // Selectivity
    d4_logSel.allocate(1,ngear,1,n_ags,syr,m_nPyr,sage,nage);
    d4_logSel.initialize();
    for( k = 1; k <= ngear; k++ )
    {
        for(int ig = 1; ig <= n_ags; ig++ )
        {
            d4_logSel(k)(ig).sub(syr,nyr) = (*mv.d4_logSel)(k)(ig);

            // Temporarily extend selectivity out to m_nPyr
            for( i = nyr+1; i <= m_nPyr; i++ )
            {
                d4_logSel(k)(ig)(i) = d4_logSel(k)(ig)(nyr);
            }
        }
    }

    // annual fishing mortality rates
    m_ft.allocate(1,n_ags,1,ngear,syr,m_nPyr);
    for(int ig = 1; ig <= n_ags; ig++ )
    {
        for(int k = 1; k <= ngear; k++ )
        {
            m_ft(ig)(k)(syr,nyr) = (*mv.d3_ft)(ig)(k);
             /* code */
        }
    }


    //cout<<"finished init member variables"<<endl;

}

void OperatingModel::conditionReferenceModel()
{
    int ig,ih;

    for( ig = 1; ig <= n_ags; ig++ )
    {
        f  = n_area(ig);
        g  = n_group(ig);
        ih = pntr_ag(f,g);

        dvector lx(sage,nage);
        dvector tr(sage,nage);
        lx(sage) = 1.0;
        for(j=sage;j< nage;j++)
        {
            lx(j+1) = lx(j) * exp( -m_M(ig)(syr)(j) );
        }
        lx(nage) /= (1.-exp(-m_M(ig)(syr,nage)));
        
        if( d_iscamCntrl(5) ) // initialize at unfished conditions.
        {
            tr =  log( m_dRo(g) ) + log(lx);
        }
        else if ( !d_iscamCntrl(5) )
        {
            tr(sage)        = ( mv.log_rbar(ih)+mv.log_rec_devs(ih)(syr));
            tr(sage+1,nage) = (mv.log_rinit(ih)+mv.init_log_rec_devs(ih));
            tr(sage+1,nage) = tr(sage+1,nage)+log(lx(sage+1,nage));
        }
        m_N(ig)(syr)(sage,nage) = 1./nsex * mfexp(tr);
        m_log_rt(ih)(syr-nage+sage,syr) = tr.shift(syr-nage+sage);

        for(i=syr;i<=nyr;i++)
        {
            if( i>syr )
            {
                m_log_rt(ih)(i) = (mv.log_rbar(ih)+mv.log_rec_devs(ih)(i));
                m_N(ig)(i,sage) = 1./nsex * mfexp( m_log_rt(ih)(i) );               
            }

            m_N(ig)(i+1)(sage+1,nage) =++elem_prod(m_N(ig)(i)(sage,nage-1)
                                                 ,m_S(ig)(i)(sage,nage-1));
            m_N(ig)(i+1,nage)        +=  m_N(ig)(i,nage)*m_S(ig)(i,nage);

            // average biomass for group in year i
            //bt(g)(i) += N(ig)(i) * d3_wt_avg(ig)(i);
        }
        m_N(ig)(nyr+1,sage) = 1./nsex * mfexp( mv.log_rbar(ih));
    }

}

void OperatingModel::setRandomVariables(const int& seed)
{
    m_nSeed = seed;
    random_number_generator rng(m_nSeed);

    m_epsilon.allocate(1,nItNobs,nyr+1,m_nPyr);
    m_epsilon.fill_randn(rng);

    m_delta.allocate(1,ngroup,nyr-sage,m_nPyr);
    m_delta.fill_randn(rng);


    // Add autocorrelation to recruitment deviations
    double rho = m_gamma_r;
    for( g = 1; g <= ngroup; g++ )
    {
        m_delta(g) += m_PDO_phase;
        for( i = nyr-sage+1; i <= m_nPyr; i++ )
        {
            m_delta(g)(i) = rho * m_delta(g)(i-1) + sqrt(1.0-square(rho)) * m_delta(g)(i);
        }
    }
}


void OperatingModel::getReferencePointsAndStockStatus(const int& iyr)
{
    switch( int(m_nAssessOpt) ) // option read in from .mpc file
    {
        case 0:
            //  set reference points to true milka values
            
            m_est_bo   = m_dBo;
            m_est_fmsy = m_fmsy;      
            m_est_msy  = m_msy;       
            m_est_bmsy = m_bmsy;      
            m_est_sbtt = m_sbt(iyr)(1,ngroup);
            m_est_btt  = m_bt(iyr)(1,ngroup);;
            
            for(int ig = 1; ig <= n_ags; ig++ )
            {
                m_est_N(ig)(sage,nage) = m_N(ig)(iyr)(sage,nage);
            }
            
            for(int ig = 1; ig <= n_ags; ig++ )
            {
                m_est_wa(ig)(sage,nage) = m_d3_wt_avg(ig)(iyr)(sage,nage);
            }       

            for(int ig = 1; ig <= n_ags; ig++ )
            {
                m_est_M(ig)(sage,nage) = m_M(ig)(iyr)(sage,nage);
            }       
            //cout<<"TIme to goto dance"<<endl;
            //cout<<"Fmsy = "<<m_est_fmsy<<endl;


            // 4darray log_sel(1,ngear,1,n_ags,syr,nyr,sage,nage);
            for(int k = 1; k <= ngear; k++ )    
            {
                for(int ig = 1; ig <= n_ags; ig++ )
                {
                    m_est_log_sel(ig)(sage,nage)= d4_logSel(k)(ig)(iyr)(sage,nage);

                }
            }
            //exit(1);
        break;

        case 1:
            // read iscam.res file to get this information.
            cifstream ifs("iSCAM.res");
            ifs >> m_est_bo;
            ifs >> m_est_fmsy;
            ifs >> m_est_msy;
            ifs >> m_est_bmsy;
            ifs >> m_est_sbtt;
            ifs >> m_est_btt;
            ifs >> m_est_N;
            ifs >> m_est_wa;
            ifs >> m_est_M;
            ifs >> m_est_log_sel;
        break;
    }

}

/**
 * @brief Calculate the Total Allowable Catch
 * @details Total Allowable Catch is based on the estimates of current
 * stock statuts, reference points and the harvest control rule.
 * 
 * Uses a switch statement for HCR.  The HCR is set in the pcf file.
 */
void OperatingModel::calculateTAC()
{
    double btmp;
    double sbt;
    double sbo;
    dvector f_rate(1,nfleet);
    m_dTAC.initialize();


    for(int g = 1; g <= ngroup; g++ )
    {
        btmp = m_est_btt(g);
        sbt  = m_est_sbtt(g);
        sbo  = m_est_bo(g);
        for(int k = 1; k <= nfleet; k++ )
        {
            if(m_est_fmsy(g,k) > m_maxf)
            {
                m_est_fmsy(g,k) = m_maxf;
            }
        }
        
        


        switch( int(m_nHCR) )
        {
            case 1: // Constant harvest rate
                 // m_dTAC(g)  = harvest_rate * btmp;
                f_rate = m_est_fmsy(g);
            break; 

            case 2: // Bthreshold:Blimit HCR.
                double status = sbt/sbo;
                f_rate = m_est_fmsy(g);
                if( status < m_dBthreshold && status >= m_dBlimit )
                {
                    f_rate *= (status-m_dBlimit)/(m_dBthreshold-m_dBlimit);
                }
                else if(status < m_dBlimit)
                {
                    f_rate = 0;
                }
                // m_dTAC(g)  = harvest_rate * btmp;
                // cout<<"Status "<<status<<endl;
                // cout<<harvest_rate<<endl;
                // exit(1);
            break;
        }
    }
    
    dvector ba(sage,nage);
    dvector va(sage,nage);
    dvector ca(sage,nage);
    dvector za(sage,nage);

    // Todo: Check this routine below, not sure if it will work for multi-sex,multiarea multifleet.
    // Working here, need to implement the Baranov
    // Catch equation to calculate the m_dTAC for group g.
    // cout<<"I'm Here "<<endl;
    ba.initialize();
    m_dTAC.initialize();
    for(int ig = 1; ig <= n_ags; ig++ )
    {
        //int f = n_area(ig);
        int g = n_group(ig);
        //int h = n_sex(ig);

        va  = exp(m_est_log_sel(ig));
        // cout<<"And made it to here"<<endl;
        ba  = elem_prod(m_est_N(ig),m_est_wa(ig));

        for( k = 1; k <= nfleet; k++ )
        {
            za  = m_est_M(ig) + f_rate(k)*va;
            ca  = elem_prod(elem_prod(elem_div(f_rate(k)*va,za),1.0-exp(-za)),ba);
        }

        m_dTAC(g) += sum(ca);
    }

    // cout<<m_nHCR<<endl;
    // exit(1);
}


void OperatingModel::allocateTAC(const int& iyr)
{
    static int irow = nCtNobs;
    //m_dCatchdata(year,gear,area,group,sex,type,value)
    int h;
    for( k = 1; k <= nfleet; k++ )
    {
        h = m_nCSex(k);
        for( f = 1; f <= narea; f++ )
        {
            if(m_nAGopen(k,f))
            { 
            for( g = 1; g <= ngroup; g++ )
            {
                if(!h)
                {
                    irow ++;
                    m_dCatchData(irow,1) = iyr;
                    m_dCatchData(irow,2) = nFleetIndex(k);
                    m_dCatchData(irow,3) = f;
                    m_dCatchData(irow,4) = g;
                    m_dCatchData(irow,5) = h;
                    m_dCatchData(irow,6) = 1;             //TODO: Fix this catch type
                    m_dCatchData(irow,7) = m_dTAC(g)(k);  // TODO: call a manager!
                }
                if(h)
                {   
                    for( h = 1; h <= nsex; h++ )
                    {
                        irow ++;
                        m_dCatchData(irow,1) = iyr;
                        m_dCatchData(irow,2) = nFleetIndex(k);
                        m_dCatchData(irow,3) = f;
                        m_dCatchData(irow,4) = g;
                        m_dCatchData(irow,5) = h;
                        m_dCatchData(irow,6) = 1;             //TODO: Fix this
                        m_dCatchData(irow,7) = m_dTAC(g)(k);  // TODO: call a manager!
                    }
                }
            }
            }
        }
    }

    
}
    
/**
 * @brief Implement spatially explicity fishery.
 * @details Implement the spatially epxlicity fishery using the Baranov catch equation
 * to determine the instantaneous fishing mortality rate in each area by each gear. This
 * routine uses the BaranovCatchEquation class object to do this.
 * 
 * Notes:
 *  m_dTAC is a vector of allocated catches assiged to each fleet.
 *  
 *  Algorithm:
 *  |- Apportion m_dTAC by area (f) for each stock (g)
 *  |- Loop over each area and allocate catch in area (f) to gear (k),
 *  |- Add implementation error to each gear and catch.
 *  |- Assemble arguments for BaranovCatchEquation Class. 
 *       -> .getFishingMortality(ct,ma,&Va,na,wa,_hCt)
 *  |- Calculate Fishing mortality rates on reference population.
 *  |- Calculate Total landed catch.
 *  |- Calculate total discards based on size-limits.
 *  |- Calculate total discards from non-retention fisheries.
 *  
 *  
 *  NOTES on Joint probability of capture & retention.
 *  Defs:
 *      - Pc = probability of capture
 *      - Pr = probability of retention
 *      - Pd = probability of discarding (1-Pr).
 *      - dm = discard mortality rate.
 *      
 *  Joint probability model:
 *   Defs: Probability of retaining a fish of a given age a is:
 *   Va = Pc*(Pr + Pd*dm) = Pc(Pr+(1-Pr)*dm)
 *  
 *  The probability of retaining a fish is a function of its length
 *  and the variance in length-at-age.  To to this we assume that length
 *  at age is normaly distributed and the cumulative distibution is 
 *  defined by the cumd_norm(z) function, where z is the 
 *  (size_limit-mu)/sd;  This function is defined as the 
 *  retention_probabilty
 *  
 */
void OperatingModel::implementFisheries(const int &iyr)
{
    dvector tac(1,narea);
    dvector  ct(1,nfleet);
    dvector  pr(sage,nage);  // probability of retention
    dvector  pd(sage,nage);  // probability of discarding
    dmatrix  ma(1,nsex,sage,nage);
    dmatrix  na(1,nsex,sage,nage);
    dmatrix  wa(1,nsex,sage,nage);
    dmatrix  mu(1,nsex,sage,nage);
    dmatrix  sd(1,nsex,sage,nage);
    dmatrix  d_allocation(1,narea,1,nfleet);
    dmatrix  _hCt(1,nsex,1,nfleet);
    dmatrix  _hDt(1,nsex,1,nfleet);             // Discarded catch
    dmatrix  _hWt(1,nsex,1,nfleet);             // Wastage = (discard mort)*(discard catch)
    d3_array d3_Va(1,nsex,1,nfleet,sage,nage);  // Retained fraction
    d3_array d3_Da(1,nsex,1,nfleet,sage,nage);  // Discard fraction
    tac.initialize();
    na.initialize();
    _hDt.initialize();
    _hWt.initialize();
    static int ft_counter = ft_count;

    BaranovCatchEquation cBCE;

    for(int f = 1; f <= narea; f++ )
    {
        for(int g = 1; g <= ngroup; g++ )
        {
            ct = m_dTAC(g);  // Catch for each fleet.
            for(int h = 1; h <= nsex; h++ )
            {
                int ig = pntr_ags(f,g,h);
                ma(h) = m_M(ig)(iyr);           // natural mortality
                na(h) = m_N(ig)(iyr);           // numbers-at-age
                wa(h) = m_d3_wt_avg(ig)(iyr);   // weight-at-age
                mu(h) = exp(log(wa(h)/d_a(ig))/d_b(ig));  // mean size-at-age
                sd(h) = 0.1 * mu(h);                      // sd in mean size-at-age
                for(int k = 1; k <= nfleet; k++ )
                { 
                    int kk = m_nGearIndex(k);
                    // Implement size limits here.
                    pr = retention_probability(m_dLslim(k),m_dUslim(k),mu(h),sd(h));
                    pd = 1.0 - pr;

                    // int kk = nFleetIndex(k);
                    d3_Va(h)(k) = exp(d4_logSel(kk)(ig)(iyr));
                    d3_Da(h)(k) = d3_Va(h)(k);

                    // Joint probability model
                    d3_Va(h)(k)=elem_prod(d3_Va(h)(k),pr + pd*m_dDiscMortRate(k));
                    d3_Da(h)(k)=elem_prod(d3_Da(h)(k),pd);
                    
                }
            }  // nsex

            // Calculate instantaneous fishing mortality rates.
            dvector ft = cBCE.getFishingMortality(ct,ma,&d3_Va,na,wa,_hCt);
            
            // cout<<"fishing rate "<<ft<<endl;

            // Fill m_dCatchData array with actual catches taken by each fleet.
            for(int k = 1; k <= nfleet; k++ )
            {
                // cout<<"DFT counter +"<<ft_count<<" "<<ft_counter<<endl;
                m_log_ft_pars(++ft_counter) = ft(k);
                // Calculate total mortality array.
                for(int h = 1; h <= nsex; h++ )
                {
                    int ig = pntr_ags(f,g,h);
                    m_F(ig)(iyr) += ft(k) * d3_Va(h)(k);
                    m_ft(ig)(k)(iyr) = ft(k);

                    // Calculate wastage based on size-limits.
                    dvector fa = ft(k)*d3_Da(h)(k);
                    dvector za = ma(h) + fa;
                    dvector ba = elem_prod(na(h),wa(h));
                    _hDt(h)   += elem_prod(elem_div(fa,za),1.0-exp(-za)) * ba;
                    _hWt(h)   += elem_prod(elem_div(fa,za),1.0-exp(-za)) 
                                 * ba*m_dDiscMortRate(k);
                }
                // cout<<"Sublegal Discards \n"<<_hDt<<endl;
                

                if( ft(k) > 0 )
                {
                    int kk = nFleetIndex(k);
                    int hh = m_nCSex(k);   // flag for sex
                    for( h = 1; h <= hh+1; h++ )
                    {
                        m_irow ++;
                        m_dCatchData(m_irow,1) = iyr;
                        m_dCatchData(m_irow,2) = kk;
                        m_dCatchData(m_irow,3) = f;
                        m_dCatchData(m_irow,4) = g;
                        m_dCatchData(m_irow,5) = hh>0?h:0;
                        m_dCatchData(m_irow,6) = 1;  //TODO: set type of catch
                        m_dCatchData(m_irow,7) = hh>0?_hCt(h,k):colsum(_hCt)(k);
                        m_dCatchData(m_irow,8) = 0.02;  //TODO: add real log_se for catch obs.

                        m_dSubLegalData(m_irow,1) = iyr;
                        m_dSubLegalData(m_irow,2) = kk;
                        m_dSubLegalData(m_irow,3) = f;
                        m_dSubLegalData(m_irow,4) = g;
                        m_dSubLegalData(m_irow,5) = hh>0?h:0;
                        m_dSubLegalData(m_irow,6) = hh>0?_hDt(h,k):colsum(_hDt)(k);
                        m_dSubLegalData(m_irow,7) = hh>0?_hWt(h,k):colsum(_hWt)(k);
                        double effsex = _hCt(h,k)/(_hCt(h,k)+_hDt(h,k));
                        double effnsx = colsum(_hCt)(k)/(colsum(_hCt)(k)+colsum(_hDt)(k));
                        m_dSubLegalData(m_irow,8) = hh>0?effsex:effnsx;
                    }
                }
            }

        }  // ngroup g
    } // narea f
    // cout<<m_dCatchData<<endl;
    // cout<<"END"<<endl;
    //cout<<"finished implementing fisheries"<<endl;

}




/**
 * @brief Calculate total mortality rates
 * @details Total mortality rates based on the sum of  natural mortality
 * fishing mortality and discard mortality, including wastage from 
 * directed fisheries that have size-limit regulations in effect.
 * 
 * @param iyr Current year.
 * 
 * TODO. Add the Discard mortality rate component.
 *     Z = M + F + D
 */
void OperatingModel::calcTotalMortality(const int& iyr)
{
    for(int ig = 1; ig <= n_ags; ig++ )
    {
        m_Z(ig)(iyr) = m_M(ig)(iyr) + m_F(ig)(iyr);
        m_S(ig)(iyr) = 1.0 - exp( -m_Z(ig)(iyr) );
    }
}

void OperatingModel::calcRelativeAbundance(const int& iyr)
{
    //m_d3SurveyData
    // Survey data header:
    // 1    2      3     4     5      6    7   8    9
    // iyr  index  gear  area  group  sex  se  pe     timing
    static int irow = 0;
    irow ++;
    int gear;
    dvector na(sage,nage);
    dvector va(sage,nage);
    dvector sa(sage,nage);
    dvector ma(sage,nage);
    dvector wa(sage,nage);
    double dV;

    for(int k = 1; k <= nItNobs; k++ )
    {
        gear = d3_survey_data(k)(1)(3);
        for( f = 1; f <= narea; f++ )
        {
            for( g = 1; g <= ngroup; g++ )
            {
                dV = 0;
                for( h = 1; h <= nsex; h++ )
                {
                    int ig = pntr_ags(f,g,h);
                    va = exp(d4_logSel(gear)(ig)(iyr));
                    wa = m_d3_wt_avg(ig)(iyr);
                    ma = m_d3_wt_mat(ig)(iyr);
                    //sa  //TODO correct for survey timing.
                    na = m_N(ig)(iyr);
                    switch(n_survey_type(k))
                    {
                        case 1: // vulnerable numbers
                            dV += na*va;
                        break;

                        case 2: // vulnerable biomass
                            dV += elem_prod(na,va)*wa;
                        break;

                        case 3: // spawning biomass
                            dV += na*ma;
                        break;
                    }
                }
                // Relative abundance index.
                double sd = m_dSigma(g);
                // cout<<"Awe man"<<endl;
                double it = m_q(k)*dV*exp(m_epsilon(k,iyr)*sd); 
                // cout<<"q "<<m_q(k)<<endl;
                // V is the population that is proportional to the index.
                m_d3SurveyData(k)(n_it_nobs(k)+irow,1) = iyr;
                m_d3SurveyData(k)(n_it_nobs(k)+irow,2) = it;
                m_d3SurveyData(k)(n_it_nobs(k)+irow,3) = gear;
                m_d3SurveyData(k)(n_it_nobs(k)+irow,4) = f;    //TODO add to MSE controls
                m_d3SurveyData(k)(n_it_nobs(k)+irow,5) = g;    //TODO add to MSE controls
                m_d3SurveyData(k)(n_it_nobs(k)+irow,6) = 0;    //TODO add to MSE controls
                m_d3SurveyData(k)(n_it_nobs(k)+irow,7) = 0.2;  //TODO add to MSE controls
                m_d3SurveyData(k)(n_it_nobs(k)+irow,8) = 0.01; //TODO add to MSE controls
                m_d3SurveyData(k)(n_it_nobs(k)+irow,9) = 0.5;  //TODO add to MSE controls
                
            }
        }
    }   // end loop over surveys

    //cout<<"finished calculating relative abundance"<<endl;

}

/**
 * @brief Composition data.
 * @details Calculate composition data for current year.
 * 
 * Loop over nAgears
 *  Loop over areas
 *    Loop over groups
 *      Loop over sex.
 *          -Determine the catch-at-age proportions.
 *           given Selectivity of the gear & abundance
 *           in area f for group g and sex h.
 *          -Sample catch with a precision specified in
 *           the MSE control file.  Can either use a 
 *           multinomial sample, or multivariate logistic.
 *          -Other ideas include 2-stage sampling with the
 *           probability of finding an aggregation, and 
 *           the composition of the aggregation is highly
 *           correlated.
 * 
 * @param iyr Current year.
 * @param m_d3_A.allocate(1,nAgears,1,m_n_A_nobs,n_A_sage-5,n_A_nage); age-comp array.
 */     
void OperatingModel::calcCompositionData(const int& iyr)
{
    int gear;
    dvector na(sage,nage);
    dvector va(sage,nage);
    dvector fa(sage,nage);
    dvector ma(sage,nage);
    dvector za(sage,nage);
    dmatrix ca(1,nsex,sage,nage);
    dmatrix pa(1,nsex,sage,nage);
    double ft;
    for(int k = 1; k <= nAgears; k++ )
    {
        if( m_n_A_nobs(k) ) 
        {
            gear = m_d3_A(k)(1)(n_A_sage(k)-5);
            for(int f = 1; f <= narea; f++ )
            {
                for(int g = 1; g <= ngroup; g++ )
                {
                    ca.initialize();

                    for(int h = 1; h <= nsex; h++ )
                    {
                        int ig = pntr_ags(f,g,h);
                        va = exp(d4_logSel(gear)(ig)(iyr));
                        na = m_N(ig)(iyr);
                        ma = m_M(ig)(iyr);
                        ft = m_ft(ig)(gear)(iyr);
                        fa = (ft>0?ft:1.0) * va;  
                        za = ma + fa;
                        ca(h) = elem_prod(elem_prod(elem_div(fa,za),1.-exp(-za)),na);
                        pa(h) = ca(h) / sum(ca(h));
                        pa(h) = rmvlogistic(pa(h),m_nATau(k),m_nSeed+iyr);
                        //rmvlogistic(pa(h),m_nATau,m_nSeed+iyr);
                    }
                    
                    int hh = m_nASex(k);   // flag for sex
                    for( h = 1; h <= hh+1; h++ )
                    {
                        // cout<<hh<<endl;
                        m_A_irow(k) ++;
                        m_d3_A(k)(n_A_nobs(k)+m_A_irow(k),n_A_sage(k)-6) = iyr;
                        m_d3_A(k)(n_A_nobs(k)+m_A_irow(k),n_A_sage(k)-5) = gear;
                        m_d3_A(k)(n_A_nobs(k)+m_A_irow(k),n_A_sage(k)-4) = f;
                        m_d3_A(k)(n_A_nobs(k)+m_A_irow(k),n_A_sage(k)-3) = g;
                        m_d3_A(k)(n_A_nobs(k)+m_A_irow(k),n_A_sage(k)-2) = hh>0?h:0;
                        m_d3_A(k)(n_A_nobs(k)+m_A_irow(k),n_A_sage(k)-1) = 1;  //age_err

                        m_d3_A(k)(n_A_nobs(k)+m_A_irow(k))(n_A_sage(k),n_A_nage(k))
                        = hh>0?pa(h)(n_A_sage(k),n_A_nage(k)):colsum(pa)(n_A_sage(k),n_A_nage(k));
                    }
                    
                }
            }
        }
    }
    
}

void OperatingModel::calcEmpiricalWeightAtAge(const int& iyr)
{
    int gear;
    
    for(int k = 1; k <= nWtTab; k++ )
    {
        if( m_nWtNobs(k) )
        {

            gear = m_d3_inp_wt_avg(k)(1)(sage-4);

            for(int f = 1; f <= narea; f++ )
            {
                for(int g = 1; g <= ngroup; g++ )
                {

                    int hh = m_nWSex(gear);   // flag for sex
                    for( h = 1; h <= hh+1; h++ )
                    {
                        m_W_irow(k) ++;
                        int ig = pntr_ags(f,g,h);
            
                        m_d3_inp_wt_avg(k)(nWtNobs(k)+m_W_irow(k))(sage-5) = iyr;
                        m_d3_inp_wt_avg(k)(nWtNobs(k)+m_W_irow(k))(sage-4) = gear;
                        m_d3_inp_wt_avg(k)(nWtNobs(k)+m_W_irow(k))(sage-3) = f;
                        m_d3_inp_wt_avg(k)(nWtNobs(k)+m_W_irow(k))(sage-2) = g;
                        m_d3_inp_wt_avg(k)(nWtNobs(k)+m_W_irow(k))(sage-1) = hh>0?h:0;
                        
                        m_d3_inp_wt_avg(k)(nWtNobs(k)+m_W_irow(k))(sage,nage) =m_d3_wt_avg(ig)(iyr)(sage,nage);
                    }
                }
            }
        }
    }    
}



void OperatingModel::updateReferenceModel(const int& iyr)
{

    // compute spawning biomass at time of spawning.
    dvector  stmp(sage,nage); stmp.initialize();

    for(int f=1;f<=narea;f++)
    {
        for(int h=1;h<=nsex;h++)
        {
            for(int g = 1; g<=ngroup; g++)
            {
                int ig = pntr_ags(f,g,h);
                    
                stmp      = mfexp(-m_Z(ig)(iyr)*d_iscamCntrl(13));
                m_sbt(iyr,g) += elem_prod(m_N(ig)(iyr),m_d3_wt_mat(ig)(iyr)) * stmp;
                m_bt(iyr,g) = elem_prod(m_N(ig)(iyr),m_d3_wt_avg(ig)(iyr)) * stmp;                 
            }
        }
    }

    dvector tmp_rec(1,narea);tmp_rec.initialize();
    dvector tmp_rec_dis(1,narea);tmp_rec_dis.initialize();
    dvector prop_rec_g(1,n_ags);prop_rec_g.initialize();


    for(int ig = 1; ig <= n_ags; ig++ )
    {
        int f  = n_area(ig);
        int g  = n_group(ig);
        int ih = pntr_ag(f,g);
        
        
        // Recruitment
        // three options : average recruitment, Beverton &Holt and Ricker
        // m_delta is the process error.
        
        double tmp_st;
        tmp_st = m_sbt(iyr-sage,g);

        switch(m_nRecType)
        {
            case 1:  // | Beverton Holt model
                m_dbeta(g) = (m_dKappa(g)-1.0)/(mv.sbo(g));
                m_N(ig)(iyr+1,sage) = mv.so(g)*tmp_st/ (1.+m_dbeta(g)*tmp_st);
            break;

            case 2:  // | Ricker model
                m_dbeta(g) = log(m_dKappa(g))/mv.sbo(g);
                m_N(ig)(iyr+1,sage) = mv.so(g)*tmp_st*exp(-m_dbeta(g)*tmp_st);
            break;

            case 3: // average recruitment
                m_N(ig)(iyr+1,sage) = m_dRbar(ih);
        }
        m_N(ig)(iyr+1,sage) = m_N(ig)(iyr+1,sage)/nsex; //* mfexp( mv.log_rbar(ih));

        // add process errors assuming each area/sex has the same deviation.
        m_N(ig)(iyr+1,sage) *= exp( m_dTau(g)*m_delta(g,iyr) - 0.5*square(m_dTau(g)) );
        // COUT(m_dTau);
        //disperse the recruits in each year 
        // assumes all groups disperse the same and prop of groups by area remain const
        // TODO allow for separate dispersal matrices for each group
        
        //1 - calculate total recruits per area: tmp_rec
        tmp_rec(f) += m_N(ig)(iyr+1,sage);
        
    }
        
    //disperse recruits 
    tmp_rec_dis = tmp_rec*m_dispersal;

    for(int ig = 1; ig <= n_ags; ig++ )
    {
        int f  = n_area(ig);
        prop_rec_g(ig)=m_N(ig)(iyr+1,sage)/tmp_rec(f);
             
        m_N(ig)(iyr+1,sage) = (tmp_rec(f)/nsex)*prop_rec_g(ig); 
    
        // Update numbers-at-age
        dvector st = exp(-(m_M(ig)(iyr)+m_F(ig)(iyr)) );
        m_N(ig)(iyr+1)(sage+1,nage) = ++ elem_prod(m_N(ig)(iyr)(sage,nage-1)
                                                   ,st(sage,nage-1));
        m_N(ig)(iyr+1,nage)        += m_N(ig)(iyr,nage) * st(nage); 
        
    }

     //cout<<"finished updatinf ref pop"<<endl;
}


void OperatingModel::writeParameterFile(const int& iyr)
{
    cout<<"Writing TRUE.par"<<endl;
    ofstream pfs("TRUE.par");
    pfs << mv.log_ro     << endl;
    pfs << mv.steepness  << endl;
    pfs << mv.m          << endl;
    pfs << mv.log_rbar   << endl;
    pfs << mv.log_rinit  << endl;
    pfs << mv.rho        << endl;
    pfs << mv.varphi     << endl;
    pfs <<"# slx_log_par" << endl;
    pfs << *mv.d3_slx_log_par << endl;
    pfs <<"# second instance until sel_par is deprecated" << endl;
    pfs << *mv.d3_log_sel_par << endl;

    pfs <<"# log_ft_pars " << endl;
    
    // fishing mortality rate parameters for each gear.
    int n = nCtNobs +(iyr-nyr)*m_nn;
    pfs << m_log_ft_pars(1,n)  <<endl;

    pfs <<"# init_log_rec_devs:" << endl;
    pfs << mv.init_log_rec_devs  << endl;

    pfs <<"# log_rec_devs:"      << endl;

    // add simulation devs m_delta * sigma_r
    for( g = 1; g <= ngroup; g++ )
    {
        pfs << mv.log_rec_devs(g)(syr,nyr-sage-1);
        for(int i=nyr-sage; i<=iyr; i++)
        {
            pfs << " " << m_delta(g)(i) * m_dTau;
        }
        pfs << endl;
    }

    pfs <<"# log_q_devs:  "      << endl;

    // now add additional qdevs
    for(k=1;k<=nItNobs;k++)
    {
        pfs << mv.log_q_devs(k);
        for(i = nyr+1; i<= iyr; i++)
        {
            pfs<<" "<< 0;
        }
        pfs << endl;
    }

    pfs <<"# log_age_tau2 "      << endl;
    pfs << mv.log_age_tau2       << endl;

    pfs << mv.phi1 <<endl;
    pfs << mv.phi2 <<endl;
    pfs << mv.log_degrees_of_freedom << endl;
}


/**
 * @brief Write a new data file for iscam
 * @details This routine writes the data file for iSCAM
 *  Note that if a prospective years is > 0, the syr will be greater
 *  than many of the data sets, and iSCAM will not work. Therfore, the data
 *  structures, have to be written from syr-iyr only.
 *  
 *  Or, for syr write ,syr - (int)d_iscamCntrl(14)
 * @param iyr the terminal year of the data in the data file.
 */
void OperatingModel::writeDataFile(const int& iyr)
{
        cout<<"WRITING SIMULATED DATA"<<endl;
        adstring sim_datafile_name = "Simulated_Data_"+str(rseed)+".dat";
        ofstream dfs(sim_datafile_name);
        dfs<<"#Model dimensions"        <<endl;
        dfs<< narea                     <<endl;
        dfs<< ngroup                    <<endl;
        dfs<< nsex                      <<endl;
        dfs<< syr-(int)d_iscamCntrl(14) <<endl;
        dfs<< iyr                       <<endl;
        dfs<< sage                      <<endl;
        dfs<< nage                      <<endl;
        dfs<< ngear                     <<endl;
         
        dfs<<"#Allocation"              <<endl;
        dfs<< dAllocation               <<endl;
        
        // Write age-schedule information
        dfs<<"#Age-schedule and population parameters"<<endl;
        dfs<< d_linf            <<endl;
        dfs<< d_vonbk           <<endl;
        dfs<< d_to              <<endl;
        dfs<< d_a               <<endl;
        dfs<< d_b               <<endl;
        dfs<< d_ah              <<endl;
        dfs<< d_gh              <<endl;
        dfs<< n_MAT             <<endl;
        dfs<< d_maturityVector  <<endl;
    
        // Write aging error definitions
        dfs<<"#Aging Error Vectors" <<endl;
        dfs<< n_age_err             <<endl;
        dfs<< age_err               <<endl;

        // Write catch array
        dfs<<"#Observed catch data"             <<endl;
        int tmp_nCtNobs = nCtNobs+(iyr-nyr)*m_nn;

        dfs<< nCtNobs + (iyr-nyr)*m_nn          <<endl; 
        dfs<< m_dCatchData.sub(1,tmp_nCtNobs)   <<endl;
    
        // Write relative abundance indices
        dfs<<"#Abundance indices"       <<endl;
        ivector tmp_n_it_nobs(1,nItNobs);
        tmp_n_it_nobs.initialize();
        d3_array tmp_d3SurveyData(1,nItNobs,1,tmp_n_it_nobs,1,8);
        tmp_d3SurveyData.initialize();
        
        for(int k=1;k<=nItNobs;k++)
        {
            tmp_n_it_nobs(k)    = n_it_nobs(k) + (iyr-nyr);
            tmp_d3SurveyData(k) = m_d3SurveyData(k).sub(1,tmp_n_it_nobs(k));
        }
        
        dfs<< nItNobs                   <<endl;
        dfs<< tmp_n_it_nobs             <<endl;
        dfs<< n_survey_type             <<endl;
        dfs<< tmp_d3SurveyData          <<endl;
            
        // Write Composition information
        dfs<<"#Age composition"     <<endl;
        ivector tmp_n_A_nobs(1,nAgears);
        tmp_n_A_nobs.initialize();

        d3_array tmp_d3_A(1,nAgears,1,tmp_n_A_nobs,n_A_sage-6,n_A_nage); //n_A_sage is a vector!!
        tmp_d3_A.initialize();
        
        for(int k=1;k<=nAgears;k++)
        {
            if(n_A_nobs(k) > 0)
            {
                tmp_n_A_nobs(k) = n_A_nobs(k) + (iyr-nyr) + (iyr-nyr) * m_nASex(k);
                tmp_d3_A(k) = m_d3_A(k).sub(1,tmp_n_A_nobs(k));    
            }
        }
        
        dfs<< nAgears               <<endl;
        dfs<< tmp_n_A_nobs          <<endl;
        dfs<< n_A_sage              <<endl;
        dfs<< n_A_nage              <<endl;
        dfs<< inp_nscaler           <<endl;
        dfs<< n_ageFlag             <<endl;
        dfs<< tmp_d3_A              <<endl;
    
        // Write Empirical weight-at-age data.
        // Issue 30.  Bug in writing empirical weight-at-age.
        // data where year in the first row is 0.
        dfs<<"#Empirical weight-at-age data"    <<endl;
        dfs<< nWtTab                    <<endl;
        ivector tmp_nWtNobs(1,nWtTab);
        for( k = 1; k <= nWtTab; k++ )
        {
            if(nWtNobs(k) > 0)
                tmp_nWtNobs(k)= nWtNobs(k) + (iyr-nyr) + (iyr-nyr) * m_nWSex(k);
        }

        d3_array tmp_d3_inp_wt_avg(1,nWtTab,1,tmp_nWtNobs,sage-5,nage);
        for(int k=1;k<=nWtTab;k++)
        {           
            
            tmp_d3_inp_wt_avg(k)= m_d3_inp_wt_avg(k).sub(1,tmp_nWtNobs(k)) ;
            if(projwt(k)<0)
            {
                for(int ii=1; ii<=abs(projwt(k)); ii++)
                {
                    tmp_d3_inp_wt_avg(k)(ii)(sage-5) = fabs(tmp_d3_inp_wt_avg(k)(1)(sage-5))*(-1);
                }
                
            }
            
        }
        
        //cout<<tmp_d3_inp_wt_avg(1)<<endl;
        dfs<< tmp_nWtNobs               <<endl;
        dfs<< tmp_d3_inp_wt_avg         <<endl; 
    
        dfs<<"#EOF" <<endl;
        dfs<< 999   <<endl;
        //exit(1);

}

void OperatingModel::runStockAssessment()
{
    /*
        Call iscam ensuring that the correct Data files and control files 
        for the estimator are specfied in the mseRUN.dat file.

        System call should reflect.

        make ARG="-ind mseRUN.dat"  on unix flavors

        system("iscam.exe -ind mseRUN.dat")  on windoze boxes.
    */

        
        ofstream rd("mseRUN.dat");
        rd<<"Simulated_Data_"+str(rseed)+".dat"<<endl;
        rd << m_controlFile        <<endl;
        rd << ProjectFileControl   <<endl;
        rd << ProcedureControlFile <<endl;
        rd << ScenarioControlFile  <<endl;
        //exit(1);

 
    switch( int(m_nAssessOpt) ) // option read in from .mpc file
    {
        case 0:
            cout<<"Perfect information scenario"<<endl;

            #if defined __APPLE__ || defined __linux
            system("./iscam -ind mseRUN.dat -maxfn 0 -nox -nohess -ainp TRUE.par -phase 4");
            #endif

        break;
        
        case 1:        
            cout<<"running stock assessment"<<endl;

            #if defined __APPLE__ || defined __linux
            // cout<<m_est_fmsy<<endl;
            system("./iscam -ind mseRUN.dat -nox > /dev/null 2>&1");

            #endif

            #if defined _WIN32 || defined _WIN64

            system("iscam.exe -ind mseRUN.dat");

            #endif
        break;
    }
}


/**
 * @brief Append true state variables to iscam.rep.
 * @details This routine is called at the end of the simulation model and appends the
 * true state variables to the iscam.rep file.
 * 
 * Need to write out simulation scenarios: low recruitment, high recruitment, size-at-age
 * 
 */
void OperatingModel::writeSimulationVariables()
{

    ofstream report("milka.rep");
    REPORT(m_PDO_phase);
    REPORT(m_SAA_flag);
    REPORT(m_nRecType);
    REPORT(m_dBo);
    REPORT(m_sbt);
    REPORT(m_dCatchData);
    REPORT(m_dSubLegalData);
    REPORT(m_ft);
    REPORT(m_yr);

}
