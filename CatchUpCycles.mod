@#ifdef dynareOBC
    @#define Estimation = 0
@#else
    @#define Estimation = 1
@#endif

@#includepath "DynareTransformationEngine"

@#include "Initialize.mod"
@#define UsingGrowthSyntax = 1

@#define ShockProcessNames = [ "Psi", "Omega", "Lambda1", "LambdaS", "LambdaT", "LambdaM", "ZS", "ZT", "beta" ]

@#for ShockProcessName in ShockProcessNames
    parameters rho_@{ShockProcessName} sigma_@{ShockProcessName};
    rho_@{ShockProcessName} = 0.9;
    sigma_@{ShockProcessName} = 0.01;
    @#define ShockProcesses = ShockProcesses + [ ShockProcessName, "0", "Inf", ShockProcessName + "_STEADY", "rho_" + ShockProcessName, "sigma_" + ShockProcessName ]
@#endfor

@#define OtherShockProcessNames = [ "Rs", "GAs", "GPdollar" ]

@#for ShockProcessName in OtherShockProcessNames
    varexo epsilon_@{ShockProcessName};
    @#define EndoVariables = EndoVariables + [ ShockProcessName, "0", "Inf", "1" ]
@#endfor

@#for Lag in 1:8
    var log_Rs_LAG_@{Lag};
@#endfor

// GAs process calibrated from http://www.frbsf.org/economic-research/indicators-data/total-factor-productivity-tfp/
// Correct sigma for GAs = 0.00824165945443227

@#define PureTrendEndoVariables = PureTrendEndoVariables + [ "As", "GAs" ]
@#define PureTrendEndoVariables = PureTrendEndoVariables + [ "Pdollar", "GPdollar" ]

@#include "CreateShocks.mod"

@#define EndoVariables = EndoVariables + [ "GYTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GBTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GkappaXTrend", "0", "Inf", "1" ]

@#define EndoVariables = EndoVariables + [ "C", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "X", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "IS", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "IT", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "KS", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "KT", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "PT", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "W", "0", "Inf", "GYTrend" ]

@#define EndoVariables = EndoVariables + [ "Bs", "-Inf", "Inf", "GBTrend" ]

@#define EndoVariables = EndoVariables + [ "A", "0", "Inf", "GAs" ]

@#define EndoVariables = EndoVariables + [ "kappaX", "0", "Inf", "GkappaXTrend" ]

@#define EndoVariables = EndoVariables + [ "sm", "-Inf", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "R", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "Jd", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "Td", "-Inf", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "LS", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "LT", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "LM", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "QS", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "QT", "0", "Inf", "1" ]

@#define EndoVariables = EndoVariables + [ "GC", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GI", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GW", "0", "Inf", "1" ]

@#include "ClassifyDeclare.mod"

// The following block of params are set in the steady-state model block.
parameters Psi_STEADY Omega_STEADY Lambda1_STEADY LambdaS_STEADY LambdaT_STEADY LambdaM_STEADY ZS_STEADY ZT_STEADY alphaS beta_STEADY Jd_STEADY lambda Rs_STEADY;

Psi_STEADY = 0;
Omega_STEADY = 0;
Lambda1_STEADY = 0;
LambdaS_STEADY = 0;
LambdaT_STEADY = 0;
LambdaM_STEADY = 0;
ZS_STEADY = 0;
ZT_STEADY = 0;
alphaS = 0;
beta_STEADY = 0;
Jd_STEADY = 0;
lambda = 0;
Rs_STEADY = 0;

// The following block of params are targets on steady-states of endogenous variables.
parameters A_STEADY U_OVER_C_STEADY LM_STEADY WTEFF_OVER_WS_STEADY WM_OVER_WS_TARGET RESEARCH_SHARE_TARGET CONSUMPTION_SHARE_TARGET FIXED_COST_SHARE_TARGET LABOUR_SHARE_TARGET;

A_STEADY = 0.75;                                 // Spain/US, BCL TFP Database http://www.longtermproductivity.com/download.html
U_OVER_C_STEADY = 0.4;                           // Smets Wouters (2003)
LM_STEADY = 0.1;                                 // ONS Annual Population Survey (via Nomis) - "corporate managers and directors" + "other managers and proprietors"
WTEFF_OVER_WS_STEADY = 2.3;                      // EUROMOD https://www.euromod.ac.uk/using-euromod/statistics Decile 9 of Original Income for Spain over Deciles 1 to 9
WM_OVER_WS_TARGET = 3.8;                         // EUROMOD https://www.euromod.ac.uk/using-euromod/statistics Decile 10 of Original Income for Spain over Deciles 1 to 9
RESEARCH_SHARE_TARGET = 0.012;                   // World Bank https://data.worldbank.org/indicator/GB.XPD.RSDV.GD.ZS
CONSUMPTION_SHARE_TARGET = 0.6 / ( 0.6 + 0.22 ); // Smets Wouters (2003)
FIXED_COST_SHARE_TARGET = 0.2;                   // http://www.beta-umr7522.fr/IMG/UserFiles/Koebel/fixed_cost_aes_new_2.pdf
LABOUR_SHARE_TARGET = 0.7;                       // Smets Wouters (2003)

// The following block of params are set here.
parameters R_STEADY GPdollar_STEADY GAs_STEADY varsigma alphaT HC_SHARE nuS gamma tau eta deltaS deltaT PhiS2 PhiT2;

R_STEADY = ( 1.03 ) ^ ( 1 / 4 );                 // Steady-state real interest rate from Smets Wouters (2007)
GPdollar_STEADY = exp( -0.00361309247175487 );   // From data used for estimation.
GAs_STEADY = exp( 0.00319206470220157 );         // http://www.frbsf.org/economic-research/indicators-data/total-factor-productivity-tfp/
varsigma = 1.4;                                  // Smets Wouters (2003)
alphaT = 0.5;                                    //
HC_SHARE = 0;                                    // Only used in calibration.
nuS = 2.5;                                       // Smets Wouters (2003)
gamma = 0.999;                                   // Jaimovich Rebelo (2009)
tau = 1;                                         //
eta = 0.5;                                       //
deltaS = 0.025;                                  // Smets Wouters (2003) calibrate to 0.025
deltaT = 0.25;                                   // Smets Wouters (2003) calibrate to 0.025
PhiS2 = 3.5;                                     // Smets Wouters (2003)
PhiT2 = 3.5;                                     // Smets Wouters (2003)

// The following parameters are estimated.
parameters theta nuT nuM;

nuT = nuS;
nuM = nuS;
theta = 0.01;

load_params_and_steady_state( 'SteadyState.txt' );

@#if Estimation
    varobs log_Rs log_GPdollar log_GAs log_GC log_GI log_GW;

    estimated_params;
        @#for ShockProcessName in ShockProcessNames
            rho_@{ShockProcessName}, rho_@{ShockProcessName}, beta_pdf, 0.5, sqrt( 1 / 20 );
            sigma_@{ShockProcessName}, sigma_@{ShockProcessName}, gamma_pdf, (1/2)*0.01+(1/2)*sqrt(0.01^2+4*1^2), 1;
        @#endfor
        nuT, nuT, gamma_pdf, (1/2)*nuS+(1/2)*sqrt(nuS^2+4*nuS^2), nuS;
        nuM, nuM, gamma_pdf, (1/2)*nuS+(1/2)*sqrt(nuS^2+4*nuS^2), nuS;
        theta, theta, gamma_pdf, (1/2)*0.01+(1/2)*sqrt(0.01^2+4*0.01^2), 0.01;
    end;
@#endif

model;

    @#include "InsertNewModelEquations.mod"
    
    #E = Pdollar * Bs - Pdollar * Rs_LAG * Bs_LAG;
    #E_LEAD = Pdollar_LEAD * Bs_LEAD - Pdollar_LEAD * Rs * Bs;
    #Y = C + IS + IT + theta / 2 * ( Pdollar * Bs ) ^ 2 / As + E;
    #Y_LEAD = C_LEAD + IS_LEAD + IT_LEAD + theta / 2 * ( Pdollar_LEAD * Bs_LEAD ) ^ 2 / As_LEAD + E_LEAD;
    #mu = lambda * eta * Jd / ( Jd - ( 1 - eta ) );
    #mu_LAG = lambda * eta * Jd_LAG / ( Jd_LAG - ( 1 - eta ) );
    #PS = A / ( 1 + mu_LAG );
    #PS_LEAD = A_LEAD / ( 1 + mu );
    #S = Y / A;
    #S_LEAD = Y_LEAD / A_LEAD;
    #T = Jd * Td;
    #T_LEAD = Jd_LEAD * Td_LEAD;
    LM = Jd * Psi;
    //#LM_LEAD = Jd_LEAD * Psi_LEAD;
    #SRS = alphaS * PS * S / KS_LAG;
    #SRS_LEAD = alphaS * PS_LEAD * S_LEAD / KS;
    #SRT = alphaT * PT * T / KT_LAG;
    #SRT_LEAD = alphaT * PT_LEAD * T_LEAD / KT;
    #omegad = Jd * ( 1 - eta ) / ( ( Jd - ( 1 - eta ) ) ^ 2 * ( 1 + mu ) );
    #sdd = 1 - omegad / ( 1 + omegad ) * ( lambda - mu ) * ( mu - eta * lambda ) / ( lambda * ( 1 - eta ) * mu );
    #U = C - X_LAG * ( Lambda1 + LambdaS * LS ^ ( 1 + nuS ) / ( 1 + nuS ) + LambdaT * LT ^ ( 1 + nuT ) / ( 1 + nuT ) + LambdaM * LM ^ ( 1 + nuM ) / ( 1 + nuM ) );
    #U_LEAD = C_LEAD - X * ( Lambda1_LEAD + LambdaS_LEAD * LS_LEAD ^ ( 1 + nuS ) / ( 1 + nuS ) + LambdaT_LEAD * LT_LEAD ^ ( 1 + nuT ) / ( 1 + nuT ) + LambdaM_LEAD * LM_LEAD ^ ( 1 + nuM ) / ( 1 + nuM ) );
    #kappaB = U ^ ( -varsigma ) - ( 1 - gamma ) * kappaX * X / C;
    #kappaB_LEAD = U_LEAD ^ ( -varsigma ) - ( 1 - gamma ) * kappaX_LEAD * X_LEAD / C_LEAD;
    #Xi_LEAD = beta * kappaB_LEAD / kappaB;
    #WS = 1 / kappaB * U ^ ( -varsigma ) * X_LAG * LambdaS * LS ^ nuS;
    #WT = 1 / kappaB * U ^ ( -varsigma ) * X_LAG * LambdaT * LT ^ nuT;
    #WM = 1 / kappaB * U ^ ( -varsigma ) * X_LAG * LambdaM * LM ^ nuM;
    
    log_Rs_LAG_1 = log_Rs(-1) - log( Rs_STEADY );
    @#for Lag in 2:8
        log_Rs_LAG_@{Lag} = log_Rs_LAG_@{Lag-1}(-1);
    @#endfor

    log_Rs - log( Rs_STEADY )             = 1.41947285017037 * log_Rs_LAG_1	- 0.801403990130289 * log_Rs_LAG_2 + 0.742975932948491 * log_Rs_LAG_3 - 0.511444343607851 * log_Rs_LAG_4 + 0.269731846458521 * log_Rs_LAG_5 - 0.180165905325992 * log_Rs_LAG_6 - 0.199767480215731 * log_Rs_LAG_7 + 0.358115595633627 * log_Rs_LAG_8 - 0.114516339974676 * log_Rs_LAG_8(-1) + sqrt( 1.58897610763807e-06 ) * epsilon_Rs;
    log_GAs - log( GAs_STEADY )           = 0.00824198231591593 * epsilon_GAs;
    log_GPdollar - log( GPdollar_STEADY ) = 0.332339490649165 * ( log_GPdollar(-1) - log( GPdollar_STEADY ) ) + sqrt( 0.00160795066254317 ) * epsilon_GPdollar;
    
    GYTrend = GAs ^ ( 1 / ( 1 - alphaS ) );
    GBTrend = GYTrend / GPdollar;
    // GSTrend = GAs ^ ( alphaS / ( 1 - alphaS ) );
    GkappaXTrend = GYTrend ^ ( -varsigma );
    
    GC = C / C_LAG;
    GI = ( IS + IT ) / ( IS_LAG + IT_LAG );
    W = ( WS * LS + WT * LT + WM * LM ) / ( LS + LT + LM );
    GW = W / W_LAG;
    
    S = KS_LAG ^ alphaS * ( ZS * LS ) ^ ( 1 - alphaS );
    T = ( KT_LAG / As ) ^ alphaT * ( ZT * LT ) ^ ( 1 - alphaT );
    WS * LS = ( 1 - alphaS ) * PS * S;
    WT * LT = ( 1 - alphaT ) * PT * T;
    0 = min( sm, Td );
    
    A = ( A_LAG ^ tau + ( As_LAG ^ tau - A_LAG ^ tau ) * Omega * Td_LAG / ( 1 + Omega * Td_LAG ) ) ^ ( 1 / tau );
    KS = ( 1 - deltaS ) * KS_LAG + ( 1 - PhiS2 / 2 * ( IS / IS_LAG - 1 ) ^ 2 ) * IS;
    KT = ( 1 - deltaT ) * KT_LAG + ( 1 - PhiT2 / 2 * ( IT / IT_LAG - 1 ) ^ 2 ) * IT;
    X = C ^ ( 1 - gamma ) * X_LAG ^ gamma;
    
    1 / Jd * mu / ( 1 + mu ) * Xi_LEAD * Y_LEAD * sdd / mu / tau * Omega_LEAD * ( As ^ tau - A ^ tau ) / A_LEAD ^ tau / ( 1 + Omega_LEAD * Td ) ^ 2 = PT * ( 1 - sm );
    1 / Jd * mu / ( 1 + mu ) * Xi_LEAD * Y_LEAD = Psi * WM + Td * PT;
    1 = Xi_LEAD * ( SRS_LEAD + QS_LEAD * ( 1 - deltaS ) ) / QS;
    1 = Xi_LEAD * ( SRT_LEAD + QT_LEAD * ( 1 - deltaT ) ) / QT;
    1 = QS * ( 1 - PhiS2 / 2 * ( IS / IS_LAG - 1 ) ^ 2 - PhiS2 * ( IS / IS_LAG - 1 ) * IS / IS_LAG ) + Xi_LEAD * QS_LEAD * PhiS2 * ( IS_LEAD / IS - 1 ) * ( IS_LEAD / IS ) ^ 2;
    1 = QT * ( 1 - PhiT2 / 2 * ( IT / IT_LAG - 1 ) ^ 2 - PhiT2 * ( IT / IT_LAG - 1 ) * IT / IT_LAG ) + Xi_LEAD * QT_LEAD * PhiT2 * ( IT_LEAD / IT - 1 ) * ( IT_LEAD / IT ) ^ 2;
    kappaX = beta * ( gamma * kappaX_LEAD * X_LEAD / X + U_LEAD ^ ( -varsigma ) * ( Lambda1_LEAD + LambdaS_LEAD * LS_LEAD ^ ( 1 + nuS ) / ( 1 + nuS ) + LambdaT_LEAD * LT_LEAD ^ ( 1 + nuT ) / ( 1 + nuT ) + LambdaM_LEAD * LM_LEAD ^ ( 1 + nuM ) / ( 1 + nuM ) ) );
    Pdollar * kappaB * ( 1 + theta * Pdollar * Bs / As ) = beta * Rs * Pdollar_LEAD * kappaB_LEAD;
    kappaB = beta * R * kappaB_LEAD;
    
end;

steady_state_model;
    // @#include "InsertNewStartSteadyStateEquations.mod"
    
    [ Jd_STEADY, lambda, alphaS ] = CalibrateModel( R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 );
    
    GPdollar_ = GPdollar_STEADY;
    Rs_STEADY = R_STEADY / GPdollar_;
    Rs_ = Rs_STEADY;
    R_ = R_STEADY;

    @#for Lag in 1:8
        log_Rs_LAG_@{Lag} = 0;
    @#endfor
    
    beta_STEADY = GAs_STEADY ^ ( varsigma / ( 1 - alphaS ) ) / R_STEADY;
                
    GAs_ = GAs_STEADY;
    beta_ = beta_STEADY;

    Jd_ = Jd_STEADY;
    A_ = A_STEADY;
    LM_ = LM_STEADY;

    GYTrend_ = GAs_ ^ ( 1 / ( 1 - alphaS ) );
    GBTrend_ = GYTrend_ / GPdollar_;
    GkappaXTrend_ = GYTrend_ ^ ( -varsigma );

    Bs_ = 0;
    sm_ = 0;

    Xi_LEAD_ = beta_ * GYTrend_ ^ ( -varsigma );

    QS_ = 1 / ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 - PhiS2 * ( GYTrend_ - 1 ) * GYTrend_ + Xi_LEAD_ * PhiS2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2 );
    QT_ = 1 / ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 - PhiT2 * ( GYTrend_ - 1 ) * GYTrend_ + Xi_LEAD_ * PhiT2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2 );

    S_ = 1;
    T_ = 1;
    
    Psi_ = LM_ / Jd_;
    
    Omega_ = A_ ^ tau * ( GAs_ ^ tau - 1 ) * Jd_ / ( 1 - A_ ^ tau * GAs_ ^ tau );
    mu_ = lambda * eta * Jd_ / ( Jd_ - ( 1 - eta ) );
    omegad_ = Jd_ * ( 1 - eta ) / ( ( Jd_ - ( 1 - eta ) ) ^ 2 * ( 1 + mu_ ) );
    sdd_ = 1 - omegad_ / ( 1 + omegad_ ) * ( lambda - mu_ ) * ( mu_ - eta * lambda ) / ( lambda * ( 1 - eta ) * mu_ );
    sfd_ = 1 - 1 / 2 * ( 1 + sdd_ / mu_ / tau ) * ( 1 - A_ ^ tau );
    PT_Tmp_ = ( 1 + ( ( Omega_ / Jd_ + sfd_ ) ^ 2 - sfd_ ^ 2 ) / A_ ^ tau ) / Omega_ / Psi_ * LM_ / sdd_ * mu_ * tau / ( A_ ^ ( -tau ) - 1 );
    PT_ = mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ / ( 1 + PT_Tmp_ );
    PS_ = A_ / ( 1 + mu_ );
    
    IS_ = Xi_LEAD_ * alphaS * PS_ * GYTrend_ / ( ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaS ) / GYTrend_ ) ) / QS_ / ( 1 - Xi_LEAD_ * ( 1 - deltaS ) );
    IT_ = Xi_LEAD_ * alphaT * PT_ * GYTrend_ / ( ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaT ) / GYTrend_ ) ) / QT_ / ( 1 - Xi_LEAD_ * ( 1 - deltaT ) );

    SRS_ = alphaS * PS_ * GYTrend_ / ( ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaS ) / GYTrend_ ) ) / IS_;
    SRT_ = alphaT * PT_ * GYTrend_ / ( ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaT ) / GYTrend_ ) ) / IT_;

    KS_ = ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaS ) / GYTrend_ ) * IS_;
    KT_ = ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaT ) / GYTrend_ ) * IT_;

    LT_Tmp_ = ( SRT_ * KT_ * HC_SHARE / alphaT + ( 1 - alphaT ) * PT_ ) / ( 1 - alphaS ) / PS_ / WTEFF_OVER_WS_STEADY;
    LT_ = LT_Tmp_ / ( 1 + LT_Tmp_ ) * ( 1 - LM_ );
    LS_ = 1 - LT_ - LM_;
    WM_OVER_WS_ = ( mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ - PT_ ) / ( 1 - alphaS ) / PS_ / LM_ * LS_;
        
    C_ = A_ - IS_ - IT_;
    X_ = C_ * GYTrend_ ^ ( -gamma / ( 1 - gamma ) );
    UHabit_Tmp_ = 1 / X_ * GYTrend_ * ( ( 1 - alphaS ) * PS_ / ( 1 + nuS ) + ( 1 - alphaT ) * PT_ / ( 1 + nuT ) + ( mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ - PT_ ) / ( 1 + nuM ) );
    kappaB_Tmp2_ = ( 1 - gamma ) * beta_ * GYTrend_ ^ ( -varsigma ) / ( 1 - beta_ * gamma * GYTrend_ ^ ( -varsigma ) * GYTrend_ ) * GYTrend_ ^ ( -gamma / ( 1 - gamma ) );
    Lambda1_ = C_ / X_ * GYTrend_ * ( 1 - U_OVER_C_STEADY ) * ( 1 + kappaB_Tmp2_ * UHabit_Tmp_ ) - UHabit_Tmp_;
    UHabit_ = C_ / X_ * GYTrend_ * ( 1 - U_OVER_C_STEADY );
    U_ = C_ * U_OVER_C_STEADY;
    kappaB_Tmp1_ = 1 - UHabit_ * kappaB_Tmp2_;

    WS_ = ( 1 - alphaS ) * PS_ / LS_;
    WT_ = ( 1 - alphaT ) * PT_ / LT_;
    WM_ = ( mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ - PT_ ) / LM_;

    snd_ = Omega_ * Psi_ * WM_ / PT_ * sdd_ / mu_ / tau * ( A_ ^ ( -tau ) - 1 );
    Check_ = -sfd_ + sqrt( sfd_ ^ 2 + A_ ^ tau * ( snd_ - 1 ) ) - Omega_ / Jd_;
    
    IsZero_ = AssertIsZero( Check_ );

    kappaB_ = U_ ^ ( -varsigma ) * kappaB_Tmp1_;
    kappaX_ = U_ ^ ( -varsigma ) * beta_ * GYTrend_ ^ ( -varsigma ) * UHabit_ / ( 1 - beta_ * gamma * GYTrend_ ^ ( -varsigma ) * GYTrend_ );

    LambdaS_ = kappaB_Tmp1_ / X_ * GYTrend_ * ( 1 - alphaS ) * PS_ / LS_ ^ ( 1 + nuS );
    LambdaT_ = kappaB_Tmp1_ / X_ * GYTrend_ * ( 1 - alphaT ) * PT_ / LT_ ^ ( 1 + nuT );
    LambdaM_ = kappaB_Tmp1_ / X_ * GYTrend_ * ( mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ - PT_ ) / LM_ ^ ( 1 + nuM );

    Y_ = A_;
    Td_ = 1 / Jd_;

    ZS_ = ( ( ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaS ) / GYTrend_ ) / GYTrend_ ) ^ alphaS * IS_ ^ alphaS * LS_ ^ ( 1 - alphaS ) ) ^ ( -1 / ( 1 - alphaS ) );
    ZT_ = ( ( ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaT ) / GYTrend_ ) / GYTrend_ ) ^ alphaT * IT_ ^ alphaT * LT_ ^ ( 1 - alphaT ) ) ^ ( -1 / ( 1 - alphaT ) );

    GC_ = GYTrend_;
    GI_ = GYTrend_;
    W_ = ( WS_ * LS_ + WT_ * LT_ + WM_ * LM_ ) / ( LS_ + LT_ + LM_ );
    GW_ = GYTrend_;
        
    Psi_STEADY = Psi_;
    Omega_STEADY = Omega_;
    Lambda1_STEADY = Lambda1_;
    LambdaS_STEADY = LambdaS_;
    LambdaT_STEADY = LambdaT_;
    LambdaM_STEADY = LambdaM_;
    ZS_STEADY = ZS_;
    ZT_STEADY = ZT_;
            
    @#include "InsertNewEndSteadyStateEquations.mod"
end;

@#define Deterministic = 0

shocks;
    @#if Deterministic
        var epsilon_TODO;
        periods 1;
        values 1;
    @#else
        @#include "InsertNewShockBlockLines.mod"
        
        @#for ShockProcessName in OtherShockProcessNames
            var epsilon_@{ShockProcessName} = 1;
        @#endfor
    @#endif
end;

options_.qz_criterium = 1 - 1e-8;

steady;
check;

@#if Estimation
    estimation( order = 1, datafile = 'EstimationData.xlsx', xls_sheet = Data, xls_range = A1:F336, plot_priors = 0, lik_init = 1, mh_replic = 0, mh_nblocks = 0, mode_check, prior_trunc = 0,
        mode_compute = 1, optim = ( 'Algorithm', 'sqp', 'Display', 'iter-detailed', 'DerivativeCheck', 'on', 'FinDiffType', 'central', 'MaxFunEvals', 1e12, 'MaxIter', 1e12, 'TolFun', 1e-16, 'TolCon', 1e-16, 'TolX', 1e-16 ),
        // mode_compute = 7, optim = ( 'Display', 'iter-detailed', 'MaxFunEvals', 1e12, 'MaxIter', 1e12, 'TolFun', 1e-16, 'TolX', 1e-16 ),
        // mode_compute = 9, optim = ( 'MaxFunEvals', 1e12, 'MaxIter', 1e12, 'TolFun', 1e-16, 'TolX', 1e-16 ),
        smoother, forecast = 400, kalman_algo = 1, keep_kalman_algo_if_singularity_is_detected, graph_format = none ) log_A level_Td log_Jd log_C log_IS log_IT log_LS log_LT log_LM log_Rs log_GPdollar log_GAs log_GC log_GI log_GW;
@#endif

@#ifndef dynareOBC
    save_params_and_steady_state( 'SteadyState.txt' );
@#endif

@#if Deterministic
    simul( periods = 10000, maxit = 1000000, tolf = 1e-8, tolx = 1e-8, stack_solve_algo = 7, solve_algo = 0 ); // endogenous_terminal_period
@#else
    @#if Estimation
        stoch_simul( order = 1, irf = 400, periods = 0, nofunctions, graph_format = none ) log_A level_Td log_Jd log_C log_IS log_IT log_LS log_LT log_LM; // k_order_solver, nocorr, nodisplay, nograph
    @#else
        stoch_simul( order = 2, irf = 400, periods = 10100, drop = 100, replic = 100, nofunctions, graph_format = none ) log_A level_Td log_Jd log_C log_IS log_IT log_LS log_LT log_LM; // k_order_solver, nocorr, nodisplay, nograph
    @#endif
@#endif
