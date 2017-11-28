function [ Jd_STEADY, lambda, alphaS ] = CalibrateModel( R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 )

    In = [ 10.6592295099414; 0.519377524696136; 0.364331847328728 ];
    
    Weights = [ 5; 12; 1; 0.1; 10; 0.0001; 0.0001; 0.0001 ];
    InTarget = [ 3; 1; 0.3 ];
    
    if false
    
        LB = [ 2; 0.1; 0 ];
        UB = [ Inf; Inf; 0.5 ];

        WarningState = warning( 'off', 'all' );

        if false
            CMAESOptions = CMAESMinimisation;
            CMAESOptions.MaxIter = Inf;
            CMAESOptions.LBounds = LB;
            CMAESOptions.UBounds = UB;
            CMAESOptions.PopSize = 80;
            CMAESOptions.TolX = eps;
            CMAESOptions.TolFun = eps;
            CMAESOptions.TolHistFun = eps;

            Sigma = ( UB - LB ) / 2;
            Sigma( ~isfinite( Sigma ) ) = 1;

            [ ~, ~, ~, ~, ~, BestEver ] = CMAESMinimisation( @( InMat, ~, ~ ) CMAESObjective( InMat, InTarget, Weights, R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 ),...
                In, Sigma, 1, CMAESOptions );
            In = BestEver.x;
        end

        if true
            In = fmincon( @( In_ ) Objective( In_, InTarget, Weights, R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 ), ...
                In, [], [], [], [], LB, UB, ...
                @( In_ ) Constraint( In_, R_STEADY, tau, eta, A_STEADY, LM_STEADY, GAs_STEADY, varsigma ), ...
                optimoptions( @fmincon, 'Algorithm', 'sqp', 'Display', 'iter-detailed', 'CheckGradients', false, 'FiniteDifferenceType', 'central', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12 ) ); %#ok<*UNRCH>
        end

        warning( WarningState );
    
    end
    
    Jd_STEADY = In( 1 );
    lambda = In( 2 );
    alphaS = In( 3 );
    
    [ FinalNorm, ~, FinalResiduals ] = Objective( In, InTarget, Weights, R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 );
    
    assert( FinalNorm < 0.86 );
    assert( all( abs( FinalResiduals( 1 : ( end - 3 ) ) ) < 0.7 ) );

    if false

        disp( 'Jd_STEADY, lambda, alphaS:' );
        disp( [ Jd_STEADY, lambda, alphaS ] );
        disp( 'Final residuals:' );
        disp( FinalResiduals( 1 : ( end - 4 ) ) );
        disp( 'Final norm:' );
        disp( FinalNorm );

        DisplayResults( In, R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 );
    
    end
    
end

function [ Minimand, DMinimand, Residuals ] = Objective( In, InTarget, Weights, R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 )

    Jd_STEADY = In( 1 );
    lambda = In( 2 );
    alphaS = In( 3 );

    GAs_ = GAs_STEADY;
    beta_ = GAs_STEADY ^ ( varsigma / ( 1 - alphaS ) ) / R_STEADY;

    Jd_ = Jd_STEADY;
    A_ = A_STEADY;
    LM_ = LM_STEADY;

    GYTrend_ = GAs_ ^ ( 1 / ( 1 - alphaS ) );

    Xi_LEAD_ = beta_ * GYTrend_ ^ ( -varsigma );

    QS_ = 1 / ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 - PhiS2 * ( GYTrend_ - 1 ) * GYTrend_ + Xi_LEAD_ * PhiS2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2 );
    QT_ = 1 / ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 - PhiT2 * ( GYTrend_ - 1 ) * GYTrend_ + Xi_LEAD_ * PhiT2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2 );

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
    
    % KS_ = ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaS ) / GYTrend_ ) * IS_;
    KT_ = ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaT ) / GYTrend_ ) * IT_;

    % SRS_ = alphaS * PS_ * GYTrend_ / ( ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaS ) / GYTrend_ ) ) / IS_;
    SRT_ = alphaT * PT_ * GYTrend_ / ( ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaT ) / GYTrend_ ) ) / IT_;
    
    HC_Income_ = SRT_ * KT_ / GYTrend_ * HC_SHARE / alphaT;
    LT_Tmp_ = ( HC_Income_ + ( 1 - alphaT ) * PT_ ) / ( 1 - alphaS ) / PS_ / WTEFF_OVER_WS_STEADY;
    LT_ = LT_Tmp_ / ( 1 + LT_Tmp_ ) * ( 1 - LM_ );
    LS_ = 1 - LT_ - LM_;
    WM_OVER_WS_ = ( mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ - PT_ ) / ( 1 - alphaS ) / PS_ / LM_ * LS_;

    WS_ = ( 1 - alphaS ) * PS_ / LS_;
    WT_ = ( 1 - alphaT ) * PT_ / LT_;
    WM_ = ( mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ - PT_ ) / LM_;

    C_ = A_ - IS_ - IT_;
    
    EntryCosts_ = Psi_ * WM_ * Jd_;
    GDP_ = A_ + PT_;
    
    IHC_ = IT_ * HC_SHARE / alphaT;

    Residuals = [ ( WM_OVER_WS_ - WM_OVER_WS_TARGET ) / WM_OVER_WS_TARGET;
                  ( ( C_ + IHC_ ) / GDP_ - CONSUMPTION_SHARE_TARGET ) / CONSUMPTION_SHARE_TARGET;
                  ( ( PT_ + IT_ ) / GDP_ - RESEARCH_SHARE_TARGET ) / RESEARCH_SHARE_TARGET;
                  ( ( EntryCosts_ + PT_ ) / ( GDP_ + EntryCosts_ ) - FIXED_COST_SHARE_TARGET ) / FIXED_COST_SHARE_TARGET;
                  ( ( WS_ * LS_ + WT_ * LT_ + WM_ * LM_ + HC_Income_ ) / GDP_ - LABOUR_SHARE_TARGET ) / LABOUR_SHARE_TARGET;
                  ( In - InTarget ) ];
    Minimand = sum( ( ( Weights .* Residuals ) .^ 2 ) .^ ( 1 / 2 ) ) + 1e8 * max( 0, 2 - WM_OVER_WS_ );
    
    if nargout > 1
        Step = sqrt( eps );
        DMinimand = zeros( 3, 1 );
        for j = 1 : 3
            In( j ) = In( j ) + Step * 1i;
            DMinimand( j ) = imag( Objective( In, InTarget, Weights, R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 ) ) / Step;
            In( j ) = In( j ) - Step * 1i;
        end
    end

end

function DisplayResults( In, R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 )

    Jd_STEADY = In( 1 );
    lambda = In( 2 );
    alphaS = In( 3 );

    GAs_ = GAs_STEADY;
    beta_ = GAs_STEADY ^ ( varsigma / ( 1 - alphaS ) ) / R_STEADY;

    Jd_ = Jd_STEADY;
    A_ = A_STEADY;
    LM_ = LM_STEADY;

    GYTrend_ = GAs_ ^ ( 1 / ( 1 - alphaS ) );

    Xi_LEAD_ = beta_ * GYTrend_ ^ ( -varsigma );

    QS_ = 1 / ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 - PhiS2 * ( GYTrend_ - 1 ) * GYTrend_ + Xi_LEAD_ * PhiS2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2 );
    QT_ = 1 / ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 - PhiT2 * ( GYTrend_ - 1 ) * GYTrend_ + Xi_LEAD_ * PhiT2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2 );

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
    
    KS_ = ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaS ) / GYTrend_ ) * IS_;
    KT_ = ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaT ) / GYTrend_ ) * IT_;

    % SRS_ = alphaS * PS_ * GYTrend_ / ( ( 1 - PhiS2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaS ) / GYTrend_ ) ) / IS_;
    SRT_ = alphaT * PT_ * GYTrend_ / ( ( 1 - PhiT2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) / ( 1 - ( 1 - deltaT ) / GYTrend_ ) ) / IT_;
    
    HC_Income_ = SRT_ * KT_ / GYTrend_ * HC_SHARE / alphaT;
    LT_Tmp_ = ( HC_Income_ + ( 1 - alphaT ) * PT_ ) / ( 1 - alphaS ) / PS_ / WTEFF_OVER_WS_STEADY;
    LT_ = LT_Tmp_ / ( 1 + LT_Tmp_ ) * ( 1 - LM_ );
    LS_ = 1 - LT_ - LM_;
    WM_OVER_WS_ = ( mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ - PT_ ) / ( 1 - alphaS ) / PS_ / LM_ * LS_;

    WS_ = ( 1 - alphaS ) * PS_ / LS_;
    WT_ = ( 1 - alphaT ) * PT_ / LT_;
    WM_ = ( mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ - PT_ ) / LM_;

    C_ = A_ - IS_ - IT_;
    
    EntryCosts_ = Psi_ * WM_ * Jd_;
    GDP_ = A_ + PT_;
    
    IHC_ = IT_ * HC_SHARE / alphaT;

    disp( 'LT_:' );
    disp( LT_ );
    disp( 'Effective delta:' );
    disp( ( deltaS * KS_ + deltaT * KT_ * ( 1 - HC_SHARE / alphaT ) ) / ( KS_ + KT_ * ( 1 - HC_SHARE / alphaT ) ) );
    disp( 'WM_OVER_WS_, WM_OVER_WS_ TARGET:' );
    disp( [ WM_OVER_WS_, WM_OVER_WS_TARGET ] );
    disp( 'CONSUMPTION_SHARE, CONSUMPTION_SHARE_TARGET:' );
    disp( [ ( C_ + IHC_ ) / GDP_, CONSUMPTION_SHARE_TARGET ] );
    disp( 'RESEARCH_SHARE, RESEARCH_SHARE_TARGET:' );
    disp( [ ( PT_ + IT_ ) / GDP_, RESEARCH_SHARE_TARGET ] );
    disp( 'FIXED_COST_SHARE, FIXED_COST_SHARE_TARGET:' );
    disp( [ ( EntryCosts_ + PT_ ) / ( GDP_ + EntryCosts_ ), FIXED_COST_SHARE_TARGET ] );
    disp( 'LABOUR_SHARE, LABOUR_SHARE_TARGET:' );
    disp( [ ( WS_ * LS_ + WT_ * LT_ + WM_ * LM_ + HC_Income_ ) / GDP_, LABOUR_SHARE_TARGET ] );

end

function [ Ineq, Eq, DIneq, DEq ] = Constraint( In, R_STEADY, tau, eta, A_STEADY, LM_STEADY, GAs_STEADY, varsigma )

    Ineq = [];
    DIneq = [];

    Jd_STEADY = In( 1 );
    lambda = In( 2 );
    alphaS = In( 3 );

    GAs_ = GAs_STEADY;
    beta_ = GAs_STEADY ^ ( varsigma / ( 1 - alphaS ) ) / R_STEADY;

    Jd_ = Jd_STEADY;
    A_ = A_STEADY;
    LM_ = LM_STEADY;

    GYTrend_ = GAs_ ^ ( 1 / ( 1 - alphaS ) );

    Xi_LEAD_ = beta_ * GYTrend_ ^ ( -varsigma );

    Psi_ = LM_ / Jd_;

    Omega_ = A_ ^ tau * ( GAs_ ^ tau - 1 ) * Jd_ / ( 1 - A_ ^ tau * GAs_ ^ tau );
    mu_ = lambda * eta * Jd_ / ( Jd_ - ( 1 - eta ) );
    omegad_ = Jd_ * ( 1 - eta ) / ( ( Jd_ - ( 1 - eta ) ) ^ 2 * ( 1 + mu_ ) );
    sdd_ = 1 - omegad_ / ( 1 + omegad_ ) * ( lambda - mu_ ) * ( mu_ - eta * lambda ) / ( lambda * ( 1 - eta ) * mu_ );
    sfd_ = 1 - 1 / 2 * ( 1 + sdd_ / mu_ / tau ) * ( 1 - A_ ^ tau );
    PT_Tmp_ = ( 1 + ( ( Omega_ / Jd_ + sfd_ ) ^ 2 - sfd_ ^ 2 ) / A_ ^ tau ) / Omega_ / Psi_ * LM_ / sdd_ * mu_ * tau / ( A_ ^ ( -tau ) - 1 );
    PT_ = mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ / ( 1 + PT_Tmp_ );

    WM_ = ( mu_ / ( 1 + mu_ ) * Xi_LEAD_ * A_ * GYTrend_ - PT_ ) / LM_;

    snd_ = Omega_ * Psi_ * WM_ / PT_ * sdd_ / mu_ / tau * ( A_ ^ ( -tau ) - 1 );
    Eq = -sfd_ + sqrt( sfd_ ^ 2 + A_ ^ tau * ( snd_ - 1 ) ) - Omega_ / Jd_;

    if nargout > 2
        Step = sqrt( eps );
        DEq = zeros( 3, 1 );
        for j = 1 : 3
            In( j ) = In( j ) + Step * 1i;
            [ ~, TEq ] = Constraint( In, R_STEADY, tau, eta, A_STEADY, LM_STEADY, GAs_STEADY, varsigma );
            DEq( j ) = imag( TEq ) / Step;
            In( j ) = In( j ) - Step * 1i;
        end
    end

end

function [ MinimandVec, PersistentState ] = CMAESObjective( InMat, InTarget, Weights, R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 ) %#ok<*DEFNU>
    PersistentState = 1;
    N = size( InMat, 2 );
    MinimandVec = zeros( 1, N );
    for n = 1 : N
        [ ~, Eq ] = Constraint( InMat( :, n ), R_STEADY, tau, eta, A_STEADY, LM_STEADY, GAs_STEADY, varsigma );
        MinimandVec( n ) = 1e8 * abs( Eq ) + Objective( InMat( :, n ), InTarget, Weights, R_STEADY, tau, eta, alphaT, HC_SHARE, deltaS, deltaT, A_STEADY, LM_STEADY, WTEFF_OVER_WS_STEADY, WM_OVER_WS_TARGET, RESEARCH_SHARE_TARGET, CONSUMPTION_SHARE_TARGET, FIXED_COST_SHARE_TARGET, LABOUR_SHARE_TARGET, GAs_STEADY, varsigma, PhiS2, PhiT2 );
    end
end
