#mu_hat: length n mean difference estimate
#sd_hat: length n standard deviation estimate of mu_hat
#n_effect: length n degree of freedom estimate of sd_hat
#alpha1, alpha2: alpha1- BH type level threshold, alpha2: p-value threshold
#delta: small quantity deducted from the empirical bayes estimate, should be an order much larger than 1/sqrt(n)
iteratedEB <- function(mu_hat, sd_hat, neffect, alpha1 = 0.05,  alpha2 = 0.01, delta = NULL){
    S1 = 1:length(mu_hat)
    n=length(S1)
    if(is.null(delta)){
        delta = sqrt(8/n)
    }
    S0 = c()
    tau2 = mean(mu_hat[S1]^2 - (1+delta)*sd_hat[S1]^2)
    if(tau2 < 0){
        tau2 = 0
    }
    while(!setequal(S0, S1)){
        S0 = S1
        t = mu_hat^2/(sd_hat^2+tau2)
        df = neffect*(tau2+sd_hat^2)^2/(sd_hat^4)
        pvalues = 1-pf(t, df1 = 1, df2 = df)
        pvalues_ordered = sort(pvalues)
        cutoffs = alpha1/n*(1:n)
        aa = which(pvalues_ordered < cutoffs)
        if(length(aa) > 0){
            k = max(aa)
            B1 = which(pvalues < pvalues_ordered[k] & pvalues < alpha2)
            if(length(B1) > 0){
                S1 = setdiff(S0, B1)
            }
            tau2 = mean(mu_hat[S1]^2-(1+delta)*sd_hat[S1]^2)
            if(tau2 < 0){
                tau2 = 0
            }
        }
    }
    return(list(sqrt(tau2), S1))
}

#mu_hat: length n mean difference estimate
#sum_vars: length n sum of variances
#m: m = m_0 - 1 if m_0 is the number of paired samples in the paired test
#   and m = m_1+m_2 - 2 if m_1, m_2 are the samples from the control and experiment groups respectively in two sample test
#   var(mu_hat) = E[sum_vars/m] 
#n_effect: length n degree of freedom estimate of var(mu_hat)
#p: proportion removed from the tails
#delta: small quantity deducted from the empirical bayes estimate, should be an order much larger than 1/sqrt(n)
truncatedMLE <- function(mu_hat, sum_vars,   m, n_effect, p, tol=10^(-3)){
    #find the minimum of the negative log likelihood by blockwise descent: tau ->mu->tau...
    tau_search <- function(tau, sigma, n_effect){
        alpha = delta0/sqrt(tau+sigma/n_effect)
        Phi = pnorm(alpha)
        a = 2*Phi - 1
        r = sum(log(a)) + sum(log(tau+sigma/n_effect))/2+sum(mean_x^2/(2*tau+2*sigma/n_effect))
        r
    }
    sigma_search <- function(sigmai, tau, sumV_xi, mean_xi, n_effecti, m){
        alpha = delta0/sqrt(tau+sigmai/n_effecti)
        Phi = pnorm(alpha)
        a = 2*Phi - 1
        r = log(a) + log(tau+sigmai/n_effecti)/2+sumV_xi/(2*sigmai)+m*log(sigmai)/2+mean_xi^2/(2*(tau+sigmai/n_effecti))
        r
    }
    delta0 = quantile(abs(mu_hat), p)
    S = which(abs(mu_hat) <= delta0)
    n = length(S)
    mean_x = mu_hat[S]
    sumV_x = sum_vars[S]
    n_effect = n_effect[S]
    tau0 = 0
    tau1 = sum(mean_x^2)/n
    sigma = sumV_x/m[S]
    while(abs(tau1 - tau0) > tol){
        tau0 = tau1
        for(i in 1:length(mean_x)){
            sigma[i] = optimize(sigma_search, interval = c(0,10), tau = tau0, sumV_xi = sumV_x[i],
            mean_xi = mean_x[i], n_effecti = n_effect[i], m = m)[[1]]
        }
        tau1 = optimize(tau_search, interval = c(0,10), sigma = sigma, n_effect = n_effect)[[1]]
    }
    return(sqrt(tau0))
}

Az_scores_ITEB <- function(mu, sigma, n_effect, alpha1, alpha2, delta = NULL){
    Nsub = length(mu)
    tau = iteratedEB(mu_hat = mu, sd_hat = sigma, n_effect = n_effect,threshold =alpha1, threshold2 = alpha2, delta = delta)[[1]]
    az = mu/sqrt(sigma^2+tau^2)
    az2 = qnorm(pt(az, df = neffect*(tau^2+sigma^2)^2/sigma^4, log.p = T),log.p = T)
    return(list(az = az2, tau = tau))
}

Az_scores_TMLE <- function(mu,sigma, sum_vars, m, n_effect, p = 0.8){
    Nsub = length(mu)
    tau = truncatedMLE(mu_hat = mu , sum_vars = sum_vars, m = m, n_effect=n_effect, p =p, tol=10^(-3))
    az = mu/sqrt(sigma^2+tau^2)
    az2 = qnorm(pt(az, df = (m)*(tau^2+sigma^2)^2/sigma^4, log.p = T),log.p = T)
    return(list(az = az2, tau = tau))
}
