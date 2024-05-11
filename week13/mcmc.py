# Serat
# Classwork, week13
# Date: 4/4/2024

import numpy as np
import matplotlib.pyplot as plt


def log_likelihood(m, b, x, y, sig):
    model = m * x + b
    chi_squared = np.sum(((y - model) / sig) ** 2)
    log_likelihood = -0.5 * chi_squared
    return log_likelihood


def log_prior(m, b):
    if (0 < m < 10) and (0 < b < 10):
        return 0.0 
    else:
        return -np.inf 


def log_posterior(m, b, x, y, sig):
    return log_likelihood(m, b, x, y, sig) + log_prior(m, b)


def mcmc(log_posterior, n_steps, step_size, x, y, sig, m, b):
    chain = np.zeros((n_steps, 2)) 
    accepted = np.zeros(n_steps, dtype=bool)
    
    m_current, b_current = m, b
    
    for i in range(n_steps):
        m_proposed = m_current + np.random.normal(scale=step_size)
        b_proposed = b_current + np.random.normal(scale=step_size)
        
        log_posterior_current = log_posterior(m_current, b_current, x, y, sig)
        log_posterior_proposed = log_posterior(m_proposed, b_proposed, x, y, sig)
        
        log_r = log_posterior_proposed - log_posterior_current
        if log_r > np.log(np.random.rand()): 
            m_current, b_current = m_proposed, b_proposed
            accepted[i] = True
        else:
            accepted[i] = False
        
        chain[i] = [m_current, b_current]
    
    return chain, accepted


def adjust_step_size(log_posterior, initial_step_size, target_acceptance_rate, x, y, sig):
    n_steps = 10000
    acceptance_rates = []

    for step_size in initial_step_size:
        chain, accepted = mcmc(log_posterior, n_steps, step_size, x, y, sig, m, b)
        acceptance_rate = np.mean(accepted)
        acceptance_rates.append(acceptance_rate)

    best_step_size = initial_step_size[np.argmin(np.abs(np.array(acceptance_rates) - target_acceptance_rate))]
    return best_step_size, acceptance_rates


def calculate_confidence_interval(samples, confidence_level=0.68):
    lower_quantile = 0.5 - confidence_level / 2
    upper_quantile = 0.5 + confidence_level / 2
    lower_bound = np.quantile(samples, lower_quantile)
    upper_bound = np.quantile(samples, upper_quantile)
    return lower_bound, upper_bound


def analyze_chain(chain):
    m_samples = chain[:, 0]
    b_samples = chain[:, 1]
    
    most_probable_m = np.median(m_samples)
    most_probable_b = np.median(b_samples)
    
    m_lower, m_upper = calculate_confidence_interval(m_samples)
    b_lower, b_upper = calculate_confidence_interval(b_samples)
    
    return most_probable_m, most_probable_b, m_lower, m_upper, b_lower, b_upper
    


if __name__ == '__main__':
    data = np.loadtxt('line.data')
    x = np.arange(10) + 0.5
    y = np.mean(data, axis=0)
    sig = np.var(data, axis=0, ddof=1)
    m = 1.0
    b = 1.0
    print("Log posterior for m =", m, "and b =", b, ":", log_posterior(m, b, x, y, sig))
    
    
    n_steps = 10000
    step_size = 0.1
    chain, accepted = mcmc(log_posterior, n_steps, step_size, x, y, sig, m, b)
    print("Acceptance rate:", np.mean(accepted))
    
    
    initial_step_size = [0.05, 0.1, 0.2, 0.3, 0.5]
    target_acceptance_rate = 0.3
    best_step_size, acceptance_rates = adjust_step_size(log_posterior, initial_step_size, target_acceptance_rate, x, y, sig)
    print("Best step size:", best_step_size)
    print("Acceptance rates for different step sizes:", acceptance_rates)
    
    most_probable_m, most_probable_b, m_lower, m_upper, b_lower, b_upper = analyze_chain(chain)
    print("Most probable value of m:", most_probable_m)
    print("Most probable value of b:", most_probable_b)
    print("68% confidence interval for m:", (m_lower, m_upper))
    print("68% confidence interval for b:", (b_lower, b_upper))
    

    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.plot(chain[:, 0], label='m chain', color='blue')
    plt.xlabel('Step Number')
    plt.ylabel('m')
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(chain[:, 1], label='b chain', color='red')
    plt.xlabel('Step Number')
    plt.ylabel('b')
    plt.legend()

    plt.tight_layout()
    plt.show()
    plt.savefig('mcmc.png')
    
    
    
    n_burn_in_steps = 5000 
    burned_in_chain = chain[n_burn_in_steps:]

    median_m = np.median(burned_in_chain[:, 0])
    median_b = np.median(burned_in_chain[:, 1])

    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.plot(chain[:, 0], label='m chain', color='blue')
    plt.axhline(y=median_m, color='green', linestyle='--', label='Median m')
    plt.xlabel('Step Number')
    plt.ylabel('m')
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(chain[:, 1], label='b chain', color='red')
    plt.axhline(y=median_b, color='orange', linestyle='--', label='Median b')
    plt.xlabel('Step Number')
    plt.ylabel('b')
    plt.legend()

    plt.tight_layout()
    plt.show()
    plt.savefig('burn.png')
    
    
    # 68% intervals for m and b
    m_lower, m_upper = calculate_confidence_interval(burned_in_chain[:, 0])
    b_lower, b_upper = calculate_confidence_interval(burned_in_chain[:, 1])

    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.hist(burned_in_chain[:, 0], bins=50, density=True, color='blue', alpha=0.7)
    plt.axvline(x=median_m, color='green', linestyle='--', label='Median m')
    plt.axvline(x=m_lower, color='red', linestyle='--', label='68% Interval')
    plt.axvline(x=m_upper, color='red', linestyle='--')
    plt.xlabel('m')
    plt.ylabel('Probability Density')
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.hist(burned_in_chain[:, 1], bins=50, density=True, color='red', alpha=0.7)
    plt.axvline(x=median_b, color='orange', linestyle='--', label='Median b')
    plt.axvline(x=b_lower, color='blue', linestyle='--', label='68% Interval')
    plt.axvline(x=b_upper, color='blue', linestyle='--')
    plt.xlabel('b')
    plt.ylabel('Probability Density')
    plt.legend()

    plt.tight_layout()
    plt.show()
    plt.savefig('final.png')

