######------ metropolis_hastings generates random samples------######
import random

def metropolis_hastings(P, proposal_distribution, initial_value, iterations):
    samples = [initial_value]
    current_value = initial_value

    for _ in range(iterations):
        # Propose a new sample from the proposal distribution
        proposed_value = proposal_distribution(current_value)
        # Ensure proposed value is greater than 0
        if proposed_value <= 0:
            continue
        # Calculate acceptance probability
        acceptance_prob = min(1, P(proposed_value) / P(current_value))
        # Accept or reject the proposed sample
        if random.random() < acceptance_prob:
            current_value = proposed_value
            samples.append(current_value)
    return samples
