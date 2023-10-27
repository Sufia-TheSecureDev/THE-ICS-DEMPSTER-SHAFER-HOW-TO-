# function for step 2 calculation
def calculate_probability_of_mass(weights, degree_of_belief):
    probability_mass = []
    m_tilde_H = []  # Initialize m_tilde_H as a list
    mH = []
    m_bar_H = []

    for i in range(len(weights)):
        p = []
        for j in range(3):
            p.append(round(weights[i] * degree_of_belief[i][j], 4))
        p.append(1 - weights[i])
        p.append(round(weights[i] * degree_of_belief[i][3], 4))
        p.insert(3, round(1 - weights[i] + p[-1], 4))
        probability_mass.append(p)
        # Calculate m_tilde_H for each attribute εi
        m_tilde_H_i = weights[i] * (1 - sum(degree_of_belief[i]))
        m_tilde_H.append(round(weights[i] * degree_of_belief[i][3], 4))

        # Calculate remaining probability mass mH for each basic attribute εi
        mH_i = 1 - sum(p[j] for j in range(3))
        mH.append(mH_i)

        # Decompose mH into ¯mH for each basic attribute εi
        m_bar_H_i = 1 - weights[i]
        m_bar_H.append(m_bar_H_i)

    return probability_mass, mH, m_bar_H, m_tilde_H


# functon for constant calculation
def calculate_constant_k(probability_masses):
    k = [None]

    for i in range(1, len(probability_masses)):
        total_k = 0
        for t in range(3):  # t represents indices for N
            for j in range(3):  # j represents indices for N
                if t == j:
                    continue
                total_k += probability_masses[i - 1][t] * probability_masses[i][j]
        total_k = 1 - total_k
        k_i = 1 / total_k
        k.append(round(k_i, 5))  # Round to four decimal places as required

    return k


# step 3 calculation
def calculate_probability_mass_aggregation(probability_masses, m_H, m_bar_H, m_tilde_H, constant):
    PoMa = [None]  # Initialize with a None placeholder

    updated_probability_masses = probability_masses.copy()
    updated_m_H = m_H.copy()
    updated_m_bar_H = m_bar_H.copy()
    updated_m_tilde_H = m_tilde_H.copy()

    for i in range(1, len(updated_probability_masses)):
        pma = []

        for j in range(3):
            p = constant[i] * (updated_probability_masses[i - 1][j] * updated_probability_masses[i][j] +
                               updated_m_H[i - 1] * updated_probability_masses[i][j] +
                               updated_probability_masses[i - 1][j] * updated_m_H[i])
            pma.append(round(p, 4))

        m_bar = constant[i] * (updated_m_bar_H[i - 1] * updated_m_bar_H[i])
        pma.append(round(m_bar, 4))

        m_tilde = constant[i] * (updated_m_tilde_H[i - 1] * updated_m_tilde_H[i] +
                                 updated_m_bar_H[i - 1] * updated_m_tilde_H[i] +
                                 updated_m_tilde_H[i - 1] * updated_m_bar_H[i])
        pma.append(round(m_tilde, 4))

        total_m = round(m_bar + m_tilde, 4)
        pma.insert(3, total_m)

        pma = [round(val, 4) for val in pma]

        PoMa.append(pma)

        # Update the values in the updated variables
        updated_probability_masses[i] = pma
        updated_m_H[i] = m_bar + m_tilde
        updated_m_bar_H[i] = m_bar
        updated_m_tilde_H[i] = m_tilde

    return PoMa


# step 4 calculation
def calculate_combined_degrees_of_belief(PoMa):
    combined_degrees_of_belief = []

    m_n_L = PoMa[-1][:-3]  # Get the last row of PoMa without the last element
    m_bar = PoMa[-1][-2]  # Get the last element of the last row of PoMa
    m_tilda = PoMa[-1][-1]  # Get the last element of the last row of PoMa

    for n in range(len(m_n_L)):
        # Calculate β_n
        beta_n = m_n_L[n] / (1 - m_bar)
        combined_degrees_of_belief.append(round(beta_n, 4))

    # Calculate β_H
    beta_H = m_tilda / (1 - m_bar)
    combined_degrees_of_belief.append(round(beta_H, 4))

    return combined_degrees_of_belief


# step 5 calculation for estimate utilities
def estimate_utilities(assessment_grades):
    """
    Estimates utilities for different assessment grades. This function returns a dictionary where the keys are the
    assessment grades and the values are the estimated utilities.

    :param assessment_grades: A dictionary where keys are assessment grades (H1, H2, H3, etc.) and values are the
                              corresponding utility values.
    :return: A dictionary of utilities for each assessment grade.
    """
    utilities = {}
    for grade, utility in assessment_grades.items():
        utilities[grade] = utility
    return utilities


# step 5 calculation for expected utility
def calculate_expected_utility(assessment_weights, utilities):
    """
    Calculate the expected utility for a complete assessment using the provided assessment weights and utilities.

    :param assessment_weights: A list of assessment weights (β1, β2, β3, etc.).
    :param utilities: A dictionary of utilities for each assessment grade.
    :return: The expected utility for the complete assessment.
    """
    expected_utility = sum(assessment_weights[i] * utilities[f'H{i + 1}'] for i in range(len(assessment_weights)))
    return expected_utility


# step 5 calculation for utility interval
def calculate_utility_interval(assessment_weights, utilities, unassigned_belief):
    """
    Calculate the utility interval for incomplete assessments.

    :param assessment_weights: A list of assessment weights (β1, β2, β3, etc.).
    :param utilities: A dictionary of utilities for each assessment grade.
    :param unassigned_belief: The unassigned belief degree (βH).
    :return: A tuple (u_min, u_max, u_avg) representing the utility interval.
    """
    u_max = sum(assessment_weights[i] * utilities[f'H{i + 1}'] for i in range(len(assessment_weights) - 1)) + \
            (assessment_weights[-1] + unassigned_belief) * utilities[f'H{len(assessment_weights)}']

    u_min = (assessment_weights[0] + unassigned_belief) * utilities['H1'] + \
            sum(assessment_weights[i] * utilities[f'H{i + 1}'] for i in range(1, len(assessment_weights)))

    u_avg = (u_max + u_min) / 2

    return u_min, u_max, u_avg


# output format
def print_results(data_type, weights, degree_of_belief, probability_masses, mH, m_bar_H, k, PoMa, combined_degrees):
    print(f"{data_type} Data Quality:")
    print("Weight =", weights)
    print("Belief =", degree_of_belief)
    print("Probability Mass =", probability_masses)
    print("Constant =", k)
    print("Probability Mass (aggregation) =", PoMa)
    print("Combined Degrees of Belief =", combined_degrees)


# Define intrinsic weights and degree of belief
intrinsic_weights = [0.35, 0.65]
intrinsic_degree_of_belief = [[0.4, 0.5, 0.0, 0.1], [0.1, 0.75, 0.15, 0.0]]

# Calculate intrinsic data quality
intrinsic_probability_masses, m_H_intrinsic, m_bar_H_intrinsic, m_tilde_H_intrinsic = calculate_probability_of_mass(
    intrinsic_weights, intrinsic_degree_of_belief)
k_intrinsic = calculate_constant_k(intrinsic_probability_masses)
PoMa_intrinsic = calculate_probability_mass_aggregation(intrinsic_probability_masses, m_H_intrinsic,
                                                        m_bar_H_intrinsic, m_tilde_H_intrinsic, k_intrinsic)
combined_degrees_intrinsic = calculate_combined_degrees_of_belief(PoMa_intrinsic)

# Define contextual weights and degree of belief
contextual_weights = [0.45, 0.25, 0.3]
contextual_degree_of_belief = [[0.6, 0.2, 0.05, 0.15], [0.25, 0.45, 0.3, 0.0], [0.55, 0.35, 0.0, 0.1]]

# Calculate contextual data quality
contextual_probability_masses, m_H_contextual, m_bar_H_contextual, m_tilde_H_contextual = calculate_probability_of_mass(
    contextual_weights, contextual_degree_of_belief)
k_contextual = [None, 1.07174, 1.0797]

PoMa_contextual = calculate_probability_mass_aggregation(contextual_probability_masses, m_H_contextual,
                                                         m_bar_H_contextual, m_tilde_H_contextual, k_contextual)
combined_degrees_contextual = calculate_combined_degrees_of_belief(PoMa_contextual)

# Define data quality's weight and degree of belief
data_weights = [0.6, 0.4]
data_degree_of_belief = [combined_degrees_intrinsic, combined_degrees_contextual]

# Calculate  data quality
data_probability_masses, m_H_data, m_bar_H_data, m_tilde_H_data = calculate_probability_of_mass(
    data_weights, data_degree_of_belief)
k_data = [None, 1.16499]

PoMa_data = calculate_probability_mass_aggregation(data_probability_masses, m_H_data,
                                                   m_bar_H_data, m_tilde_H_data, k_data)
combined_degrees_data = calculate_combined_degrees_of_belief(PoMa_data)

# Print the results for intrinsic data quality
print_results("Intrinsic", intrinsic_weights, intrinsic_degree_of_belief, intrinsic_probability_masses,
              m_H_intrinsic, m_bar_H_intrinsic, k_intrinsic, PoMa_intrinsic, combined_degrees_intrinsic)

# Print the results for contextual data quality
print_results("Contextual", contextual_weights, contextual_degree_of_belief, contextual_probability_masses,
              m_H_contextual, m_bar_H_contextual, k_contextual, PoMa_contextual, combined_degrees_contextual)

print_results("Final", data_weights, data_degree_of_belief, data_probability_masses,
              m_H_data, m_bar_H_data, k_data, PoMa_data, combined_degrees_data)

print("----------------------------------------------")

# Define utilities for assessment grades (H1, H2, H3, etc.)
assessment_utilities = {'H1': 0, 'H2': 0.5, 'H3': 1}

# Calculate expected utility for complete assessments
complete_assessment_weights = combined_degrees_data[:3]

# Define unassigned belief (βH) for incomplete assessments
unassigned_belief = combined_degrees_data[-1]

expected_utility = calculate_expected_utility(complete_assessment_weights, assessment_utilities)
print(f"Expected Utility for Complete Assessment: {expected_utility:.4f}")

# Calculate utility interval for incomplete assessments
incomplete_assessment_weights = complete_assessment_weights  # Assuming the same weights
u_min, u_max, u_avg = calculate_utility_interval(incomplete_assessment_weights, assessment_utilities, unassigned_belief)
print(f"Utility Interval (u_min): {u_min:.4f}")
print(f"Utility Interval (u_max): {u_max:.4f}")
print(f"Utility Interval (u_avg): {u_avg:.4f}")

# Define your thresholds
Threshold_Good = 0.8  # Adjust this value as needed
Threshold_Poor = 0.4  # Adjust this value as needed

# Make a decision based on the thresholds
decision = ""
if u_avg >= Threshold_Good:
    decision = "Good"
elif u_avg < Threshold_Poor:
    decision = "Poor"
else:
    decision = "Average"

# Print the decision
print(f"-------------------\nData Quality Decision: {decision}")