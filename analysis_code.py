"""
ECON 5193: Applied Development Economics
Problem Set 2: Cluster Randomization and Spillover Effects
Empirical Analysis (Questions 4-7)
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Load data
df = pd.read_csv('/mnt/user-data/uploads/water_program_data.csv')

print("="*80)
print("DATASET OVERVIEW")
print("="*80)
print(f"Number of schools: {len(df)}")
print(f"Treatment schools: {df['treatment'].sum()}")
print(f"Control schools: {(df['treatment']==0).sum()}")
print(f"\nVariables: {list(df.columns)}")
print(f"\nDescriptive statistics:")
print(df.describe().round(4))

# ============================================================================
# QUESTION 4: Balance Checks
# ============================================================================
print("\n" + "="*80)
print("QUESTION 4: BALANCE CHECKS")
print("="*80)

balance_vars = ['num_students', 'num_teachers', 'baseline_absence_rate', 
                'pct_female_students', 'has_latrine']
var_labels = {
    'num_students': 'Number of Students',
    'num_teachers': 'Number of Teachers',
    'baseline_absence_rate': 'Baseline Absence Rate',
    'pct_female_students': '% Female Students',
    'has_latrine': 'Has Functioning Latrine'
}

treat = df[df['treatment'] == 1]
control = df[df['treatment'] == 0]

balance_rows = []
for var in balance_vars:
    t_mean = treat[var].mean()
    c_mean = control[var].mean()
    diff = t_mean - c_mean
    # t-test for difference in means
    t_stat, p_val = stats.ttest_ind(treat[var], control[var])
    balance_rows.append({
        'Variable': var_labels[var],
        'Treatment Mean': round(t_mean, 3),
        'Control Mean': round(c_mean, 3),
        'Difference': round(diff, 3),
        'p-value': round(p_val, 3)
    })

balance_df = pd.DataFrame(balance_rows)
print("\n(a) & (b) Balance Table:")
print("-"*80)
print(balance_df.to_string(index=False))
print("-"*80)

print(f"\nN (Treatment): {len(treat)}")
print(f"N (Control): {len(control)}")

# ============================================================================
# QUESTION 5: Direct Treatment Effects
# ============================================================================
print("\n" + "="*80)
print("QUESTION 5: DIRECT TREATMENT EFFECTS")
print("="*80)

# (a) Simple difference in means
t_mean_end = treat['endline_absence_rate'].mean()
c_mean_end = control['endline_absence_rate'].mean()
diff_means = t_mean_end - c_mean_end

# SE of difference in means
se_diff = np.sqrt(treat['endline_absence_rate'].var()/len(treat) + 
                  control['endline_absence_rate'].var()/len(control))

t_stat_simple, p_val_simple = stats.ttest_ind(treat['endline_absence_rate'], 
                                               control['endline_absence_rate'])

print(f"\n(a) Simple Difference in Means:")
print(f"  Treatment mean (endline): {t_mean_end:.4f}")
print(f"  Control mean (endline):   {c_mean_end:.4f}")
print(f"  Difference:               {diff_means:.4f}")
print(f"  Standard Error:           {se_diff:.4f}")
print(f"  t-statistic:              {t_stat_simple:.4f}")
print(f"  p-value:                  {p_val_simple:.4f}")

# (b) Regression with controls
df['district_jinja'] = (df['district'] == 2).astype(int)
model_b = smf.ols('endline_absence_rate ~ treatment + baseline_absence_rate + district_jinja', 
                   data=df).fit()

print(f"\n(b) OLS Regression with Controls:")
print(model_b.summary().tables[1])

# (c) Robust standard errors
model_c = smf.ols('endline_absence_rate ~ treatment + baseline_absence_rate + district_jinja', 
                   data=df).fit(cov_type='HC1')

print(f"\n(c) OLS Regression with Robust (HC1) Standard Errors:")
print(model_c.summary().tables[1])

# Compare
print("\nComparison of Standard Errors:")
print(f"  {'Variable':<25} {'Default SE':>12} {'Robust SE':>12}")
for var_name in model_b.params.index:
    print(f"  {var_name:<25} {model_b.bse[var_name]:>12.4f} {model_c.bse[var_name]:>12.4f}")

# (d) Substantive interpretation
avg_students = df['num_students'].mean()
treatment_effect = model_c.params['treatment']
print(f"\n(d) Substantive Interpretation:")
print(f"  Treatment effect (coefficient): {treatment_effect:.4f}")
print(f"  Average students per school: {avg_students:.1f}")
print(f"  Student-days gained per school per month: {abs(treatment_effect) * avg_students * 20:.1f}")
print(f"    (using avg school of {avg_students:.0f} students × 20 school days × {abs(treatment_effect):.4f})")
print(f"  Student-days gained per school per year (10 months): {abs(treatment_effect) * avg_students * 200:.1f}")

# ============================================================================
# QUESTION 6: Spillover Effects
# ============================================================================
print("\n" + "="*80)
print("QUESTION 6: SPILLOVER EFFECTS")
print("="*80)

control_df = df[df['treatment'] == 0].copy()

# (a) Scatter plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(control_df['distance_nearest_treat'], control_df['endline_absence_rate'], 
           alpha=0.6, color='steelblue', edgecolors='white', s=60)

# Fitted line
slope, intercept, r_value, p_value_slope, std_err_slope = stats.linregress(
    control_df['distance_nearest_treat'], control_df['endline_absence_rate'])
x_line = np.linspace(control_df['distance_nearest_treat'].min(), 
                      control_df['distance_nearest_treat'].max(), 100)
ax.plot(x_line, intercept + slope * x_line, color='red', linewidth=2, 
        label=f'Fitted line (slope={slope:.4f})')

ax.set_xlabel('Distance to Nearest Treatment School (km)', fontsize=12)
ax.set_ylabel('Endline Absence Rate', fontsize=12)
ax.set_title('Control Schools: Absence Rate vs. Distance to Nearest Treatment School', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('/home/claude/scatter_spillover.png', dpi=150)
plt.close()
print(f"\n(a) Scatter plot saved.")
print(f"  Simple correlation: r = {r_value:.4f}")
print(f"  Slope of fitted line: {slope:.4f}")
print(f"  Pattern: {'Positive' if slope > 0 else 'Negative'} relationship")

# (b) Spillover regression
control_df['district_jinja'] = (control_df['district'] == 2).astype(int)
model_spill = smf.ols('endline_absence_rate ~ distance_nearest_treat + baseline_absence_rate + district_jinja', 
                       data=control_df).fit(cov_type='HC1')

print(f"\n(b) Spillover Regression (Control Schools Only, Robust SE):")
print(model_spill.summary().tables[1])
print(f"\n  Interpretation: A 1 km increase in distance to nearest treatment school")
print(f"  is associated with a {model_spill.params['distance_nearest_treat']:.4f} change in absence rate.")

# (c) Predicted difference
coef_dist = model_spill.params['distance_nearest_treat']
diff_1_vs_10 = coef_dist * (10 - 1)
print(f"\n(c) Predicted Difference (1 km vs 10 km):")
print(f"  Coefficient on distance: {coef_dist:.4f}")
print(f"  Predicted difference = {coef_dist:.4f} × (10 - 1) = {diff_1_vs_10:.4f}")
print(f"  A control school 1 km from treatment has an absence rate")
print(f"  approximately {abs(diff_1_vs_10):.4f} {'lower' if diff_1_vs_10 > 0 else 'higher'} than one 10 km away.")

# (d) Compare to direct effect
print(f"\n(d) Comparison of Direct and Spillover Effects:")
print(f"  Direct treatment effect: {treatment_effect:.4f}")
print(f"  Spillover effect (1km vs 10km): {diff_1_vs_10:.4f}")
print(f"  Ratio (spillover/direct): {abs(diff_1_vs_10/treatment_effect):.2f}")

# ============================================================================
# QUESTION 7: Cost-Effectiveness
# ============================================================================
print("\n" + "="*80)
print("QUESTION 7: COST-EFFECTIVENESS")
print("="*80)

# (a) Costs
avg_cost = treat['program_cost_usd'].mean()
total_cost = treat['program_cost_usd'].sum()
n_treat = len(treat)

print(f"\n(a) Program Costs:")
print(f"  Average cost per treatment school: ${avg_cost:,.2f}")
print(f"  Total cost across {n_treat} treatment schools: ${total_cost:,.2f}")

# (b) Student-days gained at treatment schools
avg_n_students = df['num_students'].mean()
direct_effect = abs(treatment_effect)  # reduction in absence rate
student_days_per_treat_school = direct_effect * avg_n_students * 200  # 200 school days/year
total_student_days_treat = student_days_per_treat_school * n_treat

print(f"\n(b) Student-Days Gained at Treatment Schools:")
print(f"  Direct effect (reduction in absence rate): {direct_effect:.4f}")
print(f"  Average students per school: {avg_n_students:.1f}")
print(f"  School days per year: 200")
print(f"  Student-days gained per treatment school: {student_days_per_treat_school:,.1f}")
print(f"  Total student-days gained ({n_treat} treatment schools): {total_student_days_treat:,.1f}")

# (c) Spillover student-days at control schools
n_control = len(control_df)
avg_dist_control = control_df['distance_nearest_treat'].mean()

# Spillover benefit: compared to counterfactual of no treatment nearby
# Use the coefficient × average distance to estimate how much worse control schools
# would be if they were infinitely far from treatment (or use the predicted benefit
# at average distance relative to no spillover)
# Alternative: use avg distance to estimate spillover magnitude
spillover_at_avg_dist = coef_dist * avg_dist_control  # this is the "penalty" from being at avg distance vs 0
# The benefit from spillover = what a school gains from being at avg_dist vs very far
# Or more simply: use the regression to estimate spillover benefit per control school
# Spillover benefit for a control school at average distance vs. one at very large distance
# Let's use the approach suggested: compare to a hypothetical with no treatment nearby

# The spillover effect is: schools closer to treatment do better
# For an average control school, the spillover benefit relative to a school 
# with no nearby treatment can be approximated as coef_dist * (max_dist - avg_dist)
# But the simplest interpretation: the coefficient tells us the per-km effect
# A control school at distance d benefits by (coef_dist * d) less than one at distance 0

# More standard: avg spillover benefit per control school relative to no spillover
# = intercept effect. Use coefficient × average distance as the magnitude.
# Since positive coefficient means farther = higher absence, the spillover 
# benefit for a control school at average distance = the amount its absence rate 
# is lower than it would be if infinitely far away.

# Practical approach: use predicted absence at avg distance vs at some "no spillover" distance
# Let's assume schools beyond some distance have no spillover.
# Or simply use: spillover benefit per school ≈ |coef_dist| * avg_distance (from the wrong direction)

# Actually: Let's be more careful.
# If coef > 0: farther from treatment = higher absence. So closer = lower absence.
# Spillover benefit for avg control school = coef_dist * (reference_dist - avg_dist)
# If we assume spillover is 0 at some maximum distance, we can compute total gain.

# But the problem says "use the average distance" so:
# Average control school gains relative to having no treatment nearby.
# If coef_dist is positive, a school at distance d has absence rate that is 
# coef_dist * d higher than if d were 0 (i.e., right next to treatment).
# So the spillover benefit OF BEING AT average distance (vs. no treatment / very far) 
# is not well-defined without a reference point.

# Simpler approach as suggested by the problem: 
# Use the difference from part (c) of Q6: diff between 1km and 10km school
# Apply the coefficient to approximate average spillover benefit per control school
# relative to a no-spillover scenario.

# Let's use: avg spillover reduction in absence = |coef_dist * avg_dist| but signed appropriately
# The "benefit" a control school gets = how much lower its absence rate is compared to no spillover
# If coef > 0: being far away = worse. Being at distance d means absence is coef*d higher than at 0.
# So at avg distance, the school's absence is coef*avg_dist higher than if right next to treatment.
# But we want: how much better off is avg control school vs. no treatment anywhere?
# If no treatment anywhere is like d → ∞, that doesn't work with linear model.

# Best practical approach: use the 1km vs 10km difference from Q6(c)
# and note that the average control school gets some benefit.
# OR just use the coefficient directly:
# Spillover benefit per control school ≈ spillover effect × students × school days
# Where "spillover effect" = the predicted reduction from being at avg distance vs. far away

# Let me use a sensible reference: max distance observed
max_dist_control = control_df['distance_nearest_treat'].max()
spillover_per_school = coef_dist * (max_dist_control - avg_dist_control)
# This is how much HIGHER absence is at max_dist vs avg_dist
# So avg school benefits by this amount compared to the farthest school.

# Actually, let me just follow the problem's hint more directly:
# "You may use the average distance to nearest treatment school among controls"
# Spillover benefit at average distance vs. counterfactual of no treatment:
# The regression predicts that being 1km closer reduces absence by |coef_dist|
# A control school at average distance d_avg benefits by approximately:
# For this, let's say spillover benefit ≈ coef_dist × (some_max - avg_dist)
# OR the problem might just want us to use the 1km vs 10km difference somehow.

# Let me just use the simplest defensible approach:
# Average control school is at avg_dist. If it were at distance 0 (max spillover),
# its absence would be lower by coef_dist * avg_dist. 
# So the REDUCTION in benefit from not being at 0 = coef_dist * avg_dist
# The benefit it DOES receive (vs no spillover at all) is harder to pin down.

# Practical: just use the coefficient and the average distance.
# Spillover reduction in absence per control school = the difference between
# the predicted absence at the max distance and the predicted absence at avg distance

print(f"\n(c) Spillover Student-Days at Control Schools:")
print(f"  Number of control schools: {n_control}")
print(f"  Average distance to nearest treatment (controls): {avg_dist_control:.2f} km")
print(f"  Coefficient on distance: {coef_dist:.4f}")

# Approach: use the spillover coefficient. For each control school, its spillover
# benefit = coef * (max_possible_distance - actual_distance). 
# But without a clear "no spillover" reference, use the simpler approach:
# The average control school's absence rate is reduced by spillover.
# Approximate the average spillover effect as the coefficient × avg_distance
# relative to some baseline. 

# SIMPLEST: The problem asks to use the Q6 estimates.
# Spillover benefit for average control school vs. no-spillover counterfactual:
# If we assume the farthest control school gets ~zero spillover,
# then average spillover = coef * (max_dist - avg_dist)
# Or: use the simple regression coefficient and note that each km closer 
# gives |coef| reduction. Average school is avg_dist km away.
# Benefit relative to a school at, say, 10 km: coef * (10 - avg_dist)

# I'll compute it relative to a school at max observed distance (conservative):
spillover_benefit_per_school = abs(coef_dist) * (max_dist_control - avg_dist_control)
student_days_spill_per_school = spillover_benefit_per_school * avg_n_students * 200
total_student_days_spill = student_days_spill_per_school * n_control

print(f"  Max distance among controls: {max_dist_control:.2f} km")
print(f"  Spillover benefit per control school")
print(f"    (relative to farthest control school): {spillover_benefit_per_school:.4f}")
print(f"  Student-days gained per control school: {student_days_spill_per_school:,.1f}")
print(f"  Total spillover student-days ({n_control} control schools): {total_student_days_spill:,.1f}")

# Also compute using a different approach: average reduction relative to "no treatment"
# Assuming no spillover at distance = max
# Total spillover = sum over control schools of coef * (max_dist - school_dist)
total_spill_v2 = 0
for _, row in control_df.iterrows():
    benefit = abs(coef_dist) * (max_dist_control - row['distance_nearest_treat'])
    total_spill_v2 += benefit * row['num_students'] * 200
print(f"\n  Alternative (using each school's actual students and distance):")
print(f"  Total spillover student-days: {total_spill_v2:,.1f}")

# (d) Cost-effectiveness
print(f"\n(d) Cost-Effectiveness:")
# Direct only
cost_per_sd_direct = total_cost / total_student_days_treat
print(f"  (i) Direct effects only:")
print(f"      Total cost: ${total_cost:,.2f}")
print(f"      Total student-days gained: {total_student_days_treat:,.1f}")
print(f"      Cost per student-day: ${cost_per_sd_direct:.2f}")

# Including spillovers
total_student_days_all = total_student_days_treat + total_student_days_spill
cost_per_sd_all = total_cost / total_student_days_all
print(f"\n  (ii) Including spillover effects:")
print(f"      Total student-days gained: {total_student_days_all:,.1f}")
print(f"      Cost per student-day: ${cost_per_sd_all:.2f}")
print(f"\n  Improvement: Cost per student-day decreases by {((cost_per_sd_direct - cost_per_sd_all)/cost_per_sd_direct)*100:.1f}%")
print(f"  when accounting for spillovers.")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
