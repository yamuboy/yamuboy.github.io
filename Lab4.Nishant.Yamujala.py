import rdata as r
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np
from scipy import stats

def r_summary(results, formula=None, data_name=None):
    """
    Generates and prints a regression summary that mimics the output of R's summary(lm()).

    Args:
        results (statsmodels.regression.linear_model.RegressionResultsWrapper):
            A fitted OLS results object from statsmodels.
        formula (str, optional):
            The formula string used to fit the model (e.g., "y ~ x").
        data_name (str, optional):
            The name of the DataFrame used.
    """

    # --- Helper function for p-value formatting ---
    def format_p_value(p):
        if p < 2e-16:
            return "<2e-16"
        elif p < 1e-4:
            return f"{p:.2e}"
        else:
            return f"{p:.4f}"

    # --- Helper function for significance codes ---
    def get_signif_codes(p_values):
        codes = []
        for p in p_values:
            if p < 0.001:
                codes.append('***')
            elif p < 0.01:
                codes.append('**')
            elif p < 0.05:
                codes.append('*')
            elif p < 0.1:
                codes.append('.')
            else:
                codes.append(' ')
        return codes

    print("##")
    print("## Call:")
    # Use formula from results object if available, otherwise use provided argument
    if hasattr(results.model, 'formula') and formula is None:
        formula = results.model.formula
    
    if formula and data_name:
        print(f"## lm(formula = {formula}, data = {data_name})")
    elif formula:
        print(f"## lm(formula = {formula})")
    else:
        print("## lm(formula = unknown)")
    print("##")

    # --- Residuals Section ---
    residuals = results.resid
    resid_quantiles = np.quantile(residuals, [0, 0.25, 0.5, 0.75, 1])
    resid_labels = ["   Min  ", "    1Q  ", " Median ", "    3Q  ", "   Max   "]
    print("## Residuals:")
    print("## " + "".join([f"{label}{q: >7.4f}" for label, q in zip(resid_labels, resid_quantiles)]))
    print("##")

    # --- Coefficients Section ---
    coeffs = pd.DataFrame({
        'Estimate': results.params,
        'Std. Error': results.bse,
        't value': results.tvalues,
        'Pr(>|t|)': results.pvalues
    })
    coeffs['Signif'] = get_signif_codes(coeffs['Pr(>|t|)'])
    
    # Formatting for printing
    max_name_len = max(len(name) for name in coeffs.index)
    
    print("## Coefficients:")
    header = f"## {' ':<{max_name_len}} {'Estimate':>10} {'Std. Error':>11} {'t value':>9} {'Pr(>|t|)':>10}"
    print(header)
    
    for i, row in coeffs.iterrows():
        p_val_formatted = format_p_value(row['Pr(>|t|)'])
        line = (f"## {i:<{max_name_len}} {row['Estimate']:>10.4f} {row['Std. Error']:>11.4f} "
                f"{row['t value']:>9.2f} {p_val_formatted:>10} {row['Signif']}")
        print(line)

    print("## ---")
    print("## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    print("##")

    # --- Model Fit Section ---
    rse = np.sqrt(results.scale)
    print(f"## Residual standard error: {rse:.3f} on {int(results.df_resid)} degrees of freedom")
    
    r_squared_line = f"## Multiple R-squared:  {results.rsquared:.5f},"
    adj_r_squared_line = f"Adjusted R-squared:  {results.rsquared_adj:.5f}"
    print(f"{r_squared_line:<43} {adj_r_squared_line}")
    
    f_stat = f"F-statistic: {results.fvalue:.2f} on {int(results.df_model)} and {int(results.df_resid)} DF,"
    f_p_val = f"p-value: {format_p_value(results.f_pvalue)}"
    print(f"## {f_stat:<48} {f_p_val}")
    print("##")

if __name__ == "__main__":

    converted_data = r.read_rda('Ch8_Exercise4_HOPE_scholarship.RData')
    # Set the number of significant digits to 3 for better readability
    pd.set_option('display.float_format', '{:.3f}'.format)
    # Load the dataset
    df = pd.DataFrame(converted_data['dta'])
    #print(df.head())
    #print(df.describe())
    #print(df.columns)
    print("Dataset loaded successfully.")

    # run a difference in difference model
    did_model = smf.ols('InCollege ~ AfterGeorgia + Georgia + After', data=df).fit()
    print("DiD Model Results:")
    r_summary(did_model, formula="FiredCoach ~ WinPct", data_name="df")

    results = pd.DataFrame({
        "dep_vars": did_model.params
    })
    print("(Intercept) \n {:0.3f}".format(results["dep_vars"].iloc[0]))
    print("(Intercept) \n {:0.3f}".format(results["dep_vars"].iloc[0] + results["dep_vars"].iloc[2]))
    print("(Intercept) \n {:0.3f}".format(results["dep_vars"].iloc[0] + results["dep_vars"].iloc[3]))
    print("(Intercept) \n {:0.3f}".format(results["dep_vars"].sum()))

    print("## {:0.3f}".format(df.loc[(df['After'] == 0) & (df['Georgia'] == 0), 'InCollege'].mean()))
    print("## {:0.3f}".format(df.loc[(df['After'] == 0) & (df['Georgia'] == 1), 'InCollege'].mean()))
    print("## {:0.3f}".format(df.loc[(df['After'] == 1) & (df['Georgia'] == 0), 'InCollege'].mean()))
    print("## {:0.3f}".format(df.loc[(df['After'] == 1) & (df['Georgia'] == 1), 'InCollege'].mean()))


    import seaborn as sns
    fig, ax = plt.subplots(figsize=(8, 6))

    sns.regplot(
        x='After',
        y='InCollege',
        data=df[df['Georgia'] == 0], # Filter data for the control group
        color='blue',
        ax=ax,
        scatter=False, # Don't show individual data points
        ci=None,       # Don't show the confidence interval
        line_kws={'linewidth': 2},
        label="Georgia's neighbors" # Label for the legend
    )

    sns.regplot(
        x='After',
        y='InCollege',
        data=df[df['Georgia'] == 1], # Filter data for the treatment group
        color='red',
        ax=ax,
        scatter=False,
        ci=None,
        line_kws={'linewidth': 2},
        label="Georgia" # Label for the legend
    )

    ax.set_xlabel("Time Period")
    ax.set_ylabel("College Enrollment Rate")
    ax.set_ylim(0, 0.5)
    ax.set_title("College Enrollment Trends: Georgia vs. Neighbors")

    # Set cleaner x-axis ticks and labels
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Before 1993', 'After 1993'])

    # Create the legend, equivalent to R's legend()
    ax.legend(loc="upper left") # "upper left" is often clearer for DiD plots

    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()

    from linearmodels.panel.model import PanelOLS

    df_indexed = df.set_index(['StateCode', 'Year'])

    model_formula = 'InCollege ~ 1 + AfterGeorgia + EntityEffects + TimeEffects'
    result_linearmodels = PanelOLS.from_formula(model_formula, data=df_indexed).fit()

    print("--- Result from linearmodels (Recommended) ---")
    print(result_linearmodels)