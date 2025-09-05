# =============================================================================
# STEP 1: SETUP AND IMPORTS
# =============================================================================
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

print("Setup and Imports Complete.")

# =============================================================================
# STEP 2: DEFINE THE SCORING FUNCTION
# =============================================================================

def calculate_boom_bust_score(df, position='WR'):
    """
    Analyzes a dataframe of player stats to identify boom-bust players.

    Args:
        df (pd.DataFrame): The comprehensive stats dataframe for a single position.
        position (str): The player position ('QB', 'RB', 'WR', 'TE') to adjust thresholds.

    Returns:
        pd.DataFrame: A summary dataframe with boom-bust metrics for each player.
    """
    print(f"\n--- Calculating Boom-Bust metrics for {position}s ---")

    # Define fantasy point thresholds for 'boom' and 'bust' weeks by position
    if position == 'QB':
        boom_threshold = 28
        bust_threshold = 15
    elif position == 'RB':
        boom_threshold = 22
        bust_threshold = 10
    elif position == 'WR':
        boom_threshold = 20
        bust_threshold = 8
    else: # TE
        boom_threshold = 15
        bust_threshold = 6

    # Group by player to calculate career stats
    # We use player_display_name as the primary key for readability
    player_summary = df.groupby('player_display_name').agg(
        games_played=('week', 'count'),
        avg_dk_points=('dk_fantasy_points', 'mean'),
        std_dev_dk_points=('dk_fantasy_points', 'std'),
        # Add aDOT and Air Yards Share for WRs/TEs if available
        avg_adot=('average_depth_of_target' if 'average_depth_of_target' in df.columns else 'dk_fantasy_points', 'mean'),
        avg_air_yards_share=('air_yards_share' if 'air_yards_share' in df.columns else 'dk_fantasy_points', 'mean')
    ).reset_index()

    # Filter for players with a meaningful sample size of games
    min_games = 16 # Roughly one full season's worth of games
    player_summary = player_summary[player_summary['games_played'] >= min_games]

    # Calculate Boom % and Bust %
    boom_weeks = df[df['dk_fantasy_points'] >= boom_threshold].groupby('player_display_name').size()
    bust_weeks = df[df['dk_fantasy_points'] <= bust_threshold].groupby('player_display_name').size()

    player_summary = pd.merge(player_summary, boom_weeks.rename('boom_weeks'), on='player_display_name', how='left')
    player_summary = pd.merge(player_summary, bust_weeks.rename('bust_weeks'), on='player_display_name', how='left')
    player_summary.fillna(0, inplace=True)

    player_summary['boom_pct'] = player_summary['boom_weeks'] / player_summary['games_played']
    player_summary['bust_pct'] = player_summary['bust_weeks'] / player_summary['games_played']

    # --- Calculate Key Boom-Bust Indicators ---
    # Coefficient of Variation (CV): a normalized measure of volatility
    player_summary['coeff_of_variation'] = player_summary['std_dev_dk_points'] / player_summary['avg_dk_points']
    
    # --- Create the Final Boom-Bust Score ---
    # We will scale each component metric from 0 to 1 and then combine them.
    scaler = MinMaxScaler()
    
    # Define the components of the score. Higher values for these metrics suggest a boom-bust profile.
    score_components = ['coeff_of_variation', 'boom_pct', 'bust_pct']
    if position in ['WR', 'TE']:
         # aDOT and Air Yards Share are leading indicators for WRs/TEs
        score_components.extend(['avg_adot', 'avg_air_yards_share'])

    # Scale the components
    player_summary_scaled = player_summary.copy()
    player_summary_scaled[score_components] = scaler.fit_transform(player_summary_scaled[score_components])
    
    # Calculate the final score by averaging the scaled components
    player_summary['boom_bust_score'] = player_summary_scaled[score_components].mean(axis=1)

    # Sort by the final score
    final_report = player_summary.sort_values(by='boom_bust_score', ascending=False)
    
    return final_report

# =============================================================================
# STEP 3: PROCESS EACH POSITION
# =============================================================================

# Define file paths for the comprehensive CSVs you created
files = {
    'QB': 'qb_comprehensive_advanced_stats_2020_2024.csv',
    'RB': 'rb_comprehensive_advanced_stats_2020_2024.csv',
    'WR': 'wr_comprehensive_advanced_stats_2020_2024.csv',
    'TE': 'te_comprehensive_advanced_stats_2020_2024.csv'
}

for position, file_path in files.items():
    try:
        # Load the data
        df = pd.read_csv(file_path)
        
        # Calculate the boom-bust report
        report = calculate_boom_bust_score(df, position=position)
        
        # Select and rename columns for a clean final output
        output_cols = [
            'player_display_name', 'games_played', 'avg_dk_points', 
            'std_dev_dk_points', 'coeff_of_variation', 
            'boom_pct', 'bust_pct', 'boom_bust_score'
        ]
        if 'avg_adot' in report.columns:
            output_cols.insert(5, 'avg_adot')
        
        # Ensure all columns exist before selecting
        report = report[[col for col in output_cols if col in report.columns]]
        
        # Define output filename
        output_filename = f'{position}_boom_bust_report.csv'
        
        # Export to CSV
        report.to_csv(output_filename, index=False, float_format='%.3f')
        
        print(f"✅ Success! Report for {position}s exported to '{output_filename}'")
        print("\n--- Top 5 Boom-Bust Players ---")
        print(report.head())
        
    except FileNotFoundError:
        print(f"\nWarning: '{file_path}' not found. Skipping {position}s.")
    except Exception as e:
        print(f"\nAn error occurred while processing {position}s: {e}")