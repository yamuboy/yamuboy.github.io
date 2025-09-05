import pandas as pd
import numpy as np
import nfl_data_py as nfl

# Load 2023 play-by-play data, which is the foundation for most advanced stats
pbp_df = nfl.import_pbp_data([2024])

# Filter for regular season plays only
pbp_df = pbp_df[pbp_df['season_type'] == 'REG']

# We'll also need roster data to link IDs to names and positions
roster_df = nfl.import_weekly_rosters([2023])

# --- Calculate Sack Yards Manually (THE FIX) ---
# Filter for sack plays
sack_plays = pbp_df[pbp_df['sack'] == 1]
# Group by the player who was sacked (the passer) and sum the yards lost
sack_yards_per_player = sack_plays.groupby('passer_player_id', as_index=False).agg(
    sack_yards=('yards_gained', 'sum')
)
# Note: 'sack_yards' will be a negative number. For ANY/A, we need to subtract it, so we'll flip the sign.
sack_yards_per_player['sack_yards'] = -sack_yards_per_player['sack_yards']


# --- Aggregate other QB stats from play-by-play data ---
qb_stats = pbp_df.groupby('passer_player_id', as_index=False).agg(
    player_name=('passer_player_name', 'first'),
    attempts=('play_id', lambda x: x[pbp_df.loc[x.index, 'pass_attempt'] == 1].count()),
    completions=('play_id', lambda x: x[pbp_df.loc[x.index, 'complete_pass'] == 1].count()),
    passing_yards=('passing_yards', 'sum'),
    passing_tds=('pass_touchdown', 'sum'),
    interceptions=('interception', 'sum'),
    sacks=('sack', 'sum')
)

# Merge our manually calculated sack_yards into the main qb_stats DataFrame
qb_stats = pd.merge(qb_stats, sack_yards_per_player, on='passer_player_id', how='left')
# Fill any NaNs with 0 for players who were not sacked
qb_stats['sack_yards'].fillna(0, inplace=True)

# Filter for QBs with a decent number of attempts
qb_stats = qb_stats[qb_stats['attempts'] > 100]

# Calculate ANY/A (this formula now works)
# We subtract sack_yards, which is now a positive value representing yards lost.
qb_stats['any_a'] = (
    (qb_stats['passing_yards'] + (20 * qb_stats['passing_tds']) - (45 * qb_stats['interceptions']) - qb_stats['sack_yards']) /
    (qb_stats['attempts'] + qb_stats['sacks'])
)

print("\n--- Top 5 QBs by ANY/A (2023) ---")
print(qb_stats.sort_values(by='any_a', ascending=False).head()[['player_name', 'any_a']])

# --- CPOE Calculation (Unchanged) ---
pass_attempts_df = pbp_df.dropna(subset=['passer_player_id', 'cp'])
cpoe_stats = pass_attempts_df.groupby('passer_player_id', as_index=False).agg(
    player_name=('passer_player_name', 'first'),
    attempts=('play_id', 'count'),
    cpoe=('cpoe', 'mean')
)
cpoe_stats = cpoe_stats[cpoe_stats['attempts'] > 100]

print("\n--- Top 5 QBs by CPOE (2023) ---")
print(cpoe_stats.sort_values(by='cpoe', ascending=False).head()[['player_name', 'cpoe']])

##############   WR ########


# Filter for plays where a player was targeted
targets_df = pbp_df.dropna(subset=['receiver_player_id', 'posteam'])

# Count total targets for each team
team_targets = targets_df.groupby('posteam')['play_id'].count().reset_index().rename(columns={'play_id': 'team_targets'})

# Count targets for each individual player
player_targets = targets_df.groupby(['receiver_player_id', 'posteam'], as_index=False).agg(
    player_name=('receiver_player_name', 'first'),
    targets=('play_id', 'count')
)

# Merge and calculate target share
target_share_df = pd.merge(player_targets, team_targets, on='posteam')
target_share_df['target_share'] = target_share_df['targets'] / target_share_df['team_targets']

print("\n--- Top 5 WRs/TEs by Target Share (2023) ---")
print(target_share_df.sort_values(by='target_share', ascending=False).head()[['player_name', 'target_share']])

# Filter for pass attempts with air yards data
air_yards_df = pbp_df.dropna(subset=['receiver_player_id', 'posteam', 'air_yards'])

# Sum total air yards for each team
team_air_yards = air_yards_df.groupby('posteam')['air_yards'].sum().reset_index().rename(columns={'air_yards': 'team_air_yards'})

# Sum air yards for each player
player_air_yards = air_yards_df.groupby(['receiver_player_id', 'posteam'], as_index=False).agg(
    player_name=('receiver_player_name', 'first'),
    air_yards=('air_yards', 'sum')
)

# Merge and calculate air yards share
air_yards_share_df = pd.merge(player_air_yards, team_air_yards, on='posteam')
air_yards_share_df['air_yards_share'] = air_yards_share_df['air_yards'] / air_yards_share_df['team_air_yards']

print("\n--- Top 5 WRs/TEs by Air Yards Share (2023) ---")
print(air_yards_share_df.sort_values(by='air_yards_share', ascending=False).head()[['player_name', 'air_yards_share']])

# Filter for plays with yards after catch data
yac_df = pbp_df.dropna(subset=['receiver_player_id', 'posteam', 'yards_after_catch'])
# Sum total yards after catch for each team
team_yac = yac_df.groupby('posteam')['yards_after_catch'].sum().reset_index().rename(columns={'yards_after_catch': 'team_yac'})
# Sum yards after catch for each player
player_yac = yac_df.groupby(['receiver_player_id', 'posteam'], as_index=False).agg(
    player_name=('receiver_player_name', 'first'),
    yac=('yards_after_catch', 'sum')
)
# Merge and calculate YAC share
yac_share_df = pd.merge(player_yac, team_yac, on='posteam')
yac_share_df['yac_share'] = yac_share_df['yac'] / yac_share_df['team_yac']
print("\n--- Top 5 WRs/TEs by YAC Share (2023) ---")
print(yac_share_df.sort_values(by='yac_share', ascending=False).head()[['player_name', 'yac_share']])
# =============================================================================

# =============================================================================
# 3. Running Backs (RBs) - CORRECTED SECTION
# =============================================================================
print("\n--- Calculating Advanced Metrics for RBs using Available Data ---")
# Filter for actual rush attempts
rb_df = pbp_df[pbp_df['rush_attempt'] == 1].dropna(subset=['rusher_player_id', 'rushing_yards'])

# Define a breakaway run (e.g., 15+ yards)
rb_df['breakaway_run'] = (rb_df['rushing_yards'] >= 15).astype(int)

# Aggregate stats for each rusher
rb_stats = rb_df.groupby('rusher_player_id', as_index=False).agg(
    player_name=('rusher_player_name', 'first'),
    carries=('play_id', 'count'),
    total_rushing_yards=('rushing_yards', 'sum'),
    breakaway_runs=('breakaway_run', 'sum')
)

# Filter for RBs with a meaningful number of carries
rb_stats = rb_stats[rb_stats['carries'] > 50].dropna(subset=['player_name'])

# Calculate Yards Per Carry (YPC)
rb_stats['ypc'] = rb_stats['total_rushing_yards'] / rb_stats['carries']

# Calculate Breakaway Run Rate
rb_stats['breakaway_rate'] = rb_stats['breakaway_runs'] / rb_stats['carries']

print("\n--- Top 5 RBs by Yards Per Carry (min 50 carries, 2023) ---")
print(rb_stats.sort_values(by='ypc', ascending=False).head()[['player_name', 'carries', 'ypc']])

print("\n--- Top 5 RBs by Breakaway Run Rate (min 50 carries, 2023) ---")
print(rb_stats.sort_values(by='breakaway_rate', ascending=False).head()[['player_name', 'carries', 'breakaway_rate']])

# =============================================================================
# 4. Defenses (DST) - (Verified to work)
# =============================================================================
dropbacks_df = pbp_df[(pbp_df['pass_attempt'] == 1) | (pbp_df['sack'] == 1)]
defense_pressure_stats = dropbacks_df.groupby('defteam', as_index=False).agg(
    total_dropbacks=('play_id', 'count'),
    sacks=('sack', 'sum'),
    hits_and_hurries=('qb_hit', 'sum')
)
defense_pressure_stats['sack_rate'] = defense_pressure_stats['sacks'] / defense_pressure_stats['total_dropbacks']
defense_pressure_stats['pressure_rate'] = (defense_pressure_stats['sacks'] + defense_pressure_stats['hits_and_hurries']) / defense_pressure_stats['total_dropbacks']

print("\n--- Top 5 Defenses by Sack Rate (2023) ---")
print(defense_pressure_stats.sort_values(by='sack_rate', ascending=False).head()[['defteam', 'sack_rate']])