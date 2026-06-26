#!/bin/bash
# [[_hpc/notify.sbatch]]
# Goal: Send a Discord notification when a monitored Slurm job completes.

# Steps:
#   1. Parse target job ID to monitor from command line argument or environment variable.
#   2. Query final execution state of the target job using sacct.
#   3. Send a Discord webhook notification with the job details and status.

#SBATCH --job-name=job_notify
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=_data/logs/notify/notify_%j.log
#SBATCH --error=_data/logs/notify/notify_%j.err

# Exit on error
set -e

# --- Input Parsing & Assertions ---
MONITOR_JOB="${1:-$MONITOR_JOB}"
if [ -z "$MONITOR_JOB" ]; then
    echo "Assertion Error: No job ID provided to monitor." >&2
    echo "Usage: sbatch [sbatch_options] _hpc/notify.sbatch <JOB_ID>" >&2
    exit 1
fi

# Assert that MONITOR_JOB is numeric
if [[ ! "$MONITOR_JOB" =~ ^[0-9]+$ ]]; then
    echo "Assertion Error: Monitored job ID must be numeric. Got: '$MONITOR_JOB'" >&2
    exit 1
fi

DISCORD_WEBHOOK="https://discord.com/api/webhooks/1518762676963119175/xzn_7Fvdsz387HOOcMBPdqSkQFjsRhevAlaURg8D1RPRfbiuZXW-bCew_4eEESLOdtBf"

# Query the job state from sacct
exit_status=$(sacct -j "$MONITOR_JOB" --format=state -n | head -n 1 | xargs)

echo "========================================================="
echo "   Slurm Job Notification Sender"
echo "   Monitoring Job: $MONITOR_JOB"
echo "   Status Resolved: $exit_status"
echo "========================================================="

# Send the post request to Discord
curl -s -H "Content-Type: application/json" -X POST \
  -d "{\"content\": \"🔔 **[HPC Notification]** Slurm Job **$MONITOR_JOB** has finished.\\n**Final Status**: ${exit_status:-UNKNOWN}\\n**Log**: \`_data/logs/notify/download_${MONITOR_JOB}.log\`\"}" \
  "$DISCORD_WEBHOOK" > /dev/null || true

echo "Notification sent."
