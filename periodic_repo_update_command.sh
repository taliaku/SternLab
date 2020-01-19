#!/usr/bin/env bash
### file to be used by crontab ###

# in crontab configuration add the following line:
# MAILTO=""
# * * * * * /sternadi/home/volume1/shared/SternLab/periodic_repo_update_command.sh

cd /sternadi/home/volume1/shared/SternLab/
pwd
git pull

# this is an update!!