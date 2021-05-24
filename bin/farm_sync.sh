#!/usr/bin/env bash
# Script to sync between a local and remote directory(s)
#
# Setup for a given project folder using the config beow
#
# Will look for global and local .rsync_exclude/rsync_include files
# .rsync_exclude in $HOME is assumed to be global prefs
# Manual include overrides all excludes
#
# Sync only desired folders based on args

## Config ##
project_name="Sequence UNET Project"
local_dir=$HOME/phd/sequence_unet
remote_dirs=( ebi:/hps/research1/beltrao/ally/sequence_unet )
folders=( "data" "figures" "models" )

## Colours for printf ##
green=$(tput setaf 2)
magenta=$(tput setaf 5)
bold=$(tput bold)
normal=$(tput sgr0)

## Check for presence of include/exclude files ##
rsync_options=( -a -u -h )

if [ -e "$HOME/.rsync_exclude" ]; then
   rsync_options+=( --exclude-from "$HOME/.rsync_exclude" )
fi

if [ -e "$local_dir/.rsync_include" ]; then
   rsync_options+=( --include-from "$local_dir/.rsync_include" --exclude"="'*')
elif [ -e "$local_dir/.rsync_exclude" ]; then
   rsync_options+=( --exclude-from "$local_dir/.rsync_exclude" )
fi

## Sync function ##
syncr () {
   rsync -v --dry-run "${rsync_options[@]}" "$1" "$2"

   read -p "Transfer? " -n 1 -r
   echo
   if [[ $REPLY =~ ^[Yy]$ ]]
   then
      rsync "${rsync_options[@]}" "$1" "$2"
   fi
}

## Override folders if argument passed ##
if [ $# -ne 0 ]; then
   folders=( "$@" )
fi

## Perform sync ##
printf "%s\n" "${magenta}${bold}Rsyncing $project_name${normal}"
for r in "${remote_dirs[@]}"; do
   read -p "Sync to remote: $r? " -n 1 -r
   echo
   if [[ $REPLY =~ ^[Yy]$ ]]; then
      for f in "${folders[@]}"; do
         printf "\n%s\n%s\n" "${green}${bold}Folder: $f${normal}" "${green}Local -> Remote${normal}"
         syncr "$local_dir/$f/" "$r/$f"
         printf "\n%s\n" "${green}Local <- Remote${normal}"
         syncr "$r/$f/" "$local_dir/$f"
      done
   fi
done


