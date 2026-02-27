# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/hep/jtafoyav/vols/dependencies/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/hep/jtafoyav/vols/dependencies/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/hep/jtafoyav/vols/dependencies/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/hep/jtafoyav/vols/dependencies/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
echo "Conda environment variables stablished"

echo "Activating the icenet conda environment"
conda activate icenet
