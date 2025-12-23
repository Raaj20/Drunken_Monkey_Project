# Drunken_Monkey_Project

# Install Homebrew 
% git clone https://github.com/Homebrew/brew homebrew
% eval "$(homebrew/bin/brew shellenv)"
brew update --force --quiet
chmod -R go-w "$(brew --prefix)/share/zsh"

# Create an environment in Conda
% brew install miniforge
% conda update -n base -c conda-forge conda

% conda --version
conda 25.11.1

conda create -n mirna_pipeline python=3.10

% echo $SHELL
/bin/zsh

conda init zsh
% conda activate mirna_pipeline

# Restart Terminal

% conda activate mirna_pipeline

