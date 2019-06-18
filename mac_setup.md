## Setting up a new macbook pro

1. Download and install [iTerm](https://iterm2.com)

2. Copy .bashrc file to $HOME dir (from previous computer)

3. Create a .bash_profile file and add the following script to initiate .bashrc at login.
```
if [ -f $HOME/.bashrc ]; then
        source $HOME/.bashrc
fi
```

4. Install [homebrew](https://brew.sh)

5. Install homebrew version of [Python as as well as virtualenv packages](https://swapps.com/blog/how-to-configure-virtualenvwrapper-with-python3-in-osx-mojave/)

6. brew cask install atom
