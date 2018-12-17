
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod ugo+x Miniconda3-latest-Linux-x86_64.sh 
bash ./Miniconda3-latest-Linux-x86_64.sh -b
source ~/.bashrc
export PATH=~/miniconda3/bin:$PATH 
conda update conda -y
echo " " >> ~/.bashrc
echo "PATH=~/miniconda3/bin:$PATH"  >> ~/.bashrc

