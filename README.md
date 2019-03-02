# RNA-seq_Plot

This is a Shiny App that plots RNA-seq data

# To Update Shiny Server, type into terminal:

1) cd /srv/shiny-server/apps/RNAseq/
2) sudo git pull
3) sudo systemctl restart shiny-server

# To Setup Repository on Server
sudo git clone https://github.com/kkovary/RNA-seq_Plot.git /srv/shiny-server/apps/RNAseq

# Useful commands for updating Shiny Server
sudo nano RNAseq-shiny-20180309-152352-41371.log
sudo rm -r RNAseq
sudo cp -r Desktop/app_180403.R /srv/shiny-server/apps/RNAseq/app.R
sudo systemctl restart shiny-server
sudo su -     -c "R -e \"install.packages(c('ggplot2'), repos='http://cran.rstudio.com/')\""
