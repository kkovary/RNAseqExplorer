# RNA-seq_Plot

This is a Shiny App that plots RNA-seq data


# Useful commands for updating Shiny Server
sudo nano RNAseq-shiny-20180309-152352-41371.log
sudo rm -r RNAseq
sudo cp -r Desktop/app_180403.R /srv/shiny-server/apps/RNAseq/app.R
sudo systemctl restart shiny-server
sudo su -     -c "R -e \"install.packages(c('ggplot2'), repos='http://cran.rstudio.com/')\""