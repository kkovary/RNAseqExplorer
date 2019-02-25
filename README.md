# RNA-seq_Plot

This is a Shiny App that plots RNA-seq data

# To Update Shiny Server, type into terminal:
1) sudo rm -r /home/teruel/RNA-seq_Plot/
1) git clone https://github.com/kkovary/RNA-seq_Plot.git
2) sudo rm -r /srv/shiny-server/apps/RNA-seq_Plot
3) sudo cp -r /home/teruel/RNA-seq_Plot/ /srv/shiny-server/apps/RNA-seq_Plot/
4) sudo systemctl restart shiny-server

# Useful commands for updating Shiny Server
sudo nano RNAseq-shiny-20180309-152352-41371.log
sudo rm -r RNAseq
sudo cp -r Desktop/app_180403.R /srv/shiny-server/apps/RNAseq/app.R
sudo systemctl restart shiny-server
sudo su -     -c "R -e \"install.packages(c('ggplot2'), repos='http://cran.rstudio.com/')\""
