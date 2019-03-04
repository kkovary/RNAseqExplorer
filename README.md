# RNAseqExplorer

This is a Shiny App that plots RNA-seq data

# Features
- Plot time courses of multiple genes of transcripts per million
- Volcano plots
- Output plots and data files

# To Update Shiny Server, type into terminal:
1) sudo rm -r /srv/shiny-server/apps/RNAseq
3) sudo git clone https://github.com/kkovary/RNA-seq_Plot.git /srv/shiny-server/apps/RNAseq
4) sudo systemctl restart shiny-server

# To Setup Repository on Server
sudo git clone https://github.com/kkovary/RNA-seq_Plot.git /srv/shiny-server/apps/RNAseq

# Useful commands for updating Shiny Server
- sudo nano RNAseq-shiny-20180309-152352-41371.log
- sudo rm -r RNAseq
- sudo cp -r Desktop/app_180403.R /srv/shiny-server/apps/RNAseq/app.R
- sudo systemctl restart shiny-server
- sudo su - -c "R -q -e \"install.packages('mypackage', repos='http://cran.rstudio.com/')\""

# Features to add
- Transcript variant plots
- Add data points to plot
- Sliders for error bar and point sizes
- Fold change plots (over D0 or siRNA)
