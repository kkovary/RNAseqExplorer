# RNAseqExplorer

This is a Shiny App that plots RNA-seq data
http://kylekovary.shinyapps.io/RNAseqExplorer/

# Features
- Plot time courses of multiple genes of transcripts per million
- Volcano plots
- Output plots and data files

# To Update Shiny Server, type into terminal:
1) sudo rm -r /srv/shiny-server/apps/RNAseqExplorer
2) sudo git clone https://github.com/kkovary/RNAseqExplorer.git /srv/shiny-server/apps/
3) sudo systemctl restart shiny-server

# Useful info for Shiny Server
- Log files are located at /var/log/shiny-server/
- List log files by newest: ls -lt
- Open log file: sudo nano

# Features to add
- Transcript variant plots
- Add data points to plot
- Sliders for error bar and point sizes
- Fold change plots (over D0 or siRNA)
- Separate program idea: Shiny interface for setting up Kallisto analysis, run local or on server
