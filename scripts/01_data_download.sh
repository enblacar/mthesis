#  __           _       _      ___  _        ___      _              _                     _                 _ 
# / _\ ___ _ __(_)_ __ | |_   / _ \/ |_     /   \__ _| |_ __ _    __| | _____      ___ __ | | ___   __ _  __| |
# \ \ / __| '__| | '_ \| __| | | | | (_)   / /\ / _` | __/ _` |  / _` |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |
# _\ \ (__| |  | | |_) | |_  | |_| | |_   / /_// (_| | || (_| | | (_| | (_) \ V  V /| | | | | (_) | (_| | (_| |
# \__/\___|_|  |_| .__/ \__|  \___/|_(_) /___,' \__,_|\__\__,_|  \__,_|\___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|
#                |_|                                                                                           

# Script: 01_data_download.sh
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: Download NA12878 datasets from SRA using Aspera Connect.
DOWNLOADPATH="./data/raw_fq"
cd $DOWNLOADPATH

# Download each of the eight datafailes for NA12878.
nohup ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR091/ERR091571/ERR091571_1.fastq.gz . &
nohup ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR091/ERR091571/ERR091571_2.fastq.gz . &
nohup ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR091/ERR091572/ERR091572_1.fastq.gz . &
nohup ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR091/ERR091572/ERR091572_2.fastq.gz . &
nohup ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR091/ERR091573/ERR091573_1.fastq.gz . &
nohup ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR091/ERR091573/ERR091573_2.fastq.gz . &
nohup ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR091/ERR091574/ERR091574_1.fastq.gz . &
nohup ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR091/ERR091574/ERR091574_2.fastq.gz . &
wait