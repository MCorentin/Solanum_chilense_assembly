# Add the barcode to the interleaved file
gunzip -c barcoded.fastq.gz | perl -ne 'chomp;$ct++;$ct=1 if($ct>4);if($ct==1){if(/(\@\S+)\sBX\:Z\:(\S{16})/){$flag=1;$head=$1."_".$2;print "$head\n";}else{$flag=0;}}else{print"$_\n" if($flag);}' > longranger_perl.fastq
