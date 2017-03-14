echo ""
echo "The file DATA/6_channel.seg has incorrect navigation, with a navigation file"
echo "and segy-change it is possible to update coordinates with correct values."
echo "The following commands update coordinates of a segy file, given a nav file."
echo "containing for each line:"
echo "record_nr lon W|E lat N|S heading sp date time"
echo ""
echo "Executing:"
echo ""
echo "cat DATA/6_channel.nav | awk '{ lon=\$2; 
       if(\$3 == \"W\") lon=-\$2; 
       lat=\$4; 
       if(\$5 == \"S\") lat=-\$4; 
       for(i=1;i<=6; i++) printf(\"%d 0 %d %lf %lf arcsec\\n\", \$1, i, lon * 3600, lat * 3600); 
     }' > OUTPUT/coords.txt 

segy-change -f DATA/6_channel.seg -add_xy OUTPUT/coords.txt,SOURCE -o OUTPUT/6_channel_with_navigation.seg
"
pause "continue."
cat DATA/6_channel.nav | awk '{ lon=$2; 
       if($3 == "W") lon=-$2; 
       lat=$4; 
       if($5 == "S") lat=-$4; 
       for(i=1;i<=6; i++) printf("%d 0 %d %lf %lf arcsec\n", $1, i, lon * 3600, lat * 3600); 
     }' > OUTPUT/coords.txt 
segy-change -f DATA/6_channel.seg -add_xy OUTPUT/coords.txt,SOURCE -o OUTPUT/6_channel_with_navigation.seg
echo "Please check that OUTPUT/6_channel_with_navigation.seg has been updated."
