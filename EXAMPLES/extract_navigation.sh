echo ""
echo "The following command extract coordinates from a segy file"
echo ""
echo "If the coordinates are in arc seconds they can be converted in degrees with"
echo "Executing:"
echo "segy-change -f DATA/6_channel.seg -dump_xy SOURCE | awk '{printf(\"%d %d %d %lf %lf\\n\", \$1, \$2, \$2, \$4/3600.00, \$5/3600.0);}'"
pause "continue."
segy-change -f DATA/6_channel.seg -dump_xy SOURCE | awk '{printf("%d %d %d %lf %lf\n", $1, $2, $2, $4/3600.00, $5/3600.0);}'


