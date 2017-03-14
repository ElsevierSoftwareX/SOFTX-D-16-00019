echo ""
echo "The following command convert the endianess of a segy file, converting values"
echo "stored inside the headers as well as traces values."
echo ""
pause "continue."
echo "Executing:"
echo "segy-change -f DATA/wrong_endianess.seg -flip_endianess -o DATA/adjusted_endianess.seg"
segy-change -f DATA/wrong_endianess.seg -flip_endianess -o OUTPUT/adjusted_endianess.seg
echo ""
echo "Now a postscript file named wrong_endianess.ps  with a density of 10 traces/cm"
echo "will be saved in the OUTPUT folder"
pause "continue."
echo ""
echo "Executing:"
echo "segy-change -f DATA/wrong_endianess.seg -flip_endianess -do_ps A4,10,0.001 > wrong_endianess.ps"
segy-change -f DATA/wrong_endianess.seg -flip_endianess -do_ps A4,10,0.001 > OUTPUT/wrong_endianess.ps

