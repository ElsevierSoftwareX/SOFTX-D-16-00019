echo ""
echo "The following command will save a postscript file named 6_channel.ps with a"
echo "density of 10 traces/cm in the OUTPUT folder"
pause "continue."
echo "Executing:"
echo "segy-change -f DATA/6_channel.seg -do_ps A4,10,0.001 > OUTPUT/6_channel.ps"
echo ""
segy-change -f DATA/6_channel.seg -do_ps A4,10,0.001 > OUTPUT/6_channel.ps
echo "Now only the second channel will be plotte ..."
pause "continue."
echo "Executing:"
echo "segy-change -f DATA/6_channel.seg -trace 2 2 -do_ps A4,10,0.001 > OUTPUT/6_channel_only_2.ps"
echo ""
segy-change -f DATA/6_channel.seg -trace 2 2 -do_ps A4,10,0.001 > OUTPUT/6_channel_only_2.ps

